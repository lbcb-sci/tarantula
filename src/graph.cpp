// Copyright (c) 2021 Cecilia Lee, Robert Vaser

#include "graph.hpp"

#include <fstream>
#include <iostream>

#include "biosoup/timer.hpp"
#include "cereal/archives/json.hpp"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace tarantula {

Graph::Graph(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool)
    : thread_pool_(thread_pool ?
        thread_pool :
        std::make_shared<thread_pool::ThreadPool>(1)) {
}

void Graph::Construct(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>& targets,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences) {
  // paramters for RAM
  uint32_t k = 21, w = 11, bandwidth = 100, chain = 2, matches = 25, gap = 100;
  double frequency = 0.01;
  int num_pair = 0;

  // sort the sequence first
  std::cerr << "[tarantula::Construct] Sorting: "
            << sequences.size()
            << std::endl;

  std::sort(sequences.begin(), sequences.end(),
    [] (const std::unique_ptr<biosoup::NucleicAcid>& lhs,
        const std::unique_ptr<biosoup::NucleicAcid>& rhs) -> bool {
      return lhs->name < rhs->name;
    });

  // map hic_reads to assembly
  ram::MinimizerEngine minimizer_engine{
    thread_pool_,
    k,
    w,
    bandwidth,
    chain,
    matches,
    gap};

  minimizer_engine.Minimize(targets.begin(), targets.end());
  minimizer_engine.Filter(frequency);

  if (sequences.empty()) {
    return;
  }

  std::vector<std::future<std::pair<std::string, std::vector<std::vector<biosoup::Overlap>>>>> futures;
  std::unordered_map<std::string, std::vector<std::vector<biosoup::Overlap>>> multiple_overlap_read_pairs;
  std::unordered_map<std::string, std::vector<std::vector<biosoup::Overlap>>> read_pairs;
  std::unordered_map<std::string, std::vector<std::vector<biosoup::Overlap>>> interchromosome_read_pairs;
  biosoup::Timer timer;

  timer.Start();
  // then create pilo-o-gram per contig
  CreateGraph(targets);
  std::cerr << "[tarantula::Construct] Graph created, number of nodes: "
            << contigs.size() << " "
            << timer.Stop() << "s"
            << std::endl;
  // filter pair + create pile-o-gram + less than 4GB
  timer.Start();
  std::size_t bytes = 0;
  for (std::uint32_t i=0; i< sequences.size()-1; i++) {
    if (sequences[i]->name.compare(sequences[i+1]->name) == 0) {
      num_pair+=1;
      Graph::Process(futures, minimizer_engine, sequences[i], sequences[i+1]);
      bytes += sequences[i]->inflated_len + sequences[i+1]->inflated_len;
      i++;
    } else {
      continue;
    }
    // less than 4GB
    if (bytes >= (1ULL << 32) || i >= sequences.size()-1) {
      // wait for futures
      for (auto& it : futures) {
        // if it == empty , discard
        auto result = it.get();
        if (result.first.find("_interchromosome") != std::string::npos) {
          interchromosome_read_pairs.insert(result);
        } else if (result.first.find("_multiple") != std::string::npos) {
          multiple_overlap_read_pairs.insert(result);
        } else if (result.first.compare("empty") != 0) {
          read_pairs.insert(result);
        }
      }
      std::cerr << "[tarantula::Construct] Number of good read pair: "
          << read_pairs.size() << " "
          << timer.Stop() << "s"
          << std::endl;

      // fill pile-o-gram
      timer.Start();
      FillPileogram(read_pairs);
      std::cerr << "[tarantula::Construct] Pile-o-gram created, number of nodes: "
                << contigs.size() << " "
                << timer.Stop() << "s"
                << std::endl;

      // discard read pair
      read_pairs.clear();
      futures.clear();
      bytes = 0;
    }
  }
  CalcualteInterChromosomeLinks(interchromosome_read_pairs);
  uint32_t sum_interchromosome_links = 0, sum_intrachromosome_links = 0;
  int num_chromosomes = targets.size();
  std::vector<std::vector<std::uint32_t>> matrix(num_chromosomes, std::vector<std::uint32_t>(num_chromosomes));
  std::unordered_map<std::uint32_t, Node>::iterator it;
  GenerateMatrix(matrix, interchromosome_read_pairs);

  std::ofstream myfile;
  myfile.open("matrix.csv");

  std::string output = ",";
  for (int i = 0; i < static_cast<int>(matrix.size()); i++) {
    output += std::to_string(i) + ",";
  }
  output += "Total\n";
  std::cerr << output << std::endl;
  myfile << output;

  for (int i = 0; i < static_cast<int>(matrix.size()); i++) {
    output ="";
    output += std::to_string(i) + ",";
    for (int r = 0; r < static_cast<int>(matrix[i].size()); r++) {
      output += std::to_string(matrix[i][r]) + ",";
    }
    it = contigs.find(i);
    if (it != contigs.end()) {
      output += std::to_string(it->second.interchromosome_links)+"\n";
    }
    myfile << output;
    std::cerr << output << std::endl;
  }
  myfile.close();

  for (const auto& contig : contigs) {
    sum_interchromosome_links += contig.second.interchromosome_links;
    sum_intrachromosome_links += contig.second.intrachromosome_links;
    std::cerr << "Contig = " << contig.first
              << " Intrachromosome links = " << contig.second.intrachromosome_links
              << " Interchromosome links = " << contig.second.interchromosome_links
              << std::endl;
  }
  std::cerr << "Total Intrachromosome links = " << sum_intrachromosome_links
            << " Total Interchromosome links = " << sum_interchromosome_links/2
            << std::endl;
  
  /*
  std::ofstream myfile;
  myfile.open ("log.txt");
  for (const auto& rp : multiple_overlap_read_pairs){
    myfile << "read pair: " << rp.first
           << " id: " << rp.second[0][0].lhs_id << "\n";
    
    for (const auto& overlaps : rp.second){
      for (const auto& overlap : overlaps){
        myfile << " rhs id: " << overlap.rhs_id
               << " begin: " << overlap.rhs_begin
               << " end: " << overlap.rhs_end << "\n"; 
      }
    }
  
  }
  myfile.close();*/
  return;
}

// technically this also can be multi-thread - just remove the for loop
void Graph::CreateGraph(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>& targets) {
  for (auto const& target : targets) {
    contigs.emplace(target->id, Node(target->id, target->inflated_len));
  }
}

void Graph::GenerateMatrix(
  std::vector<std::vector<std::uint32_t>> &matrix, 
  std::unordered_map<std::string, std::vector<std::vector<biosoup::Overlap>>>& interchromsome_read_pairs) {
  for (const auto& rp : interchromsome_read_pairs) {
    auto id_1 = rp.second[0][0].rhs_id;
    auto id_2 = rp.second[1][0].rhs_id;
    matrix[id_1][id_2] += 1;
    matrix[id_2][id_1] += 1;
  }

  for (const auto& contig : contigs) {
    auto id = contig.first;
    matrix[id][id] = contig.second.intrachromosome_links;
  }
}

void Graph::CalcualteInterChromosomeLinks(
  std::unordered_map<std::string, std::vector<std::vector<biosoup::Overlap>>>& interchromsome_read_pairs) {
  std::unordered_map<std::uint32_t, Node>::iterator found;
  for (const auto& rp : interchromsome_read_pairs) {
    found = contigs.find(rp.second[0][0].rhs_id);
    if (found == contigs.end()) {
      std::cerr << "ERROR contig not found" << std::endl;
    } else {
      found->second.interchromosome_links++;
    }
    found = contigs.find(rp.second[1][0].rhs_id);
    if (found == contigs.end()) {
      std::cerr << "ERROR contig not found" << std::endl;
    } else {
      found->second.interchromosome_links++;
    }
  }
}

void Graph::FillPileogram(std::unordered_map<std::string, std::vector<std::vector<biosoup::Overlap>>>& read_pairs) {  // NOLINT
  std::unordered_map<std::uint32_t, Node>::iterator found;
  std::pair<std::uint32_t, std::uint32_t> overlap;
  std::uint32_t min_overlap = UINT32_MAX, max_overlap = 0, average_overlap = 0;
  std::unordered_map<std::uint32_t, std::vector<std::pair<std::uint32_t, std::uint32_t>>> overlap_map;
  std::unordered_map<std::uint32_t, std::vector<std::pair<std::uint32_t, std::uint32_t>>>::iterator overlap_map_iter;
  // find the contig, then add layer to pileogram
  for (const auto& rp : read_pairs) {
    found = contigs.find(rp.second[0][0].rhs_id);
    if (found == contigs.end()) {
      // not found
      std::cerr << "ERROR contig not found" << std::endl;
    } else {
      found->second.intrachromosome_links++;
      overlap = GetOverlap(rp.second[0][0], rp.second[1][0]);
      auto length = std::get<1>(overlap) - std::get<0>(overlap);
      min_overlap = (length < min_overlap) ? length:min_overlap;
      max_overlap = (length > max_overlap) ? length:max_overlap;
      average_overlap += length;

      // throw into overlap map
      overlap_map_iter = overlap_map.find(rp.second[0][0].rhs_id);
      if (overlap_map_iter == overlap_map.end()) {
        // not found, create map
        std::vector<std::pair<std::uint32_t, std::uint32_t>> temp;
        temp.emplace_back(overlap);
        overlap_map.insert(std::make_pair(rp.second[0][0].rhs_id, temp));
      } else {
        overlap_map_iter->second.emplace_back(overlap);
      }
    }
  }
  for (overlap_map_iter = overlap_map.begin(); overlap_map_iter != overlap_map.end(); overlap_map_iter++) {
    found = contigs.find(overlap_map_iter->first);
    found->second.pileogram.AddLayer(overlap_map_iter->second);
  }
  average_overlap /= read_pairs.size();
  std::cerr << "[tarantula::Construct] Stats: " << std::endl;
  std::cerr << "[tarantula::Construct] min overlap: " << min_overlap << std::endl;
  std::cerr << "[tarantula::Construct] max overlap: " << max_overlap << std::endl;
  std::cerr << "[tarantula::Construct] average overlap: " << average_overlap << std::endl;
}

std::pair<std::uint32_t, std::uint32_t> Graph::GetOverlap(
    biosoup::Overlap ol1,
    biosoup::Overlap ol2) {
  return std::make_pair(
      std::min(ol1.rhs_begin, ol2.rhs_begin),
      std::max(ol1.rhs_end,   ol2.rhs_end));
}

void Graph::Process(
  std::vector<std::future<std::pair<std::string,
  std::vector<std::vector<biosoup::Overlap>>>>>& futures,
  ram::MinimizerEngine& minimizer_engine,
  std::unique_ptr<biosoup::NucleicAcid>& sequence1,
  std::unique_ptr<biosoup::NucleicAcid>& sequence2) {
  futures.emplace_back(thread_pool_->Submit([&] (
    ram::MinimizerEngine& minimizer_engine,
    const std::unique_ptr<biosoup::NucleicAcid>& sequence1,
    const std::unique_ptr<biosoup::NucleicAcid>& sequence2)
    -> std::pair<std::string, std::vector<std::vector<biosoup::Overlap>>> {
      std::vector<std::vector<biosoup::Overlap>> minimizer_result;
      minimizer_result.emplace_back(minimizer_engine.Map(sequence1, false, false));
      minimizer_result.emplace_back(minimizer_engine.Map(sequence2, false, false));
      auto long_read_len = sequence1->inflated_len*0.8;

      if (minimizer_result[0].size() < 1 || minimizer_result[1].size() < 1) {
        // no overlap, discard
        minimizer_result.clear();
        return {"empty", minimizer_result};
      }

      for (auto& rp : minimizer_result) {
        if (rp.size() == 1) {
          continue;
        }
        std::vector<int> discard;
        biosoup::Overlap temp;
        auto numLR = 0;
        for (auto const& ol : rp) {
          if (ol.rhs_end-ol.rhs_begin > long_read_len) {
            numLR++;
            temp = ol;
          }
          if (numLR > 1) {
            // if there are more than 1 long read discard
            return {sequence1->name + "_multiple", minimizer_result};
          }
        }
        if (numLR == 0) {
          // dicard too because all are short reads
          return {sequence1->name + "_multiple", minimizer_result};
        }

        rp.clear();
        rp.emplace_back(temp);
      }

      if (minimizer_result[0][0].rhs_id != minimizer_result[1][0].rhs_id) {
        // interchromosome
        return {sequence1->name + "_interchromosome", minimizer_result};
      }

      // return pair
      return {sequence1->name, minimizer_result};
    },
    std::ref(minimizer_engine),
    std::cref(sequence1),
    std::cref(sequence2)));
}

void Graph::PrintJson(const std::string& path) const {
  if (path.empty()) {
    return;
  }

  std::ofstream os(path);
  cereal::JSONOutputArchive archive(os);
  for (const auto& it : contigs) {
    archive(cereal::make_nvp(std::to_string(it.second.id), it.second.pileogram));  // NOLINT
  }
}

}  // namespace tarantula
