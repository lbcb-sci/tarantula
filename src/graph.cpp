// Copyright (c) 2021 Cecilia Lee, Robert Vaser

#include "graph.hpp"

#include <fstream>
#include <iostream>
#include <queue>

#include "biosoup/timer.hpp"
#include "cereal/archives/json.hpp"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace tarantula {

Graph::Graph(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool, 
    int window_size)
    : thread_pool_(thread_pool ?
        thread_pool :
        std::make_shared<thread_pool::ThreadPool>(1)), 
      window_size(window_size) {}

void Graph::Construct(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>& targets,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences) {
  // paramters for RAM
  uint32_t k = 21, w = 11, bandwidth = 100, chain = 2, matches = 25, gap = 100;
  double frequency = 0.00001;

  // statistics
  int num_pair = 0;
  int num_filtered_pair = 0;
  int discard_due_to_multiple_sr = 0;
  int discard_due_to_1_multiple_lr = 0;
  int discard_due_to_2_multiple_lr = 0;
  int discard_due_to_no_overlap = 0;
  int only_1_pair_overlap = 0;

  int multiple_sr_mt_4 = 0;
  int multiple_sr_mt_0 = 0;
  int no_ol_1_mt_4 = 0;
  int no_ol_1_mt_0 = 0;
  int no_ol_2_mt_4 = 0;
  int no_ol_2_mt_0 = 0;

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

  std::vector<std::future<std::vector<std::pair<std::string, std::vector<biosoup::Overlap>>>>> futures;
  std::unordered_map<std::string, std::vector<biosoup::Overlap>> multiple_overlap_read_pairs;
  //std::unordered_map<std::string, std::vector<biosoup::Overlap>> read_pairs;
  //std::unordered_map<std::string, std::vector<biosoup::Overlap>> interchromosome_read_pairs;
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
  std::ofstream myfile, myfile2, myfile3;
  std::ofstream inter, intra;
  // inter.open("interchromosome_strand.txt");
  // intra.open("intrachromosome_strand.txt");
  std::vector<int> window_id_map = GenerateMapWindowID();
  int total_windows = GetNumWindows();
  std::cerr << "[tarantula::Construct] Number of windows = " << total_windows << std::endl;
  std::vector<std::vector<std::uint32_t>> window_matrix(total_windows, std::vector<std::uint32_t>(total_windows, 0));
  std::cerr << "[tarantula::Construct] Total number of sequence: "
            << sequences.size() << std::endl;
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
        auto result_vector = it.get();
        for (auto result : result_vector ) {
          if (result.first.find("_interchromosome") != std::string::npos) {
            interchromosome_read_pairs.insert(result);
          } else if (result.first.find("_multiple_short_reads") != std::string::npos) {
            discard_due_to_multiple_sr++;
            multiple_overlap_read_pairs.insert(result);

            if (result.first.find("_multiple_short_reads_mt_4") != std::string::npos) {
              multiple_sr_mt_4++;
              multiple_sr_mt_0++;
            } else if (result.first.find("_multiple_short_reads_mt_0") != std::string::npos)
              multiple_sr_mt_0++;

          } else if (result.first.find("_multiple_long_reads_1_pair") != std::string::npos) {
            discard_due_to_1_multiple_lr++;
            multiple_overlap_read_pairs.insert(result);
          } else if (result.first.find("_multiple_long_reads_2_pair") != std::string::npos) {
            discard_due_to_2_multiple_lr++;
            multiple_overlap_read_pairs.insert(result);
          } else if (result.first.find("1_overlap") != std::string::npos) {
            only_1_pair_overlap++;
            discard_due_to_no_overlap++;

            if (result.first.find("1_overlap_mt_4") != std::string::npos) {
              no_ol_1_mt_4++;
              no_ol_1_mt_0++;
            } else if (result.first.find("1_overlap_mt_0") != std::string::npos) {
              no_ol_1_mt_0++; 
            }
              

          } else if (result.first.find("empty") != std::string::npos) {
            discard_due_to_no_overlap++;
            if (result.first.find("empty_mt_4") != std::string::npos) {
              no_ol_2_mt_4++;
              no_ol_2_mt_0++;
            } else if (result.first.find("empty_mt_0") != std::string::npos) {
              no_ol_2_mt_0++; 
            }
          } else if (result.first.compare("discard") != 0) {
            read_pairs.insert(result);
          }
        }
        
      }
      num_filtered_pair += read_pairs.size();
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

      // interwindow links
      GenerateMatrixWindowIntraLinks(window_id_map, window_matrix, read_pairs);

      // discard read pair
      read_pairs.clear();
      futures.clear();
      bytes = 0;
    }
  }
  // inter.close();
  // intra.close();
  std::cerr << "[tarantula::Construct]-----------FILTERING STATISTICS-----------------" << std::endl;
  std::cerr << "[tarantula::Construct] Total number of reads: "
            << sequences.size() << std::endl;
  std::cerr << "[tarantula::Construct] Number of reads that are not a pair: "
            << sequences.size() - num_pair*2 << std::endl;
  std::cerr << "[tarantula::Construct] Number of pairs that are sucessfully filtered: "
            << num_filtered_pair << std::endl;
  std::cerr << "[tarantula::Construct] Number of pairs that are interchomosome: "
            << interchromosome_read_pairs.size() << std::endl;
  std::cerr << "[tarantula::Construct] Number of pairs that discard as there is multiple short reads: "
            << discard_due_to_multiple_sr << std::endl;
  std::cerr << "[tarantula::Construct] Number of pairs that discard as there is multiple short reads, filtered size > 0: "
            << multiple_sr_mt_0 << std::endl;
  std::cerr << "[tarantula::Construct] Number of pairs that discard as there is multiple short reads, filtered size > 4: "
            << multiple_sr_mt_4 << std::endl << std::endl;
  std::cerr << "[tarantula::Construct] Number of pairs that discard as there is 1 pair with multiple long reads: "
            << discard_due_to_1_multiple_lr << std::endl;
  std::cerr << "[tarantula::Construct] Number of pairs that discard as there is 2 pair with multiple long reads: "
            << discard_due_to_2_multiple_lr << std::endl << std::endl;
  std::cerr << "[tarantula::Construct] Number of pairs that discard as there is only 1 pair with overlaps: "
            << only_1_pair_overlap << std::endl;
  std::cerr << "[tarantula::Construct] Number of pairs that discard as there is only 1 pair with overlaps, filtered size > 0: "
            << no_ol_1_mt_0 << std::endl;
  std::cerr << "[tarantula::Construct] Number of pairs that discard as there is only 1 pair with overlaps, filtered size > 4: "
            << no_ol_1_mt_4 << std::endl << std::endl;
  std::cerr << "[tarantula::Construct] Number of pairs that discard as both pair hav no overlaps: "
            << discard_due_to_no_overlap - only_1_pair_overlap << std::endl;
  std::cerr << "[tarantula::Construct] Number of pairs that discard as both pair hav no overlaps, filtered size > 0: "
            << no_ol_2_mt_0 << std::endl;
  std::cerr << "[tarantula::Construct] Number of pairs that discard as both pair hav no overlaps, filtered size > 4: "
            << no_ol_2_mt_4 << std::endl;
  std::cerr << "[tarantula::Construct]-------------------------------------------------" << std::endl;

  // save overlaps
  Store();

  CalcualteInterChromosomeLinks(interchromosome_read_pairs);
  uint32_t sum_interchromosome_links = 0, sum_intrachromosome_links = 0;
  int num_chromosomes = targets.size();
  std::vector<std::vector<std::uint32_t>> matrix(num_chromosomes, std::vector<std::uint32_t>(num_chromosomes, 0));
  std::unordered_map<std::uint32_t, Node>::iterator it;
  GenerateMatrix(matrix, interchromosome_read_pairs);

  // write matrix of interchromosome links to csv
  /*
  myfile.open("matrix.csv");

  std::string output = ",";
  for (int i = 0; i < static_cast<int>(matrix.size()); i++) {
    output += std::to_string(i) + ",";
  }
  output += "Total, 0011 number, 0110 number\n";
  myfile << output;

  for (int i = 0; i < static_cast<int>(matrix.size()); i++) {
    output ="";
    output += std::to_string(i) + ",";
    for (int r = 0; r < static_cast<int>(matrix[i].size()); r++) {
      output += std::to_string(matrix[i][r]) + ",";
    }
    it = contigs.find(i);
    if (it != contigs.end()) {
      output += std::to_string(it->second.interchromosome_links)+"," +
                std::to_string(it->second.link_0011)+","+
                std::to_string(it->second.link_0110)+"\n";
    }
    myfile << output;
  }
  myfile.close();*/

  
  // contig window links text file
  /*
  int sum = 0;
  myfile.open("contig_window_links.txt");
  for (auto contig : contigs) {
    for (auto window : contig.second.windows) {
      myfile <<  "contig = " << contig.first << " , window = " << window.id << ", links = " << window.interchromosome_links << "\n";
      sum+=window.interchromosome_links;
    }
  }
  myfile.close();*/

  // window matrix csv
  /* 
  GenerateMatrixWindow(window_id_map, window_matrix, interchromosome_read_pairs);
  std::cerr << "[tarantula::Construct] Generation of window matrix done" << std::endl;
  myfile.open("window_matrix.csv");
  
  output = ",";
  for (int i = 0; i < static_cast<int>(window_matrix.size()); i++) {
    output += std::to_string(i) + ",";
  }
  output += "Total\n";
  myfile << output;

  for (int i = 0; i < static_cast<int>(window_matrix.size()); i++) {
    auto sum = 0;
    output ="";
    output += std::to_string(i);
    for (int r = 0; r < static_cast<int>(window_matrix[i].size()); r++) {
      output += "," + std::to_string(window_matrix[i][r]);
      sum += window_matrix[i][r];
    }
    output += "," + std::to_string(sum)+"\n";
    myfile << output;
  }
  myfile.close(); */

  for (const auto& contig : contigs) {
    sum_interchromosome_links += contig.second.interchromosome_links;
    sum_intrachromosome_links += contig.second.intrachromosome_links;
    std::cerr << "[tarantula::Construct] Contig = " << contig.first
              << " Intrachromosome links = " << contig.second.intrachromosome_links
              << " Interchromosome links = " << contig.second.interchromosome_links
              << std::endl;
  }
  std::cerr << "[tarantula::Construct] Total Intrachromosome links = "
            << sum_intrachromosome_links
            << " Total Interchromosome links = " << sum_interchromosome_links/2
            << std::endl;

  int start, end;

  for (std::uint32_t i = 0; i < window_id_map.size(); i++) {
    myfile.open("contig_" + std::to_string(i) + ".txt");
    start = window_id_map[i];
    if (i+1 == window_id_map.size()) {
      end = window_matrix.size();
    } else {
      end = window_id_map[i+1];
    }
    for (int r = start; r < end; r++) {
      for (int x = start; x < r; x++) {
        if (window_matrix[r][x] != 0)
          myfile << r << "--" << x << "," << window_matrix[r][x] << "\n";
      }
    }
    myfile.close();
  }

  // directed force
  /*
  for (std::uint32_t i = 0; i < window_id_map.size(); i++) {
    std::string input = "contig_" + std::to_string(i) + ".txt";
    std::string output = "contig_" + std::to_string(i) + "_output.txt";
    directedforce::GenerateGraphFromDirectedForceAlgorithm(input, output);
  }*/

  


  // graph -- 1 contig == 1 node
  /*
  // find component & only draw components that have >= 3 nodes
  std::vector<std::vector<std::uint32_t>> components = GetComponents(matrix);
  for (int it = 0; it < static_cast<int>(components.size()); it++) {
    if (components[it].size() <= 2)
      continue;
    myfile.open("component_" + std::to_string(it) + ".txt");
    for (int i = 0; i < static_cast<int>(components[it].size()); i++) {
      for (int r = 0; r < i; r++) {
        if (matrix[components[it][i]][components[it][r]] != 0)
          myfile << i << "--" << r << "," << matrix[components[it][i]][components[it][r]] << "\n";
      }
    }
    myfile.close();
  }

  // create text file to visualise graph
  myfile.open("graph_max.txt");
  myfile2.open("graph_inter_max.txt");
  myfile3.open("graph.txt");
  std::uint32_t max_contig = MaxInMatrix(matrix);
  std::uint32_t max_contig_inter = MaxInter(matrix);
  for (int i = 0; i < static_cast<int>(matrix.size()); i++) {
    for (int r = 0; r < i; r++) {
      if (matrix[i][r] == 0)
        continue;
      double weight = (double) matrix[i][r]/ (double) max_contig;
      if (weight > 0.0001)
        myfile << i << "--" << r << "," << weight << "\n";
      double weight_inter = (double) matrix[i][r]/ (double) max_contig_inter;
      if (weight_inter > 0.0001)
        myfile2 << i << "--" << r << "," << weight_inter << "\n";
      myfile3 << i << "--" << r << "," << matrix[i][r] << "\n";
    }
  }
  myfile.close();
  myfile2.close();
  myfile3.close();*/

  // find component & only draw components that have >= 3 nodes (BY WINDOW)  -- error in get components, create new method for it
  /*
  for (int it = 0; it < static_cast<int>(components.size()); it++) {
    if (components[it].size() <= 2)
      continue;
    myfile.open("window_component_" + std::to_string(it) + ".txt");

    // assume contig id is from 0 - contigs.size()-1

    int window_in_contig = 0;
    int counter = 1;
    for (int i = 0; i < window_matrix.size(); i++) {
      if (i >= window_id_map[counter]) {
        window_in_contig++;
        counter++;
      }
      myfile << i << " " << window_in_contig << "\n";
    }
    myfile << "%\n";
    for (int i = 0; i < static_cast<int>(components[it].size()); i++) {
      for (int r = 0; r < i; r++) {
        int contig_id_1 = components[it][i];
        int contig_id_2 = components[it][r];
        for (int y = window_id_map[contig_id_1]; y < window_id_map[contig_id_1+1]; y++) {
          for (int x = window_id_map[contig_id_2]; x < window_id_map[contig_id_2+1]; x++) {
            if (window_matrix[y][x]!=0)
              myfile << i << "--" << r << "," << window_matrix[y][x] << "\n";
          }
        }
      }
    }
    // deal with all the links between windows in each contig
    for (int i = 0; i < static_cast<int>(components[it].size()); i++) {
      int contig_id = components[it][i];
      for (int r = window_id_map[contig_id]; r < window_id_map[contig_id+1]-1; r++){
        if (window_matrix[r][r+1]!=0)
          myfile << i << "--" << r << "," << window_matrix[r][r+1] << "\n";
      }
    }
    myfile.close();
  }*/

// window graph -- inter and intra links
/*
  int window_in_contig = 0;
  int counter = 1;
  myfile.open("window_graph_max.txt");
  myfile2.open("window_graph_max_inter.txt");
  myfile3.open("window_graph.txt");
  for (int i = 0; i < window_matrix.size(); i++) {
    if (i >= window_id_map[counter]) {
      window_in_contig++;
      counter++;
    }
    myfile << i << " " << window_in_contig << "\n";
    myfile2 << i << " " << window_in_contig << "\n";
    myfile3 << i << " " << window_in_contig << "\n";
  }
  myfile << "%\n";
  myfile2 << "%\n";
  myfile3 << "%\n";

  // scale inter and intra the same way
  
  std::uint32_t max = MaxInMatrix(window_matrix);
  std::uint32_t max_inter = MaxInter(window_matrix);
  std::cerr << "max: " << max << std::endl;
  for (int i = 0; i < static_cast<int>(window_matrix.size()); i++) {
    for (int r = 0; r < i; r++) {
      if (window_matrix[i][r] == 0)
        continue;
      double weight = (double) window_matrix[i][r]/ (double) max;
      if (weight > 0.0001 ) {
        myfile << i << "--" << r << "," << weight << "\n";
      }
      double weight_inter = (double) window_matrix[i][r]/ (double) max_inter;
      if (weight > 0.0001 ) {
        myfile2 << i << "--" << r << "," << weight << "\n";
      }
      myfile3 << i << "--" << r << "," << window_matrix[i][r] << "\n";
    }
  }

  myfile.close();
  myfile2.close();
  myfile3.close();*/

  // get pairs are discarded because there is multiple overlap
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
    contigs.emplace(target->id, Node(target->id, target->inflated_len, window_size));
  }
}

void Graph::GenerateMatrix(
  std::vector<std::vector<std::uint32_t>> &matrix,
  std::unordered_map<std::string, std::vector<biosoup::Overlap>>& interchromsome_read_pairs) {
  for (const auto& rp : interchromsome_read_pairs) {
    auto id_1 = rp.second[0].rhs_id;
    auto id_2 = rp.second[1].rhs_id;
    matrix[id_1][id_2] += 1;
    matrix[id_2][id_1] += 1;
  }

  for (const auto& contig : contigs) {
    auto id = contig.first;
    matrix[id][id] = contig.second.intrachromosome_links;
  }
}

std::vector<int> Graph::GenerateMapWindowID() {
  std::vector<int> window_id_map(contigs.size(), 0);
  std::unordered_map<std::uint32_t, Node>::iterator found;
  int sum = 0;

  std::cerr << "[tarantula::Construct] Number of contigs: " << contigs.size() << std::endl;
  for (std::uint32_t i = 0; i < contigs.size(); i++) {
    found = contigs.find(i);
    window_id_map[found->first] = sum;
    std::cerr << "[tarantula::Construct] contig "<< i << " : window index start from = " << sum << std::endl;
    sum += found->second.windows.size();
  }
  return window_id_map;
}

void Graph::GenerateMatrixWindowIntraLinks(
  std::vector<int>& window_id_map,
  std::vector<std::vector<std::uint32_t>> &matrix,
  std::unordered_map<std::string, std::vector<biosoup::Overlap>>& read_pairs) {
  for (const auto& rp : read_pairs) {
    int window_index_begin_1 = rp.second[0].rhs_begin/window_size;
    int window_index_begin_2 = rp.second[1].rhs_begin/window_size;
    int window = window_id_map[rp.second[1].rhs_id];
    int window_id_1 = window + window_index_begin_1;
    int window_id_2 = window + window_index_begin_2;
    matrix[window_id_1][window_id_2]+=1;
    matrix[window_id_2][window_id_1]+=1;
  }
}

std::uint32_t Graph::MaxInMatrix(std::vector<std::vector<std::uint32_t>> &matrix) {
  std::uint32_t max = 0;
  for (std::uint32_t i = 0; i < matrix.size(); i++) {
    for (std::uint32_t r = 0; r <= i; r++) {
      if (matrix[i][r] > max) {
        max = matrix[i][r];
      }
    }
  }
  return max;
}

std::uint32_t Graph::MaxInter(std::vector<std::vector<std::uint32_t>> &matrix) {
  std::uint32_t max = 0;
  for (std::uint32_t i = 0; i < matrix.size(); i++) {
    for (std::uint32_t r = 0; r < i; r++) {
      if (i == r)
        continue;
      if (matrix[i][r] > max) {
        max = matrix[i][r];
      }
    }
  }
  return max;
}

std::uint32_t MaxIntra(std::vector<std::vector<std::uint32_t>> &matrix) {
  std::uint32_t max = 0;
  for (std::uint32_t i = 0; i < matrix.size(); i++) {
    if (matrix[i][i] > max) {
      max = matrix[i][i];
    }
  }
  return max;
}

void Graph::GenerateMatrixWindow(
  std::vector<int>& window_id_map,
  std::vector<std::vector<std::uint32_t>> &matrix,
  std::unordered_map<std::string, std::vector<biosoup::Overlap>>& interchromsome_read_pairs) {
  std::unordered_map<std::uint32_t, Node>::iterator found;
  int extra = 0;

  for (const auto& rp : interchromsome_read_pairs) {
    int window_index_begin_1 = rp.second[0].rhs_begin/window_size;
    int window_index_end_1 = rp.second[0].rhs_end/window_size;
    int window_index_begin_2 = rp.second[1].rhs_begin/window_size;
    int window_index_end_2 = rp.second[1].rhs_end/window_size;
    int id_1 = rp.second[0].rhs_id;
    int id_2 = rp.second[1].rhs_id;

    int window_id_1_begin = window_index_begin_1+window_id_map[id_1];
    int window_id_2_begin = window_index_begin_2+window_id_map[id_2];

    if (window_index_begin_1 != window_index_end_1) {
      int window_id_1_end = window_index_end_1 + window_id_map[id_1];
      if (window_index_begin_2 != window_index_end_2) {
        int window_id_2_end = window_index_end_2 + window_id_map[id_2];
        matrix[window_id_1_end][window_id_2_end] += 1;
        matrix[window_id_2_end][window_id_1_end] += 1;
      } else {
        matrix[window_id_1_end][window_id_2_begin] += 1;
        matrix[window_id_2_begin][window_id_1_end] += 1;
      }
      extra++;
    } else if (window_index_begin_2 != window_index_end_2) {
      int window_id_2_end = window_index_end_2 + window_id_map[id_2];
      matrix[window_id_1_begin][window_id_2_end] += 1;
      matrix[window_id_2_end][window_id_1_begin] += 1;
      extra++;
    }
    matrix[window_id_1_begin][window_id_2_begin] += 1;
    matrix[window_id_2_begin][window_id_1_begin] += 1;
  }
  
  // link between windows in the same contig -- take the max interchromosome link
  /*
  for (int i = 0; i < window_id_map.size(); i++) {
    int id = window_id_map[i];
    int max = *max_element(std::begin(matrix[id]), std::end(matrix[id]));
    int next = i+1;
    int end;
    if (next != window_id_map.size()) {
      end = window_id_map[next]-1;
    } else {
      end = matrix.size()-1;
    }

    for (int r = id; r < end; r++) {
      matrix[r][r+1] = max;
      matrix[r+1][r] = max;
    }
  }*/

  // intrachromosome links
  /*
  int counter = 0;
  for (int i = 0; i < contigs.size(); i++) {
    found = contigs.find(i);
    int window_id = window_id_map[i];
    for (int r = 0; r < found->second.windows.size(); r++) {
      matrix[window_id+r][window_id+r] = found->second.windows[r].intrachromosome_links;
    }
  }*/
}

void Graph::CalcualteInterChromosomeLinks(
  std::unordered_map<std::string, std::vector<biosoup::Overlap>>& interchromsome_read_pairs) {
  std::unordered_map<std::uint32_t, Node>::iterator found;
  int window_index_begin, window_index_end;
  int strand_1, strand_2;
  int extra = 0;
  for (const auto& rp : interchromsome_read_pairs) {
    strand_1 = rp.second[0].strand;
    strand_2 = rp.second[1].strand;
    window_index_begin = rp.second[0].rhs_begin/window_size;
    window_index_end = rp.second[0].rhs_end/window_size;
  
    found = contigs.find(rp.second[0].rhs_id);
    if (found == contigs.end()) {
      std::cerr << "ERROR contig not found" << std::endl;
    } else {
      if (window_index_begin != window_index_end) {
        found->second.windows[window_index_end].interchromosome_links++;
        extra++;
      }
      found->second.windows[window_index_begin].interchromosome_links++;
      found->second.interchromosome_links++;

      if (strand_1 == strand_2) {
        found->second.link_0011++;
      } else {
        found->second.link_0110++;
      }
    }
    if (found->second.windows.size() <= window_index_end) {
      std::cerr << "window size:" << found->second.windows.size() <<"window index: " << window_index_end << std::endl;
    }

    window_index_begin = rp.second[1].rhs_begin/window_size;
    window_index_end = rp.second[1].rhs_end/window_size;
    found = contigs.find(rp.second[1].rhs_id);
    if (found == contigs.end()) {
      std::cerr << "ERROR contig not found" << std::endl;
    } else {
      if (window_index_begin != window_index_end) {
        found->second.windows[window_index_end].interchromosome_links++;
        extra++;
      }
      found->second.windows[window_index_begin].interchromosome_links++;
      found->second.interchromosome_links++;

      if (strand_1 == strand_2) {
        found->second.link_0011++;
      } else {
        found->second.link_0110++;
      }
    }
    if (found->second.windows.size() <= window_index_end) {
      std::cerr << "window size:" << found->second.windows.size() <<"window index: " << window_index_end << std::endl;
    }
  }
}

void Graph::FillPileogram(std::unordered_map<std::string, std::vector<biosoup::Overlap>>& read_pairs) {  // NOLINT
  std::unordered_map<std::uint32_t, Node>::iterator found;
  std::pair<std::uint32_t, std::uint32_t> overlap;
  std::uint32_t min_overlap = UINT32_MAX, max_overlap = 0, average_overlap = 0;
  std::unordered_map<std::uint32_t, std::vector<std::pair<std::uint32_t, std::uint32_t>>> overlap_map;
  std::unordered_map<std::uint32_t, std::vector<std::pair<std::uint32_t, std::uint32_t>>>::iterator overlap_map_iter;
  // find the contig, then add layer to pileogram
  int extra = 0;
  for (const auto& rp : read_pairs) {
    found = contigs.find(rp.second[0].rhs_id);
    if (found == contigs.end()) {
      // not found
      std::cerr << "ERROR contig not found" << std::endl;
    } else {
      // update intrachromosome links
      found->second.intrachromosome_links++;
      
      // overlap
      overlap = GetOverlap(rp.second[0], rp.second[1]);
      auto length = std::get<1>(overlap) - std::get<0>(overlap);
      min_overlap = (length < min_overlap) ? length:min_overlap;
      max_overlap = (length > max_overlap) ? length:max_overlap;
      average_overlap += length;

      // throw into overlap map
      overlap_map_iter = overlap_map.find(rp.second[0].rhs_id);
      if (overlap_map_iter == overlap_map.end()) {
        // not found, create map
        std::vector<std::pair<std::uint32_t, std::uint32_t>> temp;
        temp.emplace_back(overlap);
        overlap_map.insert(std::make_pair(rp.second[0].rhs_id, temp));
      } else {
        overlap_map_iter->second.emplace_back(overlap);
      }
    }
  }
  std::cerr << "[tarantula::Construct] Extra intrachromosome links between windows: " << extra << std::endl;
  
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


// map hi-c reads to contig & filter
void Graph::Process(
  std::vector<std::future<std::vector<std::pair<std::string,std::vector<biosoup::Overlap>>>>>& futures,
  ram::MinimizerEngine& minimizer_engine,
  std::unique_ptr<biosoup::NucleicAcid>& sequence1,
  std::unique_ptr<biosoup::NucleicAcid>& sequence2) {

  futures.emplace_back(thread_pool_->Submit([&] (
    ram::MinimizerEngine& minimizer_engine,
    const std::unique_ptr<biosoup::NucleicAcid>& sequence1,
    const std::unique_ptr<biosoup::NucleicAcid>& sequence2)
    -> std::vector<std::pair<std::string, std::vector<biosoup::Overlap>>>{
      std::vector<std::uint32_t> filtered1, filtered2;
      std::vector<std::vector<biosoup::Overlap>> minimizer_result, minimizer_result_long_read;
      minimizer_result.emplace_back(minimizer_engine.Map(sequence1, false, false, false, &filtered1));
      minimizer_result.emplace_back(minimizer_engine.Map(sequence2, false, false, false, &filtered2));
      auto long_read_len = sequence1->inflated_len*0.75;

      // no overlap for 1 / both
      if (minimizer_result[0].size() < 1 || minimizer_result[1].size() < 1) {
        // no overlap, discard
        if (minimizer_result[0].size() > 0 || minimizer_result[1].size() > 0) {
          // minimizer_result.clear();
          
          if (filtered1.size() > 4 || filtered2.size() > 4)
            return {{"1_overlap_mt_4", {}}};
          if (filtered1.size() > 0 || filtered2.size() > 0)
            return {{"1_overlap_mt_0", {}}};

          //minimizer_result.clear();
          return {{"1_overlap", {}}};
        } else {
          // minimizer_result.clear();
          
          if (filtered1.size() > 4 || filtered2.size() > 4)
            return {{"empty_mt_4", {}}};
          if (filtered1.size() > 0 || filtered2.size() > 0)
            return {{"empty_mt_0", {}}};
          //minimizer_result.clear();
          return {{"empty", {}}};
        }
      }
      
      // filter out those len > 0.75
      for (auto& minimizer : minimizer_result) {
        std::vector<biosoup::Overlap> filtered_overlaps;
        for (auto& overlap: minimizer) {
          if (overlap.rhs_end - overlap.rhs_begin > long_read_len) {
            filtered_overlaps.push_back(overlap);
          }
        }
        minimizer_result_long_read.push_back(filtered_overlaps);
      }

      // no long reads for 1 / both
      if (minimizer_result_long_read[0].size() == 0 || minimizer_result_long_read[1].size() == 0) {
        if (filtered1.size() > 4 || filtered2.size() > 4)
          return {{sequence1->name + "__multiple_short_reads_mt_4", {}}};
        if (filtered1.size() > 0 || filtered2.size() > 0)
          return {{sequence1->name + "_multiple_short_reads_mt_0", {}}};
        return {{sequence1->name + "_multiple_short_reads", {}}};
      }
        
      // straight forward interchromosome / intrachromosome
      if (minimizer_result_long_read[0].size() == 1 && minimizer_result_long_read[1].size() == 1) {
        std::vector<biosoup::Overlap> result;
        result.emplace_back(minimizer_result_long_read[0][0]);
        result.emplace_back(minimizer_result_long_read[1][0]);
        if (minimizer_result_long_read[0][0].rhs_id == minimizer_result_long_read[1][0].rhs_id) {
          return {{sequence1->name, result}};
        } else {
          return {{sequence1->name + "_interchromosome", result}};
        }
      }


      // multiple long reads
      std::unordered_map<std::uint32_t, std::tuple<std::vector<biosoup::Overlap>, std::vector<biosoup::Overlap>>> intra_pairs;
      std::unordered_map<std::uint32_t, std::tuple<std::vector<biosoup::Overlap>, std::vector<biosoup::Overlap>>>::iterator iter;
      int pair_count = 0;
      for (auto& minimizer_result_each_pair : minimizer_result_long_read) {
        for (auto& overlap : minimizer_result_each_pair) {
          iter = intra_pairs.find(overlap.rhs_id);
          if (iter == intra_pairs.end()) {
            if (pair_count == 1) {
              continue;
            }
            std::vector<biosoup::Overlap> vec1 = {overlap};
            std::vector<biosoup::Overlap> vec2 = {};
            std::tuple<std::vector<biosoup::Overlap>, std::vector<biosoup::Overlap>> temp_tuple = std::make_tuple(vec1, vec2);
            intra_pairs.insert({overlap.rhs_id, temp_tuple});
          } else {
            if (pair_count == 0)
              std::get<0>(iter->second).push_back(overlap);
            else
              std::get<1>(iter->second).push_back(overlap);
          }
        }
        pair_count++;
      }

      
      std::vector<std::pair<std::string, std::vector<biosoup::Overlap>>> result;
      bool there_is_result = false;
      for (auto const& pair : intra_pairs) {
        if (std::get<0>(pair.second).size() != 0 && std::get<1>(pair.second).size() != 0) {
          there_is_result = true;
          for (auto const& ol1 : std::get<0>(pair.second)) {
            for (auto const& ol2 : std::get<1>(pair.second)) {
              result.push_back({sequence1->name, {ol1, ol2 }});
            }
          }
        }
      }
      if (there_is_result)
        return result;
      
      int num_pair_w_multi_overlap = 0;
      if (minimizer_result_long_read[0].size() > 1)
        num_pair_w_multi_overlap++;
      if (minimizer_result_long_read[1].size() > 1)
        num_pair_w_multi_overlap++;

      if (num_pair_w_multi_overlap == 1) {
        return {{sequence1->name + "_multiple_long_reads_1_pair", {}}};
      } else if (num_pair_w_multi_overlap == 2) {
        return {{sequence1->name + "_multiple_long_reads_2_pair", {}}};
      }

      
      //std::cerr <<"still have discard" << std::endl;
      return {{"discard", {}}};


    },
    std::ref(minimizer_engine),
    std::cref(sequence1),
    std::cref(sequence2)));
  }

std::vector<std::vector<uint32_t>> Graph::GetComponents(std::vector<std::vector<std::uint32_t>> &matrix) {
  std::cerr << "start get components" << std::endl;
  std::vector<std::vector<uint32_t>> components;
  std::queue<uint32_t> q;
  std::vector<bool> is_visited(contigs.size(), false);
  for (int r = 0; r < static_cast<int>(contigs.size()); r++) {
    std::vector<uint32_t> temp;
    if (!is_visited[r]) {
      q.push(r);
    }
    while (!q.empty()) {
      uint32_t node = q.front();
      q.pop();
      if (is_visited[node]) {
        continue;
      }

      for (int i = 0; i < static_cast<int>(matrix[node].size()); i++) {
        if ((matrix[node][i] != 0 && static_cast<int>(node) != i) && !is_visited[i])
          q.push(i);
      }
      is_visited[node] = true;
      temp.push_back(node);
    }
    if (temp.size()!= 0)
      components.push_back(temp);
  }
  std::cerr << "end get components" << std::endl;
  return components;
}

int Graph::GetNumWindows() {
  int num = 0;
  for (auto contig : contigs) {
    num += contig.second.windows.size();
  }
  return num;
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


void Graph::Store() const {
  std::ofstream os("overlaps.cereal");
  try {
    cereal::BinaryOutputArchive archive(os);
    archive(*this);
  } catch (std::exception&) {
    throw std::logic_error(
        "[tarantula::Construct] error: unable to store archive");
  }
}

void Graph::Load() {
  std::ifstream is("overlaps.cereal");
  try {
    cereal::BinaryInputArchive archive(is);
    archive(*this);
  } catch (std::exception&) {
    throw std::logic_error(
        "[tarantula::Construct] error: unable to load archive");
  }
}

}  // namespace tarantula
