// Copyright (c) 2021 Cecilia Lee, Robert Vaser

#include "graph.hpp"

#include <cstdlib> 

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
  std::vector<int> window_id_map;
  int numRuns = 0;
  std::ofstream myfile, myfile2;
 
  // if (false) {
  // paramters for RAM
  uint32_t k = 21, w = 11, bandwidth = 100, chain = 2, matches = 25, gap = 100;
  double frequency = 0.0001;

  // statistics
  int num_pair = 0;
  int num_filtered_links = 0;

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

  std::vector<std::future<std::vector<Link>>> futures;
  biosoup::Timer timer;

  timer.Start();
  // then create pilo-o-gram per contig
  CreateGraph(targets);
  std::cerr << "[tarantula::Construct] Graph created, number of nodes: "
            << contigs.size() << " "
            << timer.Stop() << "s"
            << std::endl;


  // filter pair + less than 4GB
  timer.Start();
  std::size_t bytes = 0;
  window_id_map = GenerateMapWindowID();
  int total_windows = GetNumWindows();
  std::cerr << "[tarantula::Construct] Number of windows = " << total_windows << std::endl;
  //std::vector<std::vector<std::uint32_t>> window_matrix(total_windows, std::vector<std::uint32_t>(total_windows, 0));

  // window matrix per contig
  std::unordered_map<std::uint32_t, ::vector<std::vector<std::uint32_t>>> window_matrix_all_contig;

  // initialise window matrix per contig
  for (int i = 1; i <= window_id_map.size(); i ++) {
    int size;
    if (i == window_id_map.size())
      size = total_windows - window_id_map[i-1];
    else
      size = window_id_map[i] - window_id_map[i-1];
    std::vector<std::vector<std::uint32_t>> temp_matrix(size, std::vector<std::uint32_t>(size, 0));
    window_matrix_all_contig.insert({i-1,temp_matrix});
  }


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
        for (auto& result : result_vector) {
          try {
            if (intra_links.capacity() == intra_links.size()) 
              intra_links.reserve(intra_links.capacity() * 1.5);
          }
          catch (...) {
            std::cerr << "RESERVATION OF MEMORY FAILED\n";
          } 
          if (isIntraLink(result)) {
            intra_links.emplace_back(result);
          } else {
            inter_links.emplace_back(result);
          }
        }        
      }
      num_filtered_links += intra_links.size();
      std::cerr << "[tarantula::Construct] Number of good inter links: "
          << inter_links.size() << std::endl;
      std::cerr << "[tarantula::Construct] Number of good intra links: "
          << intra_links.size() << " "
          << timer.Stop() << "s"
          << std::endl;

      // interwindow links
      GenerateMatrixWindowIntraLinks(window_matrix_all_contig);
      std::cerr << "[tarantula::Construct] Generate matrix window intra links" << std::endl;
      // save overlaps
      Store(numRuns);
      // discard read pair
      intra_links.clear();
      futures.clear();
      bytes = 0;
      numRuns++;
    }
  }
  // inter.close();
  // intra.close();
  std::cerr << "[tarantula::Construct]-----------FILTERING STATISTICS-----------------" << std::endl;
  std::cerr << "[tarantula::Construct] Total number of reads: "
            << sequences.size() << std::endl;
  std::cerr << "[tarantula::Construct] Number of reads that are not a pair: "
            << sequences.size() - num_pair*2 << std::endl;
  std::cerr << "[tarantula::Construct] Number of intra links that are sucessfully filtered: "
            << num_filtered_links << std::endl;
  std::cerr << "[tarantula::Construct] Number of inter links that are sucessfully filtered: "
            << inter_links.size() << std::endl;
  std::cerr << "[tarantula::Construct]-------------------------------------------------" << std::endl;

  CalcualteInterChromosomeLinks();
  uint32_t sum_interchromosome_links = 0, sum_intrachromosome_links = 0;
  int num_chromosomes = targets.size();
  std::vector<std::vector<std::uint32_t>> matrix(num_chromosomes, std::vector<std::uint32_t>(num_chromosomes, 0));
  std::unordered_map<std::uint32_t, Node>::iterator it;
  GenerateMatrix(matrix);

  for (const auto& window_matrix : window_matrix_all_contig) {
    string file_name = "window_matrix_contig_" + std::to_string(window_matrix.first) + ".csv";
    myfile.open(file_name);
  
    string output = ",";
    for (int i = 0; i < static_cast<int>(window_matrix.second.size()); i++) {
      output += std::to_string(i) + ",";
    }
    output += "Total\n";
    myfile << output;

    for (int i = 0; i < static_cast<int>(window_matrix.second.size()); i++) {
      auto sum = 0;
      output ="";
      output += std::to_string(i);
      for (int r = 0; r < static_cast<int>(window_matrix.second[i].size()); r++) {
        output += "," + std::to_string(window_matrix.second[i][r]);
        sum += window_matrix.second[i][r];
      }
      output += "," + std::to_string(sum)+"\n";
      myfile << output;
    }
    myfile.close(); 
  }


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

  std::unordered_map<std::uint32_t, Node>::iterator contigs_iter;
  for (const auto& window_matrix : window_matrix_all_contig) {
    myfile.open("contig_" + std::to_string(window_matrix.first) + ".txt");
    myfile2.open("contig_" + std::to_string(window_matrix.first) + "_scaled.txt");
    contigs_iter = contigs.find(window_matrix.first);
    for (int r = 0; r < window_matrix.second.size(); r++) {
      for (int x = 0; x < r; x++) {
        if (window_matrix.second[r][x] != 0) {
          int site_count = contigs_iter->second.restriction_sites[r] + contigs_iter->second.restriction_sites[x];
          myfile << r << "--" << x << "," << window_matrix.second[r][x] << "\n";
          myfile2 << r << "--" << x << "," << (double)window_matrix.second[r][x]/site_count << "\n";
        }       
      }
    }
    myfile.close();
    myfile2.close();
  }
  /*} else {
    // skip stage
    CreateGraph(targets);
    window_id_map = GenerateMapWindowID();
    std:cerr << "Number of contigs " << window_id_map.size() << std::endl;
    numRuns = 1;
  }*/
  // directed force
  //window_id_map.size()
  //std::cerr << "whats happeniing here? " << std::endl;
  std::vector<std::future<void>> void_futures;
  for (std::uint32_t i = 0; i < window_id_map.size(); i++) {
      // be sure to used std::ref() or std::cref() for references
    void_futures.emplace_back(thread_pool_->Submit([&] (
      std::vector<int>& window_id_map, 
      int contig_id, 
      int numRuns) -> void {
      std::ofstream myfile;
      std::string input = "contig_" + std::to_string(contig_id) + "_scaled.txt";
      std::string output = "contig_" + std::to_string(contig_id) + "_scaled_output.txt";
      std::vector<std::shared_ptr<directedforce::Vertex>> vertices_unmaped;
      std::vector<std::vector<double>> edges;
      std::unordered_map<std::string, int> map_table;
      std::unordered_map<int, std::shared_ptr<directedforce::Vertex>> vertices;
      std::unordered_map<int, std::shared_ptr<directedforce::Vertex>>::iterator iter1;
      std::unordered_map<int, std::shared_ptr<directedforce::Vertex>>::iterator iter2;

      directedforce::GenerateGraphFromDirectedForceAlgorithm(input, output, vertices_unmaped, edges, map_table);

      // no split
      return;
      
      if (vertices_unmaped.size() == 0)
        return;
      // map first
      int max_index = 0;
      int min_index = INT_MAX;
      //std::cerr << "directed force map table size: " << map_table.size() << endl;
      for (auto& index : map_table) { 
        vertices.insert({std::stoi(index.first),vertices_unmaped[index.second]});
        if (std::stoi(index.first) > max_index) {
          max_index = std::stoi(index.first);
        }
        if (std::stoi(index.first) < min_index) {
          min_index = std::stoi(index.first);
        }
      }
      //std::cerr << "min index: " << min_index << std::endl;
      //std::cerr << "max index: " << max_index << std::endl; 
      std::vector<std::tuple<double, int, int>> conseq_edges_length;
      std::vector<std::tuple<double, int, int>> long_edges;
      for (int r = min_index; r < max_index; r++) {
        iter1 = vertices.find(r);
        iter2 = vertices.find(r+1);
        if (iter1 != vertices.end() && iter2 != vertices.end()) {
          // edge exist
          // calculate edge length
          directedforce::MathVector edge = iter1->second->pos - iter2->second->pos;
          conseq_edges_length.emplace_back(std::make_tuple(edge.abs(), r, r+1));
        }
      }

      //std::cerr << "crash where?" << std::endl;
      
      std::sort(conseq_edges_length.begin(), conseq_edges_length.end(),
      [] (const std::tuple<double, int, int>& e1,
          const std::tuple<double, int, int>& e2) -> bool {
        return get<0>(e1) < get<0>(e2);
      });

      //std::cerr << "sorting?" << std::endl;
      int num_conseq_edges = conseq_edges_length.size();
      //std::cerr << "size: " << num_conseq_edges << std::endl;
      double median_conseq_edges;
      // calculate median
      if (num_conseq_edges % 2 == 0) {
        // even
        //std::cerr << "size/2 " << get<0>(conseq_edges_length[num_conseq_edges/2]) << std::endl;
        median_conseq_edges = get<0>(conseq_edges_length[num_conseq_edges/2]) +
                              get<0>(conseq_edges_length[num_conseq_edges/2 - 1]);
        median_conseq_edges /= 2;
      } else {
        //std::cerr << "size/2 " << get<0>(conseq_edges_length[num_conseq_edges/2]) << std::endl;
        median_conseq_edges = get<0>(conseq_edges_length[num_conseq_edges/2]);
      }
      //std::cerr << "calculate median" << std::endl;
      // contig id, pair: begin. end
      std::unordered_map<int, std::pair<int, int>> contigs_id;
      std::unordered_map<int, std::pair<int, int>>::iterator contigs_itr;
      // if edge > 2x median

      std::unordered_map<int, int> window_id;
      std::unordered_map<int, int>::iterator window_id_iter;
      // change to if edge > 2x median --> window id

      for (auto& conseq_edge : conseq_edges_length) {
        if (get<0>(conseq_edge) >= median_conseq_edges) {
          long_edges.emplace_back(conseq_edge);
          window_id_iter = window_id.find(get<1>(conseq_edge));
          if (window_id_iter == window_id.end()) {
            window_id.insert({get<1>(conseq_edge), 0});
          }
          window_id_iter = window_id.find(get<2>(conseq_edge));
          if (window_id_iter == window_id.end()) {
            window_id.insert({get<2>(conseq_edge), 0});
          }
        }
      }

      for (auto& w : window_id) {
        //std::cerr << w.first << ", " << w.second << std::endl;
      }

     // std::cerr << "get long edge" << std::endl;

      if (long_edges.size() == 0) {
        std::cerr << "contig " << contig_id  << " have no long edges" << std::endl;
        return;
      }


      //std::cerr << "crash where?" << std::endl;

      std::unordered_map<int, pair<int, bool>> window_split_map;
      std::unordered_map<int, pair<int, bool>>::iterator window_split_map_iter;

      int num_windows_this_contig = window_id_map[contig_id+1] - window_id_map[contig_id];
      int start = 0;
      int end = window_id_map[contig_id+1] - window_id_map[contig_id];
      //std::cerr << "start: " << start <<" , end: " << end << std::endl; 
      int starting_node_number = 0;
      int window_id_in_each_contig = 0;
      for (int r = start; r < end; r++) {
        window_id_iter = window_id.find(r);
        if (window_id_iter != window_id.end()) {
          window_split_map.insert({window_id_in_each_contig, {starting_node_number, true}});
          starting_node_number += 5;
        } else {
          window_split_map.insert({window_id_in_each_contig, {starting_node_number, false}});
          starting_node_number += 1;
        }
        window_id_in_each_contig++;
      }

      //std::cerr << "matrix size: " << starting_node_number << std::endl;

      // print window split map
      for (auto& window_split : window_split_map) {
        //std::cerr << "window: " << window_split.first << ", map: " << window_split.second.first << " ," << window_split.second.second <<  std::endl;
      }

      //std::cerr << "crash where?" << std::endl;

      //std::cerr << "number of windows that require splitting: " << window_id.size() << std::endl;
      //std::cerr << "original matrix size: " << window_matrix.size() << std::endl; 
      //std::cerr << "current matrix size aft splitting: " << starting_node_number << std::endl; 
      std::vector<std::vector<std::uint32_t>> window_split_matrix(starting_node_number, std::vector<std::uint32_t>(starting_node_number, 0));
      for (int r = 0; r < numRuns; r++) {
        Load(r);
        //std::cerr << "load sucess? " << intra_links.size() << std::endl;
        GenerateMatrixAftSplitWindow(contig_id, window_split_map, window_split_matrix);
      }

      //std::cerr << "after generate matrix, split matrix size: " << window_split_matrix.size() << std::endl;

        // window matrix csv
      
      myfile.open("contig_" + std::to_string(contig_id) + "_window_matrix.csv");
      
      output = ",";
      for (int i = 0; i < static_cast<int>(window_split_matrix.size()); i++) {
        output += std::to_string(i) + ",";
      }
      output += "Total\n";
      myfile << output;

      for (int i = 0; i < static_cast<int>(window_split_matrix.size()); i++) {
        auto sum = 0;
        output ="";
        output += std::to_string(i);
        for (int r = 0; r < static_cast<int>(window_split_matrix[i].size()); r++) {
          output += "," + std::to_string(window_split_matrix[i][r]);
          sum += window_split_matrix[i][r];
        }
        output += "," + std::to_string(sum)+"\n";
        myfile << output;
      }
      myfile.close();

      // split into window into 1kbp
     // std::cerr << "csv done" << std::endl;

      myfile.open("contig_" + std::to_string(contig_id) + "_split.txt");
      for (int r = 0; r < window_split_matrix.size(); r++) {
        for (int x = 0; x < r; x++) {
          if (window_split_matrix[r][x] != 0)
            myfile << r << "--" << x << "," << window_split_matrix[r][x] << "\n";
        } 
      }
      myfile.close();

      //std::cerr << "txt file done" << std::endl;

      directedforce::GenerateGraphFromDirectedForceAlgorithm("contig_" + std::to_string(contig_id) + "_split.txt", "contig_" + std::to_string(contig_id) + "_split_output.txt");
    }, std::ref(window_id_map), i, numRuns));
  }
  for (const auto& it : void_futures) {
    it.wait();
  }
  return;
}

/*
uint32_t Graph::CalculateRestrictionSites(std::unique_ptr<biosoup::NucleicAcid> const& target) {
  std::string contig_seq = target->InflateData();
  uint32_t gatc = KMPSearch("GATC", contig_seq);
  uint32_t gantc = KMPSearch("GANTC", contig_seq);
  return gantc + gatc;
}*/

std::vector<uint32_t> Graph::CalculateRestrictionSites(std::unique_ptr<biosoup::NucleicAcid> const& target, int window_size) {
  std::vector<uint32_t> restriction_sites;
  std::string contig_seq = target->InflateData();
  std::cerr << "contig seq size: " << contig_seq.size() << std::endl;
  int size = ceil(contig_seq.size()/window_size);
  int start = 0, end = window_size;
  uint32_t gatc, gantc;
  while (start < contig_seq.size()) {
    if (end > contig_seq.size()) {
      end = contig_seq.size();
    }
    gatc = KMPSearch("GATC", contig_seq.substr(start, end));
    gantc = KMPSearch("GANTC", contig_seq.substr(start, end));
    restriction_sites.push_back(gantc + gatc);
    start += window_size;
    end += window_size;
    //std::cerr << "start: " << start << " ,end: " << end << std::endl;
  }
  
  return restriction_sites;
}

uint32_t Graph::KMPSearch(std::string pat, std::string txt)
{
  int M = pat.size();
  int N = txt.size();
  uint32_t count = 0;
  // Preprocess the pattern (calculate lps[] array)
  int i = 0; // index for txt[]
  int j = 0; // index for pat[]
  while (i < N) {
    if (pat[j] == txt[i]) {
      j++;
      i++;
    }
    if (j == M) {
      count++;
      j = 0;
    }
    // mismatch after j matches
    else if (i < N && pat[j] != txt[i]) {
        // Do not match lps[0..lps[j-1]] characters,
        // they will match anyway
      if (j != 0)
          j = 0;
      else
          i = i + 1;
    }
  }
  return count;
}


bool Graph::isIntraLink(Link &link) {
  if (link.rhs_id == link.lhs_id) 
    return true;
  else 
    return false;
}

int Graph::FindContigID(int window_id, std::vector<int>& window_id_map, int *begin, int *end) {
  for (int i = 0; i < window_id_map.size(); i++) {
    if (window_id < window_id_map[i]) {
      *begin = (window_id - window_id_map[i-1])*window_size;
      *end = (window_id - window_id_map[i-1])*window_size + window_size;
      return i-1;
    }
      
  }
}

// technically this also can be multi-thread - just remove the for loop
void Graph::CreateGraph(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>& targets) {
  for (auto const& target : targets) {
    std::vector<uint32_t> restriction_sites = CalculateRestrictionSites(target, window_size);
    //std::cerr << "contig : " << target->id << std::endl;
    /*
    for (auto const& rs : restriction_sites) {
      std::cerr << rs << " ";
    }
    std::cerr << "\n";*/
    contigs.emplace(target->id, Node(target->id, target->inflated_len, window_size, restriction_sites));    
  }
}

void Graph::GenerateMatrix(
  std::vector<std::vector<std::uint32_t>> &matrix) {
  for (const auto& inter_link : inter_links) {
    auto id_1 = inter_link.rhs_id;
    auto id_2 = inter_link.lhs_id;
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
  std::unordered_map<std::uint32_t, ::vector<std::vector<std::uint32_t>>> &matrix) {
    std::unordered_map<std::uint32_t, ::vector<std::vector<std::uint32_t>>>::iterator matrix_iter;
    for (const auto& intra_link : intra_links) {
      matrix_iter = matrix.find(intra_link.rhs_id);
      if (matrix_iter == matrix.end())
        std::cerr << "ERROR IN GENERATING INTRA LINKS WINDOW MATRIX" << std::endl;
      int window_id_1= intra_link.rhs_begin/window_size;
      int window_id_2 = intra_link.lhs_begin/window_size;
      matrix_iter->second[window_id_1][window_id_2]+=1;
      matrix_iter->second[window_id_2][window_id_1]+=1;
    }
  }

/*
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
}*/

// windows split in order
void Graph::GenerateMatrixAftSplitWindow(
  int contig_id,
  std::unordered_map<int, pair<int, bool>>& window_split_map,
  std::vector<std::vector<std::uint32_t>> &matrix) {
    std::unordered_map<int, pair<int, bool>>::iterator window_split_map_iter;

    int window_id_1 = -1, window_id_2 = -1;
    for (const auto& intra_link : intra_links) {
      if (intra_link.rhs_id == contig_id) {
        // find window 
        int window_index_begin_1 = intra_link.rhs_begin/window_size;
        int window_index_begin_2 = intra_link.lhs_begin/window_size;

        window_split_map_iter = window_split_map.find(window_index_begin_1);
        if (window_split_map_iter != window_split_map.end()) {
          if (window_split_map_iter->second.second) {
            // split window
            int split_window_id = (intra_link.rhs_begin - window_index_begin_1*window_size)/5000;
            //std::cerr << "split window id 1: " << split_window_id << std::endl;
            window_id_1 = window_split_map_iter->second.first + split_window_id;
          } else {
            window_id_1 = window_split_map_iter->second.first;
          }
        } else {
          std::cerr << "ERROR IN MAPPING" << std::endl;
        }

        window_split_map_iter = window_split_map.find(window_index_begin_2);
        if (window_split_map_iter != window_split_map.end()) {
          if (window_split_map_iter->second.second) {
            // split window
            int split_window_id = (intra_link.lhs_begin - window_index_begin_2*window_size)/5000;
            //std::cerr << "split window id 2: " << split_window_id << std::endl;
            window_id_2 = window_split_map_iter->second.first + split_window_id;
          } else {
            window_id_2 = window_split_map_iter->second.first;
          }
        } else {
          std::cerr << "ERROR IN MAPPING" << std::endl;
        }
        //std::cerr << "window id 1: " << window_id_1 << " ,window_id_2: " << window_id_2 << std::endl;
        if (window_id_1 != -1 && window_id_2 != -1) {
          matrix[window_id_1][window_id_2] += 1;
          matrix[window_id_2][window_id_1] += 1;
        } else {
          std::cerr << "window index not assigned" << std::endl;
        }
        

        window_id_1 = -1;
        window_id_2 = -1;
        
      }
    }
  }

/*
void Graph::GenerateMatrixAftSplitWindow(
  std::vector<int>& window_id_map,
  std::unordered_map<int, int>& window_split_map,
  std::vector<std::vector<std::uint32_t>> &matrix,
  std::unordered_map<std::string, std::vector<biosoup::Overlap>>& read_pairs) {
    std::unordered_map<int, int>::iterator window_split_map_iter1, window_split_map_iter2;
    for (const auto& rp : read_pairs) {
      int window_index_begin_1 = rp.second[0].rhs_begin/window_size;
      int window_index_begin_2 = rp.second[1].rhs_begin/window_size;
      int window = window_id_map[rp.second[1].rhs_id];
      int window_id_1 = window + window_index_begin_1;
      int window_id_2 = window + window_index_begin_2;
      //std::cerr << "window id 1: " << window_id_1 << ", window id 2: " << window_id_2 << std::endl; 

      window_split_map_iter1 = window_split_map.find(window_id_1);
      if (window_split_map_iter1 != window_split_map.end()) {
        //std::cerr << "+ " << (rp.second[0].rhs_begin - window_id_1*window_size)/1000 << std::endl;
        int window_split_begin_1 = (rp.second[0].rhs_begin - window_index_begin_1*window_size)/1000 + window_split_map_iter1->second;
        //std::cerr << "window split 1: " << window_split_begin_1 << std::endl;
        window_split_map_iter2 = window_split_map.find(window_id_2);
        if (window_split_map_iter2 != window_split_map.end()) {
          //std::cerr << "+ " << (rp.second[1].rhs_begin - window_index_begin_2*window_size)/1000 << std::endl;
          int window_split_begin_2 = (rp.second[1].rhs_begin - window_index_begin_2*window_size)/1000+ window_split_map_iter2->second;
          matrix[window_split_begin_1][window_split_begin_2] += 1;
          matrix[window_split_begin_2][window_split_begin_1] += 1;
          //std::cerr << "window split 2: " << window_split_begin_2 << std::endl;
        } else {
          matrix[window_split_begin_1][window_id_2] += 1;
          matrix[window_id_2][window_split_begin_1] += 1;
        }
        // if no window split, then just use window_id
      } else {
        window_split_map_iter2 = window_split_map.find(window_id_2);
        if (window_split_map_iter2 != window_split_map.end()) {
          //std::cerr << "+ " << (rp.second[1].rhs_begin - window_index_begin_2*window_size)/1000 << std::endl;
          int window_split_begin_2 = (rp.second[1].rhs_begin - window_index_begin_2*window_size)/1000 + window_split_map_iter2->second;
          //std::cerr << "window split 2: " << window_split_begin_2 << std::endl;
          matrix[window_id_1][window_split_begin_2] += 1;
          matrix[window_split_begin_2][window_id_1] += 1;
        } else {
          matrix[window_id_1][window_id_2] += 1;
          matrix[window_id_2][window_id_1] += 1;
        }
      }
    }
  }*/

/*
void Graph::GenerateMatrixWindowIntraLinks(
  std::vector<int>& window_id_map,
  std::vector<std::vector<std::uint32_t>> &matrix,
  std::unordered_map<std::string, std::vector<biosoup::Overlap>>& read_pairs,
  int window_size) {
  for (const auto& rp : read_pairs) {
    int window_index_begin_1 = rp.second[0].rhs_begin/window_size;
    int window_index_begin_2 = rp.second[1].rhs_begin/window_size;
    int window = window_id_map[rp.second[1].rhs_id];
    int window_id_1 = window + window_index_begin_1;
    int window_id_2 = window + window_index_begin_2;
    matrix[window_id_1][window_id_2]+=1;
    matrix[window_id_2][window_id_1]+=1;
  }
} */

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

// for interlinks
void Graph::GenerateMatrixWindow(
  std::vector<int>& window_id_map,
  std::vector<std::vector<std::uint32_t>> &matrix) {
  std::unordered_map<std::uint32_t, Node>::iterator found;
  int extra = 0;

  for (const auto& inter_link : inter_links) {
    int window_index_begin_1 = inter_link.rhs_begin/window_size;
    int window_index_begin_2 = inter_link.lhs_begin/window_size;
    int id_1 = inter_link.rhs_id;
    int id_2 = inter_link.lhs_id;

    int window_id_1_begin = window_index_begin_1+window_id_map[id_1];
    int window_id_2_begin = window_index_begin_2+window_id_map[id_2];

    matrix[window_id_1_begin][window_id_2_begin] += 1;
    matrix[window_id_2_begin][window_id_1_begin] += 1;
  }
}

void Graph::CalcualteInterChromosomeLinks() {
  std::unordered_map<std::uint32_t, Node>::iterator contig_iter;
  int window_index_begin;
  int extra = 0;
  for (const auto& inter_link : inter_links) {
    window_index_begin = inter_link.rhs_begin/window_size;
  
    contig_iter = contigs.find(inter_link.rhs_id);
    if (contig_iter == contigs.end()) {
      std::cerr << "ERROR contig not found" << std::endl;
    } else {
      contig_iter->second.windows[window_index_begin].interchromosome_links++;
      contig_iter->second.interchromosome_links++;
    }
    window_index_begin = inter_link.lhs_begin/window_size;
    contig_iter = contigs.find(inter_link.lhs_id);
    if (contig_iter == contigs.end()) {
      std::cerr << "ERROR contig not found" << std::endl;
    } else {
      contig_iter->second.windows[window_index_begin].interchromosome_links++;
      contig_iter->second.interchromosome_links++;
    }
  }
}

/*
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
}*/

/*
void Graph::FillPileogram() {
  std::unordered_map<std::uint32_t, Node>::iterator contig_iter;
  std::pair<std::uint32_t, std::uint32_t> overlap;
  std::uint32_t average_overlap = 0;
  std::unordered_map<std::uint32_t, std::vector<std::pair<std::uint32_t, std::uint32_t>>> overlap_map;
  std::unordered_map<std::uint32_t, std::vector<std::pair<std::uint32_t, std::uint32_t>>>::iterator overlap_map_iter;

  for (const auto& intra_link : intra_links) {
    contig_iter = contigs.find(intra_link.rhs_id);
    if (contig_iter == contigs.end()) {
      std::cerr << "ERROR CONTIG NOT FOUND" << std::endl;
    } else {
      // update intrachromosome links
      contig_iter->second.intrachromosome_links++;
      
      // overlap

      auto length = intra_link.rhs_begin - intra_link.lhs_begin;
      if (length < -1) {
        length *= -1;
      }
      average_overlap += length;

      // throw into overlap map
      overlap_map_iter = overlap_map.find(intra_link.rhs_id);
      if (overlap_map_iter == overlap_map.end()) {
        // not found, create map
        std::vector<std::pair<std::uint32_t, std::uint32_t>> temp;
        temp.emplace_back(overlap);
        overlap_map.insert(std::make_pair(intra_link.rhs_id, temp));
      } else {
        overlap_map_iter->second.emplace_back(overlap);
      }
    }
  }
  for (overlap_map_iter = overlap_map.begin(); overlap_map_iter != overlap_map.end(); overlap_map_iter++) {
    contig_iter = contigs.find(overlap_map_iter->first);
    contig_iter->second.pileogram.AddLayer(overlap_map_iter->second);
  }
  average_overlap /= intra_links.size();
  
  std::cerr << "[tarantula::Construct] Stats: " << std::endl;
  std::cerr << "[tarantula::Construct] average overlap: " << average_overlap << std::endl;
}*/

/*
void Graph::FillPileogram() {  // NOLINT
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
} */

/*
std::pair<std::uint32_t, std::uint32_t> Graph::GetOverlap(const Link& link) {
  return std::make_pair(
    std::min(link.rhs_begin, link.lhs_begin),
    std::max(link.rhs_end,   link.lhs_end));
}*/

/*
std::pair<std::uint32_t, std::uint32_t> Graph::GetOverlap(
    biosoup::Overlap ol1,
    biosoup::Overlap ol2) {
  return std::make_pair(
      std::min(ol1.rhs_begin, ol2.rhs_begin),
      std::max(ol1.rhs_end,   ol2.rhs_end));
}*/


void Graph::Process(
  std::vector<std::future<std::vector<Link>>>& futures,
  ram::MinimizerEngine& minimizer_engine,
  std::unique_ptr<biosoup::NucleicAcid>& sequence1,
  std::unique_ptr<biosoup::NucleicAcid>& sequence2) {
    futures.emplace_back(thread_pool_->Submit([&] (
      ram::MinimizerEngine& minimizer_engine,
      const std::unique_ptr<biosoup::NucleicAcid>& sequence1,
      const std::unique_ptr<biosoup::NucleicAcid>& sequence2)
      -> std::vector<Link>{
        std::vector<std::vector<biosoup::Overlap>> minimizer_result, minimizer_result_long_read;
        minimizer_result.emplace_back(minimizer_engine.Map(sequence1, false, false, false));
        minimizer_result.emplace_back(minimizer_engine.Map(sequence2, false, false, false));
        auto long_read_len = sequence1->inflated_len*0.75;

        // no overlap for 1/both --> filter out
        if (minimizer_result[0].size() < 1 || minimizer_result[1].size() < 1) {
          return {};
        }
        
        // keep reads > 0.75*read_length
        for (auto& minimizer : minimizer_result) {
          std::vector<biosoup::Overlap> filtered_overlaps;
          for (auto& overlap: minimizer) {
            if (overlap.rhs_end - overlap.rhs_begin > long_read_len) {
              filtered_overlaps.push_back(overlap);
            }
          }
          minimizer_result_long_read.push_back(filtered_overlaps);
        }

        // no long reads for 1 / both --> discard
        if (minimizer_result_long_read[0].size() == 0 || minimizer_result_long_read[1].size() == 0) {
          return {};
        }

        // straight forward interchromosome / intrachromosome
        if (minimizer_result_long_read[0].size() == 1 && minimizer_result_long_read[1].size() == 1) {
          std::vector<Link> result;
          Link link = {
            static_cast<uint16_t>(minimizer_result_long_read[0][0].rhs_id),
            static_cast<uint16_t>(minimizer_result_long_read[1][0].rhs_id),
            minimizer_result_long_read[0][0].rhs_begin, 
            minimizer_result_long_read[1][0].rhs_begin,
            };
          
          result.emplace_back(link);
          return result;
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

        
        std::vector<Link> result;
        bool there_is_result = false;
        for (auto const& pair : intra_pairs) {
          if (std::get<0>(pair.second).size() != 0 && std::get<1>(pair.second).size() != 0) {
            there_is_result = true;
            for (auto const& ol1 : std::get<0>(pair.second)) {
              for (auto const& ol2 : std::get<1>(pair.second)) {
                Link temp_link = {
                  static_cast<uint16_t>(ol1.rhs_id),
                  static_cast<uint16_t>(ol2.rhs_id),
                  ol1.rhs_begin,
                  ol2.rhs_begin};

                result.push_back(temp_link);
              }
            }
          }
        }
        if (there_is_result)
          return result;
        else 
          return {};

        
        std::cerr << "ERROR IN FILTERING" << std::endl;
        return {};


      },
    std::ref(minimizer_engine),
    std::cref(sequence1),
    std::cref(sequence2)));
}

/*
// map hi-c reads to contig & filter
void Graph::Process(
  std::vector<std::future<std::vector<std::pair<std::string, std::vector<biosoup::Overlap>>>>>& futures,
  ram::MinimizerEngine& minimizer_engine,
  std::unique_ptr<biosoup::NucleicAcid>& sequence1,
  std::unique_ptr<biosoup::NucleicAcid>& sequence2) {

  futures.emplace_back(thread_pool_->Submit([&] (
    ram::MinimizerEngine& minimizer_engine,
    const std::unique_ptr<biosoup::NucleicAcid>& sequence1,
    const std::unique_ptr<biosoup::NucleicAcid>& sequence2)
    -> std::vector<std::pair<std::string, std::vector<biosoup::Overlap>>>{
      //std::vector<std::uint32_t> filtered1, filtered2;
      std::vector<std::vector<biosoup::Overlap>> minimizer_result, minimizer_result_long_read;
      minimizer_result.emplace_back(minimizer_engine.Map(sequence1, false, false, false));
      minimizer_result.emplace_back(minimizer_engine.Map(sequence2, false, false, false));
      auto long_read_len = sequence1->inflated_len*0.75;

      // no overlap for 1 / both
      if (minimizer_result[0].size() < 1 || minimizer_result[1].size() < 1) {
        // no overlap, discard
        if (minimizer_result[0].size() > 0 || minimizer_result[1].size() > 0) {
          // minimizer_result.clear();
          /*
          if (filtered1.size() > 4 || filtered2.size() > 4)
            return {{"1_overlap_mt_4", {}}};
          if (filtered1.size() > 0 || filtered2.size() > 0)
            return {{"1_overlap_mt_0", {}}};*/

          //minimizer_result.clear();
/*          return {{"1_overlap", {}}};
        } else {
          // minimizer_result.clear();
          /*
          if (filtered1.size() > 4 || filtered2.size() > 4)
            return {{"empty_mt_4", {}}};
          if (filtered1.size() > 0 || filtered2.size() > 0)
            return {{"empty_mt_0", {}}};*/
          //minimizer_result.clear();
/*          return {{"empty", {}}};
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
        /*
        if (filtered1.size() > 4 || filtered2.size() > 4)
          return {{sequence1->name + "__multiple_short_reads_mt_4", {}}};
        if (filtered1.size() > 0 || filtered2.size() > 0)
          return {{sequence1->name + "_multiple_short_reads_mt_0", {}}};*/
/*        return {{sequence1->name + "_multiple_short_reads", {}}};
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
  }*/

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


void Graph::Store(int i) const {
  std::ofstream os("overlaps_" + to_string(i) + ".cereal");
  try {
    cereal::BinaryOutputArchive archive(os);
    archive(*this);
  } catch (std::exception&) {
    throw std::logic_error(
        "[tarantula::Construct] error: unable to store archive");
  }
}

void Graph::Load(int i) {
  std::ifstream is("overlaps_" + to_string(i) + ".cereal");
  intra_links.clear();
  try {
    cereal::BinaryInputArchive archive(is);
    archive(*this);
  } catch (std::exception&) {
    throw std::logic_error(
        "[tarantula::Construct] error: unable to load archive");
  }
}

}  // namespace tarantula
