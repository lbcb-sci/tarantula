// Copyright (c) 2021 Cecilia Lee, Robert Vaser

#include "graph.hpp"

#include <cstdlib> 

#include <fstream>
#include <iostream>
#include <queue>
#include <numeric>

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
  std::ofstream myfile, myfile2, myfile3, myfile4, myfile5, myfile6, myfile7;

 // change to false if contig_?.txt is already generated --> will skip to directed force algorithm immediately
  if (false) {
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

  // create graph - initialise contigs & windows for each contig
  std::cerr << "[tarantula::Construct]-----------------CREATE GRAPH-----------------" << std::endl;
  timer.Start();
  CreateGraph(targets);
  std::cerr << "[tarantula::Construct] Graph created, number of nodes: "
            << contigs.size() << " "
            << timer.Stop() << "s"
            << std::endl;

  // map window id such that all windows in all the contigs have unique id
  timer.Start();
  std::size_t bytes = 0;
  window_id_map = GenerateMapWindowID();
  int total_windows = GetNumWindows();
  std::cerr << "[tarantula::Construct] Number of windows = " << total_windows << std::endl;

  // window matrix per contig - key is contig id, value is window matrix
  std::unordered_map<std::uint32_t, ::vector<std::vector<std::uint32_t>>> window_matrix_all_contig;

  // initialise window matrix per contig
  for (std::uint32_t i = 1; i <= window_id_map.size(); i ++) {
    int size;
    if (i == window_id_map.size())
      size = total_windows - window_id_map[i-1];
    else
      size = window_id_map[i] - window_id_map[i-1];
    std::vector<std::vector<std::uint32_t>> temp_matrix(size, std::vector<std::uint32_t>(size, 0));
    window_matrix_all_contig.insert({i-1, temp_matrix});
  }


  std::cerr << "[tarantula::Construct] Total number of sequence: "
            << sequences.size() << std::endl;

  std::cerr << "[tarantula::Construct]-----------------FILTERING-----------------" << std::endl;
  

  // filter pair + in batches of 4gb of sequence
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
              intra_links.reserve(intra_links.capacity() * 1.1);
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
  }
  CalculateAllRestrictionSite(targets, window_size);

  // write graph to text file for directed force algorithm
  std::unordered_map<std::uint32_t, Node>::iterator contigs_iter;
  for (const auto& window_matrix : window_matrix_all_contig) {
    myfile.open("contig_" + std::to_string(window_matrix.first) + ".txt");
    myfile2.open("contig_" + std::to_string(window_matrix.first) + "_scaled_m.txt");
    myfile3.open("contig_" + std::to_string(window_matrix.first) + "_scaled_m_all_nodes.txt");
    myfile5.open("contig_" + std::to_string(window_matrix.first) + "_all_nodes.txt");
    myfile4.open("contig_" + std::to_string(window_matrix.first) + "_restriction_site_analysis.txt");
    contigs_iter = contigs.find(window_matrix.first);
    int node_count = 0;
    bool node_aval = false;
    double weight;
    int count = 0;
    for (std::uint32_t r = 0; r < window_matrix.second.size(); r++) {
      for (std::uint32_t x = 0; x < r; x++) {
        if (x == r-1) {
          //std::cerr << "did this happen at all...." << std::endl;
          node_aval = true;
          int site_count = contigs_iter->second.restriction_sites[r] + contigs_iter->second.restriction_sites[x];
          weight = static_cast<double>(window_matrix.second[r][x])/site_count;
          myfile5 << r << "--" << x << "," << 1 + window_matrix.second[r][x] << "\n";
          myfile3 << r << "--" << x << "," << 1 + weight*100 << "\n";
          if (window_matrix.second[r][x] != 0) {
            myfile2 << r << "--" << x << "," << weight*100 << "\n";
            myfile << r << "--" << x << "," << window_matrix.second[r][x] << "\n";
          }
        } else if (window_matrix.second[r][x] != 0) {
          
          node_aval = true;
          int site_count = contigs_iter->second.restriction_sites[r] + contigs_iter->second.restriction_sites[x];
          weight = static_cast<double>(window_matrix.second[r][x])/site_count;
          myfile << r << "--" << x << "," << window_matrix.second[r][x] << "\n";
          myfile2 << r << "--" << x << "," << weight*100 << "\n"; 
          myfile3 << r << "--" << x << "," << weight*100 << "\n";
          myfile5 << r << "--" << x << "," << window_matrix.second[r][x] << "\n";
          myfile4 << r << "--" << x << "," << window_matrix.second[r][x] << " / (" 
                  << contigs_iter->second.restriction_sites[r] << " + "
                  << contigs_iter->second.restriction_sites[x] << ") = "
                  << window_matrix.second[r][x] << " / " << site_count << " = "
                  << weight << "\n";
          if (weight > 1) {
            count++;
          }
        }
      }
      if (node_aval) {
        node_aval = false;
        node_count++;
      }
    }
    
    /*
    std::cerr << "contig: " << std::to_string(window_matrix.first) 
    << " num nodes: " << window_matrix.second.size() 
    << " node count in graph: " << node_count 
    << " number of edges with edge > 1 = "  << count << std::endl;*/
    myfile.close();
    myfile2.close();
    myfile3.close();
    myfile4.close();
    myfile5.close();
  }
  } else {
    // skip stage
    CreateGraph(targets);
    window_id_map = GenerateMapWindowID();
    std::cerr << "[tarantula::Construct] Number of contigs " << window_id_map.size() << std::endl;
    numRuns = 1;
  }
  std::cerr << "[tarantula::Construct]-----------DIRECTED FORCE ALGORITHM-----------------" << std::endl;
  // directed force algorithm to get position of nodes
  std::vector<std::future<std::string>> futures;
  for (std::uint32_t i = 0; i < window_id_map.size(); i++) {
      // be sure to used std::ref() or std::cref() for references
    futures.emplace_back(thread_pool_->Submit([&] (
      std::vector<int>& window_id_map, 
      int contig_id, 
      int numRuns) -> std::string {
      std::string output_string = "";
      std::string input = "contig_" + std::to_string(contig_id) + ".txt";
      std::string output = "contig_" + std::to_string(contig_id) + "_output.txt";
      std::vector<std::shared_ptr<directedforce::Vertex>> vertices_unmaped;
      std::vector<std::vector<double>> edges;
      std::unordered_map<std::string, std::uint32_t> map_table;
      std::unordered_map<int, std::shared_ptr<directedforce::Vertex>> vertices;
      std::unordered_map<int, std::shared_ptr<directedforce::Vertex>>::iterator iter1;
      std::unordered_map<int, std::shared_ptr<directedforce::Vertex>>::iterator iter2;

      directedforce::GenerateGraphFromDirectedForceAlgorithm(input, output);

      //directedforce::GenerateGraphFromDirectedForceAlgorithm(input, output, vertices_unmaped, edges, map_table);

      input = "contig_" + std::to_string(contig_id) + "_scaled_m.txt";
      output = "contig_" + std::to_string(contig_id) + "_scaled_m_output.txt";

      directedforce::GenerateGraphFromDirectedForceAlgorithm(input, output);

      input = "contig_" + std::to_string(contig_id) + "_scaled_m_all_nodes.txt";
      output = "contig_" + std::to_string(contig_id) + "_scaled_m_all_nodes_output.txt";

      directedforce::GenerateGraphFromDirectedForceAlgorithm(input, output);

      input = "contig_" + std::to_string(contig_id) + "_all_nodes.txt";
      output = "contig_" + std::to_string(contig_id) + "_all_nodes_output.txt";

      directedforce::GenerateGraphFromDirectedForceAlgorithm(input, output);
      
      // no split
      return "";
      
      // find nearest +- 10 nodes for long edges
      if (vertices_unmaped.size() == 0)
        return output_string;

      // map first because the vertices is unmapped to the id
      int max_index = 0;
      int min_index = INT_MAX;
      for (auto& index : map_table) { 
        vertices.insert({std::stoi(index.first), vertices_unmaped[index.second]});
        if (std::stoi(index.first) > max_index) {
          max_index = std::stoi(index.first);
        }
        if (std::stoi(index.first) < min_index) {
          min_index = std::stoi(index.first);
        }
      }
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

      std::sort(conseq_edges_length.begin(), conseq_edges_length.end(),
      [] (const std::tuple<double, int, int>& e1,
          const std::tuple<double, int, int>& e2) -> bool {
        return get<0>(e1) < get<0>(e2);
      });

      int num_conseq_edges = conseq_edges_length.size();
      double median_conseq_edges;
      // calculate median
      if (num_conseq_edges % 2 == 0) {
        // even
        median_conseq_edges = get<0>(conseq_edges_length[num_conseq_edges/2]) +
                              get<0>(conseq_edges_length[num_conseq_edges/2 - 1]);
        median_conseq_edges /= 2;
      } else {
        median_conseq_edges = get<0>(conseq_edges_length[num_conseq_edges/2]);
      }

      
      // contig id, pair: begin. end
      std::unordered_map<int, std::pair<int, int>> contigs_id;
      std::unordered_map<int, std::pair<int, int>>::iterator contigs_itr;
      std::unordered_map<int, int> window_id;
      std::unordered_map<int, int>::iterator window_id_iter;
      std::unordered_map<std::uint32_t , std::uint32_t> nodes_with_long_edges;
      std::unordered_map<std::uint32_t , std::uint32_t>::iterator nodes_iter;
      
      // if edge > *2 median, count the +-10 nodes distance & write to file

      // if edge > *2 median, get all the nodes
      for (auto& conseq_edge : conseq_edges_length) {
        if (get<0>(conseq_edge) >= median_conseq_edges*2) {
          nodes_iter = nodes_with_long_edges.find(get<1>(conseq_edge));
          if (nodes_iter == nodes_with_long_edges.end()) {
            nodes_with_long_edges.insert({get<1>(conseq_edge), 0});
          }
          nodes_iter = nodes_with_long_edges.find(get<2>(conseq_edge));
          if (nodes_iter == nodes_with_long_edges.end()) {
            nodes_with_long_edges.insert({get<2>(conseq_edge), 0});
          }
        }
      }
      output_string += "contig " + to_string(contig_id) + " nodes with long edges: " + to_string(nodes_with_long_edges.size()) + "\n";
      std::unordered_map<uint32_t, double> neighbours;
      std::unordered_map<uint32_t, double>::iterator neighbours_iter;
      std::unordered_map<int, std::shared_ptr<directedforce::Vertex>>::iterator vert_iter;
      std::shared_ptr<directedforce::Vertex> neighbour_vertex, cur_vertex;
      directedforce::MathVector edge;
      std::uint32_t first_nbour, last_nbour;
      uint32_t num_vertices = vertices.size();
      uint32_t multiplier = 10;
      uint32_t node_multiplied;
      while (num_vertices != 0) {
        num_vertices /= 10;
        multiplier*=10;
      }
      for (const auto& node : nodes_with_long_edges) {
        first_nbour = node.first - 10;
        last_nbour = node.first + 10;
        if (last_nbour >= vertices.size()) {
          last_nbour = vertices.size()-1;
        }
        node_multiplied = node.first*multiplier;
        vert_iter = vertices.find(node.first);
        if (vert_iter != vertices.end()) {
          cur_vertex = vert_iter->second;
        }
        
        for (std::uint32_t r = first_nbour; r <= last_nbour; r++) {
          if (r != node.first) {
            // check if its already calculated
            neighbours_iter = neighbours.find(node_multiplied+r);
            if (neighbours_iter == neighbours.end()) {
              // calculate distance between r & node.first
              vert_iter = vertices.find(r);
              if (vert_iter != vertices.end()) {
                neighbour_vertex = vert_iter->second; 
                edge = cur_vertex->pos - neighbour_vertex->pos;
                neighbours.insert({node_multiplied+r, edge.abs()});
              }
            }
          } 
        }
      }
      output_string += "contig " + to_string(contig_id) + " number of neighbours = " + to_string(neighbours.size()) + "\n";
      // get median & deviation for this contig
      std::vector<double> neighbour_distances;
      if (neighbours.size() != 0) {
        for (const auto& n : neighbours)
          neighbour_distances.push_back(n.second);
        
        std::sort(neighbour_distances.begin(), neighbour_distances.end());
        auto median = neighbour_distances[neighbour_distances.size()/2];
        output_string += "contig " + to_string(contig_id) + " Median: " + to_string(median) + "\n";
        double sum = std::accumulate(neighbour_distances.begin(), neighbour_distances.end(), 0.0);
        double mean = sum / neighbour_distances.size();

        double sq_sum = std::inner_product(neighbour_distances.begin(), neighbour_distances.end(), neighbour_distances.begin(), 0.0);
        double stdev = std::sqrt(sq_sum / neighbour_distances.size() - mean * mean);

        output_string += "contig " + to_string(contig_id) + " sd: " + to_string(stdev) + "\n";
      }
      output_string += "\n";
      return output_string;


      // split into smaller bp (25 -> 5)
      for (auto& conseq_edge : conseq_edges_length) {
        if (get<0>(conseq_edge) >= median_conseq_edges*2) {
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

      if (long_edges.size() == 0) {
        std::cerr << "contig " << contig_id  << " have no long edges" << std::endl;
        return output_string;
      }

      std::unordered_map<int, pair<int, bool>> window_split_map;
      std::unordered_map<int, pair<int, bool>>::iterator window_split_map_iter;

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

      std::vector<std::vector<std::uint32_t>> window_split_matrix(starting_node_number, std::vector<std::uint32_t>(starting_node_number, 0));
      for (int r = 0; r < numRuns; r++) {
        Load(r);
        GenerateMatrixAftSplitWindow(contig_id, window_split_map, window_split_matrix);
      }

      
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

      myfile.open("contig_" + std::to_string(contig_id) + "_split.txt");
      for (std::uint32_t r = 0; r < window_split_matrix.size(); r++) {
        for (std::uint32_t x = 0; x < r; x++) {
          if (window_split_matrix[r][x] != 0)
            myfile << r << "--" << x << "," << window_split_matrix[r][x] << "\n";
        } 
      }
      myfile.close();

      directedforce::GenerateGraphFromDirectedForceAlgorithm("contig_" + std::to_string(contig_id) + "_split.txt", "contig_" + std::to_string(contig_id) + "_split_output.txt");
    }, std::ref(window_id_map), i, numRuns));
  }
  std::string output_string;
  myfile.open("median_sd.txt");
  for (auto& it : futures) {
    output_string = it.get();
    myfile << output_string;
  }
  return;
}

uint32_t Graph::KMPSearch(std::string pat, std::string txt) {
  int n = txt.length();
  int m = pat.length();
  uint32_t count = 0;
  if (m == 0)
    return (n == 0);
  int i = 0, j = 0, index_txt = -1;
  while (i < n) {
    if (j < m && txt[i] == pat[j]) {
      i++;
      j++;
    }
    else if (j < m && pat[j] == '?') {
      index_txt = i;
      i++;
      j++;
    }
    else if (index_txt != -1) {
      j = 0;
      i = index_txt;
      index_txt = -1;
    } 
    else {
      j = 0;
      i++;
    }
    if (j == m) {
      count++;
      j = 0;
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

// initialise vector for contig & number of windows in each contig
void Graph::CreateGraph(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>& targets) {
  for (auto const& target : targets) {
    std::cerr << "[tarantula::Construct] contig id: " << target->id << ", len: " << target->inflated_len <<  ", num windows: " << std::ceil((double)target->inflated_len/25000) << std::endl;
    contigs.emplace(target->id, Node(target->id, target->inflated_len, window_size)); 
  }
}

void Graph::CalculateAllRestrictionSite(std::vector<std::unique_ptr<biosoup::NucleicAcid>>& targets, int window_size) {
  std::vector<std::future<void>> futures;
  for (auto const& target : targets) {
    futures.emplace_back(thread_pool_->Submit([&](std::unique_ptr<biosoup::NucleicAcid> const& target, int window_size) 
    ->void{
      std::vector<uint32_t> restriction_sites;
      std::string contig_seq = target->InflateData();
      std::uint32_t start = 0;
      uint32_t gatc, gantc;
      while (start < contig_seq.size()) {
        gatc = KMPSearch("GATC", contig_seq.substr(start, window_size));
        gantc = KMPSearch("GA?TC", contig_seq.substr(start, window_size));
        try {
          if (restriction_sites.capacity() == restriction_sites.size()) 
            restriction_sites.reserve(restriction_sites.capacity() * 1.1);
        }
        catch (...) {
          std::cerr << "RESERVATION OF MEMORY FAILED\n";
        } 
        restriction_sites.push_back(gantc + gatc);
        start += window_size;
      }
      contigs.find(target->id)->second.restriction_sites = restriction_sites;

      std::ofstream myfile;
      //std::cerr << "contig_" + std::to_string(target->id) +"number seg restriction sites: " << restriction_sites.size() << std::endl;
      myfile.open("contig_" + std::to_string(target->id) + "_restriction_sites.txt");
      for (auto const& rs : restriction_sites) {
        myfile << rs << "\n";
      }
      myfile.close();
    }, std::cref(target), window_size));
  }
  for (const auto& it : futures) {
    it.wait();
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

void Graph::CalcualteInterChromosomeLinks() {
  std::unordered_map<std::uint32_t, Node>::iterator contig_iter;
  int window_index_begin;
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

        // comment if do not want multiple long reads 
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
        
        // uncomment if you do not want multiple long read
        //return {};

        std::cerr << "ERROR IN FILTERING" << std::endl;
        return {};
      },
    std::ref(minimizer_engine),
    std::cref(sequence1),
    std::cref(sequence2)));
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

// store overlap, i is the batch num
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

// load overlap, i is the batch num
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
