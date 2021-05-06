// Copyright (c) 2020 Robert Vaser

#ifndef TARANTULA_GRAPH_HPP_
#define TARANTULA_GRAPH_HPP_

#include <math.h>

#include <memory>
#include <unordered_map>
#include <tuple>
#include <algorithm>
#include <utility>
#include <string>
#include <vector>

#include "biosoup/overlap.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "thread_pool/thread_pool.hpp"
#include "ram/minimizer_engine.hpp"
#include "cereal/access.hpp"
#include "cereal/types/unordered_map.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/types/string.hpp"
#include "cereal/archives/binary.hpp"

#include "pileogram.hpp"
#include "algorithm.h"

namespace biosoup {
  template<class Archive>
  void serialize(Archive & archive, Overlap& m) {
    archive(m.alignment, m.lhs_begin, m.lhs_end, m.lhs_id, m.rhs_begin, m.rhs_end, m.rhs_id, m.score, m.strand);
  }
}

namespace tarantula {

struct Window{
  uint32_t id;
  uint32_t interchromosome_links;
  uint32_t intrachromosome_links;
  Window(uint32_t id)
    : id(id),
      interchromosome_links(0),
      intrachromosome_links(0) {}
  
  Window(const Window&) = default;
  Window& operator=(const Window&) = default;

  Window(Window&&) = default;
  Window& operator=(Window&&) = default;

  ~Window() = default;
};

struct Node{
  uint32_t id;
  uint32_t len;
  Pileogram pileogram;
  uint32_t interchromosome_links;
  uint32_t intrachromosome_links;
  uint32_t link_0011;
  uint32_t link_0110;
  std::vector<Window> windows;
  Node(uint32_t id, uint32_t len, int window_size)
      : id(id),
        len(len),
        pileogram(id, len),
        interchromosome_links(0), 
        intrachromosome_links(0),
        link_0011(0),
        link_0110(0) {
          int size = ceil(len/window_size);
          for (int i = 0; i <= size; i++) {
            Window w(i);
            windows.push_back(w);
          }
        }

  Node(const Node&) = default;
  Node& operator=(const Node&) = default;

  Node(Node&&) = default;
  Node& operator=(Node&&) = default;

  ~Node() = default;
};

class Graph {
 public:

  explicit Graph(std::shared_ptr<thread_pool::ThreadPool> thread_pool = nullptr, int window_size = 100000);  // NOLINT

  Graph(const Graph&) = delete;
  Graph& operator=(const Graph&) = delete;

  Graph(Graph&&) = default;
  Graph& operator=(Graph&&) = default;

  ~Graph() = default;

  void Construct(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>& targets,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences);
  
  void PrintJson(const std::string& path) const;

 private:
  friend class cereal::access; 
  std::shared_ptr<thread_pool::ThreadPool> thread_pool_;
  std::unordered_map<std::uint32_t, Node> contigs;
  std::vector<std::vector<std::uint32_t>> adjMatrix;
  std::unordered_map<std::string, std::vector<biosoup::Overlap>> read_pairs;
  std::vector<std::unordered_map<std::string, std::vector<biosoup::Overlap>>> all_read_pairs;
  std::unordered_map<std::string, std::vector<biosoup::Overlap>> interchromosome_read_pairs;
  uint32_t window_size;
    
  template<class Archive>
  void serialize(Archive & archive) {
    archive(CEREAL_NVP(read_pairs), CEREAL_NVP(interchromosome_read_pairs)); // serialize things by passing them to the archive
  }

  void Process(
    std::vector<std::future<std::vector<std::pair<std::string, std::vector<biosoup::Overlap>>>>>& futures,
    ram::MinimizerEngine& minimizer_engine,
    std::unique_ptr<biosoup::NucleicAcid>& sequence1,
    std::unique_ptr<biosoup::NucleicAcid>& sequence2);

  std::pair<std::uint32_t, std::uint32_t> GetOverlap(
    biosoup::Overlap ol1,
    biosoup::Overlap ol2); 

  void FillPileogram(std::unordered_map<std::string, std::vector<biosoup::Overlap>>& read_pairs);
  void CreateGraph(std::vector<std::unique_ptr<biosoup::NucleicAcid>>& targets);

  int GetNumWindows();

  void CalcualteInterChromosomeLinks(
  std::unordered_map<std::string, std::vector<biosoup::Overlap>>& interchromsome_read_pairs);

  void GenerateMatrix(
    std::vector<std::vector<std::uint32_t>> &matrix, 
    std::unordered_map<std::string, std::vector<biosoup::Overlap>>& interchromsome_read_pairs);
  
  void GenerateMatrixWindowIntraLinks(
    std::vector<int>& window_id_map,
    std::vector<std::vector<std::uint32_t>> &matrix,
    std::unordered_map<std::string, std::vector<biosoup::Overlap>>& read_pairs);
  
  void GenerateMatrixAftSplitWindow(
  std::vector<int>& window_id_map,
  std::unordered_map<int, int>& window_split_map,
  std::vector<std::vector<std::uint32_t>> &matrix,
  std::unordered_map<std::string, std::vector<biosoup::Overlap>>& read_pairs);

  void GenerateMatrixAftSplitWindow(
    int contig_id,
    std::unordered_map<int, pair<int, bool>>& window_split_map,
    std::vector<std::vector<std::uint32_t>> &matrix,
    std::unordered_map<std::string, std::vector<biosoup::Overlap>>& read_pairs);

  // REMOVE
  void GenerateMatrixWindowIntraLinks(
    std::vector<int>& window_id_map,
    std::vector<std::vector<std::uint32_t>> &matrix,
    std::unordered_map<std::string, std::vector<biosoup::Overlap>>& read_pairs,
    int window_size);

  void GenerateMatrixWindow(
    std::vector<int>& window_id_map,
    std::vector<std::vector<std::uint32_t>> &matrix,
    std::unordered_map<std::string, std::vector<biosoup::Overlap>>& interchromsome_read_pairs);
  
  std::vector<std::vector<uint32_t>> GetComponents(std::vector<std::vector<std::uint32_t>> &matrix);

  std::vector<int> GenerateMapWindowID();
  std::uint32_t MaxInMatrix(std::vector<std::vector<std::uint32_t>> &matrix);
  std::uint32_t MaxInter(std::vector<std::vector<std::uint32_t>> &matrix);
  std::uint32_t MaxIntra(std::vector<std::vector<std::uint32_t>> &matrix);
  int FindContigID(int window_id, std::vector<int>& window_id_map, int *begin, int *end);
  
  void Store(int i) const;
  void Load(int i);
};

}  // namespace tarantula

#endif  // TARANTULA_GRAPH_HPP_
