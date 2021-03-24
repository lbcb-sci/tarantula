// Copyright (c) 2020 Robert Vaser

#ifndef TARANTULA_GRAPH_HPP_
#define TARANTULA_GRAPH_HPP_

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

#include "pileogram.hpp"

namespace tarantula {

struct Node{
  uint32_t id;
  uint32_t len;
  Pileogram pileogram;
  uint32_t interchromosome_links;
  uint32_t intrachromosome_links;
  Node(uint32_t id, uint32_t len)
      : id(id),
        len(len),
        pileogram(id, len),
        interchromosome_links(0), 
        intrachromosome_links(0) {}

  Node(const Node&) = default;
  Node& operator=(const Node&) = default;

  Node(Node&&) = default;
  Node& operator=(Node&&) = default;

  ~Node() = default;
};

class Graph {
 public:
  explicit Graph(std::shared_ptr<thread_pool::ThreadPool> thread_pool = nullptr);  // NOLINT

  Graph(const Graph&) = delete;
  Graph& operator=(const Graph&) = delete;

  Graph(Graph&&) = default;
  Graph& operator=(Graph&&) = default;

  ~Graph() = default;

  void Construct(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>& targets,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences);

  void Process(
    std::vector<std::future<std::pair<std::string, std::vector<std::vector<biosoup::Overlap>>>>>& futures,
    ram::MinimizerEngine& minimizer_engine,
    std::unique_ptr<biosoup::NucleicAcid>& sequence1,
    std::unique_ptr<biosoup::NucleicAcid>& sequence2);

  std::pair<std::uint32_t, std::uint32_t> GetOverlap(
    biosoup::Overlap ol1,
    biosoup::Overlap ol2); 

  void FillPileogram(std::unordered_map<std::string, std::vector<std::vector<biosoup::Overlap>>>& read_pairs);
  void CreateGraph(std::vector<std::unique_ptr<biosoup::NucleicAcid>>& targets);

  void PrintJson(const std::string& path) const;

  void CalcualteInterChromosomeLinks(
  std::unordered_map<std::string, std::vector<std::vector<biosoup::Overlap>>>& interchromsome_read_pairs);

 private:
  std::shared_ptr<thread_pool::ThreadPool> thread_pool_;
  std::unordered_map<std::uint32_t, Node> contigs;
};

}  // namespace tarantula

#endif  // TARANTULA_GRAPH_HPP_
