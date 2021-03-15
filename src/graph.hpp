// Copyright (c) 2020 Robert Vaser

#ifndef TARANTULA_GRAPH_HPP_
#define TARANTULA_GRAPH_HPP_

#include <memory>
#include <unordered_map>
#include <tuple>

#include "biosoup/overlap.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "thread_pool/thread_pool.hpp"
#include "ram/minimizer_engine.hpp"

namespace tarantula {

class Graph {
 public:
  explicit Graph(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool = nullptr);
    
  Graph(const Graph&) = delete;
  Graph& operator=(const Graph&) = delete;

  Graph(Graph&&) = default;
  Graph& operator=(Graph&&) = default;

  void Construct(std::vector<std::unique_ptr<biosoup::NucleicAcid>>& targets ,std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences); 
  void Process(std::vector<std::future<void>>& futures, ram::MinimizerEngine& minimizer_engine,std::unique_ptr<biosoup::NucleicAcid>& sequence1, std::unique_ptr<biosoup::NucleicAcid>& sequence2, std::unordered_map<std::string, std::vector<std::vector<biosoup::Overlap>>>& readPairs); 
  std::tuple<std::uint32_t, std::uint32_t> GetOverlap(biosoup::Overlap ol1, biosoup::Overlap ol2);

  ~Graph() = default;

 private:
  std::shared_ptr<thread_pool::ThreadPool> thread_pool_;
};

}  // namespace tarantula

#endif  // TARANTULA_GRAPH_HPP_
