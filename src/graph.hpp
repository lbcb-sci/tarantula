// Copyright (c) 2020 Robert Vaser

#ifndef TARANTULA_GRAPH_HPP_
#define TARANTULA_GRAPH_HPP_

#include <memory>

#include "biosoup/sequence.hpp"
#include "thread_pool/thread_pool.hpp"

namespace tarantula {

class Graph {
 public:
  explicit Graph(
      std::shared_ptr<thread_pool::ThreadPool> thread_pool = nullptr);

  Graph(const Graph&) = delete;
  Graph& operator=(const Graph&) = delete;

  Graph(Graph&&) = default;
  Graph& operator=(Graph&&) = default;

  ~Graph() = default;

 private:
  std::shared_ptr<thread_pool::ThreadPool> thread_pool_;
};

}  // namespace tarantula

#endif  // TARANTULA_GRAPH_HPP_
