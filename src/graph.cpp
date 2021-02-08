// Copyright (c) 2020 Robert Vaser

#include "graph.hpp"

namespace tarantula {

Graph::Graph(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool)
    : thread_pool_(thread_pool ?
        thread_pool :
        std::make_shared<thread_pool::ThreadPool>(1)) {
}

}  // namespace tarantula
