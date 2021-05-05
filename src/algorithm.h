// Copyright (c) 2021 Cecilia Lee, Robert Vaser

#ifndef TARANTULA_ALGORITHN_HPP_
#define TARANTULA_ALGORITHM_HPP_


#include <vector>

#include "tree.h"

namespace directedforce {

void GenerateGraphFromDirectedForceAlgorithm(
  std::string input,
  std::string output);

void GenerateGraphFromDirectedForceAlgorithm(
  std::string input,
  std::string output,
  std::vector<std::shared_ptr<Vertex>>& vertices,
  std::vector<std::vector<double>>& edges,
  std::unordered_map<std::string, int>& map_table);

}  // namespace directedforce

#endif  // TARANTULA_ALGORITHM_HPP_