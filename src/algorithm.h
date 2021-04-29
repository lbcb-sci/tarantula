#include <vector>

#include "tree.h"

namespace directedforce {

void initVerticesPosition(std::vector<std::shared_ptr<Vertex>>& vertices); 

void directedForceAlgorithm(std::vector<std::shared_ptr<Vertex>>& vertices, std::vector<std::vector<double>>& adjMax);

void GenerateGraphFromDirectedForceAlgorithm(std::string input, std::string output);

}