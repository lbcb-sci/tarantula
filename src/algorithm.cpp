// Copyright (c) 2021 Cecilia Lee, Robert Vaser

#include <math.h>  

#include <random>
#include <memory>
#include <cmath>  
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <utility>
#include <string>

#include "algorithm.h"
#include "progressBar.h"

using namespace std; 

namespace directedforce {
// k is const
double af(double k, double x) {
  return x*x/k; 
}

double af(double k, double x, double weight) {
  return x*x/k*weight;
}

// k is const
double rf(double k, double z) {
  return -k*k/z; 
}

double cool(double t) {
  if (t > 0.001)
    return t*0.99; 
  else 
    return 0.001; 
}

double fRand(double fMin, double fMax) {
  double f = static_cast<double>(rand()) / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

void initVerticesPosition(vector<shared_ptr<Vertex>>& vertices, double xMax, double yMax, bool random) {
  if (random) {
    for (int i = 0; i < vertices.size(); i++) {
      vertices[i]->pos.x = fRand(0, xMax); 
      vertices[i]->pos.y = fRand(0, yMax); 
    }
  } else {
    int numV = vertices.size(); 
    double angle;

    angle = 2.0 * M_PI / numV;
    for (int i = 0; i < numV; i++) {
      vertices[i]->pos.x = cos(angle*i); 
      vertices[i]->pos.y = sin(angle*i); 
    }
  }
}

void calculateAttrativeForce(vector<shared_ptr<Vertex>>& vertices, vector<vector<double>>& adjMax, double k) {
  MathVector diff; 
  double diffABS; 
  int numVertices = vertices.size(); 
  for (int i = 0; i < numVertices; i++) {
    for (int r = 0; r < i; r++) {
      if (adjMax[i][r] > 0) {
        diff = vertices[i]->pos - vertices[r]->pos; 
        diffABS = diff.abs(); 
        if (diffABS != 0) {
          vertices[i]->disp -= diff/diffABS*af(k, diffABS, adjMax[i][r]); 
          vertices[r]->disp += diff/diffABS*af(k, diffABS, adjMax[i][r]); 
        }
      }
    }
  }
}

void calculateForceBruteForce(vector<shared_ptr<Vertex>>& vertices, vector<vector<double>>& adjMax, double k) {
  MathVector diff; 
  double diffABS; 
  int numVertices = vertices.size(); 
  for (int i = 0; i < numVertices; i++) {
    vertices[i]->disp = 0; 
    for (int r = 0; r < numVertices; r++) {
      if (i == r)
        continue; 
      diff = vertices[i]->pos - vertices[r]->pos; 
      diffABS = diff.abs(); 
      if (diffABS != 0) {
        vertices[i]->disp -= (diff/diffABS)*rf(k, diffABS); 
      } else {
        vertices[i]->disp -= rf(k, 10); 
      }
    }
  }

  for (int i = 0; i< numVertices; i++) {
    for (int r = 0; r < i; r++) {
      if (adjMax[i][r] > 0) {
        diff = vertices[i]->pos - vertices[r]->pos; 
        diffABS = diff.abs(); 
        if (diffABS != 0) {
          vertices[i]->disp -= diff/diffABS*af(k, diffABS, adjMax[i][r]); 
          vertices[r]->disp += diff/diffABS*af(k, diffABS, adjMax[i][r]); 
        }
      }
    }
  }
}

bool insert(shared_ptr<Node>& node, shared_ptr<Vertex>& particle) {
  if (node->box.in(particle->pos)) {
    if (node->noParticles()) {
      node->n = particle;  
    } else {
      if (node->n != nullptr) { 
        if (node->n->pos.x == particle->pos.x && node->n->pos.y == particle->pos.y) {
          particle->pos += (node->box.c2 - node->box.c1)*2; 
          return false;
        }
        shared_ptr<Node>& nQuadrant = node->getQuadrant(node->n->pos); 
        if (nQuadrant->noParticles()) {
          nQuadrant->n = node->n; 
        } else {
          insert(nQuadrant, node->n); 
        }  
        node->n = nullptr;  
      }
      shared_ptr<Node>& quadrant = node->getQuadrant(particle->pos);
      if (quadrant->noParticles()) {
        quadrant->n = particle; 
      } else {
        insert(quadrant, particle); 
      }
    }
  } else {
    cerr << "[GraphVisualisation::WARNING] Increase width/length" << endl; 
  }
  return true; 
}

Box getBoundingBox(vector<shared_ptr<Vertex>>& vertices) {
  double xMin = 0, yMin = 0, xMax = 0, yMax = 0; 
  for (int i = 0; i < vertices.size(); i++) {
    if (vertices[i]->pos.x < xMin) {
      xMin = vertices[i]->pos.x; 
    }
    if (vertices[i]->pos.x > xMax) {
      xMax = vertices[i]->pos.x; 
    }

    if (vertices[i]->pos.y < yMin) {
      yMin = vertices[i]->pos.y; 
    }
    
    if (vertices[i]->pos.y > yMax) {
      yMax = vertices[i]->pos.y; 
    }
  }
  Box box(xMin, yMin, xMax, yMin, xMax, yMax, xMin, yMax);
  return box; 
}

void initTree(shared_ptr<Node>& root, double width, double length, vector<shared_ptr<Vertex>>& vertices, bool dynamic) {
  root->first = nullptr; 
  root->second = nullptr;
  root->third = nullptr; 
  root->fourth = nullptr; 
  root->n = nullptr;
  
  if (dynamic) {
    root->box = getBoundingBox(vertices);     
  } else {
    Box box(0, 0, width, 0, width, length, 0, length); 
    root->box = box;
  }
}

void generateTree(vector<shared_ptr<Vertex>>& vertices, double width, double length, shared_ptr<Node>& root, bool dynamic) {
 
  int numVertices = vertices.size(); 
  // init tree
  initTree(root, width, length, vertices, dynamic); 
  for (int i = 0; i < numVertices; i++) {
    auto result = insert(root, vertices[i]); 
    while (!result) {
      vertices[i]->pos.x = fRand(0, width); 
      vertices[i]->pos.y = fRand(0, length); 
      result = insert(root, vertices[i]); 
    }
  } 
}

void computeMassDistribution(shared_ptr<Node>& node, double mass = 1) {
  if (node->numChild() == 0 && node->n != nullptr) {
    node->mass = mass;
    node->centreOfMass.x = node->n->pos.x; 
    node->centreOfMass.y = node->n->pos.y; 
  } else {
    if (node->first != nullptr && !node->first->noParticles()) {
      computeMassDistribution(node->first, mass); 
      node->mass += node->first->mass; 
      node->centreOfMass += node->first->centreOfMass; 
    }
    if (node->second != nullptr && !node->second->noParticles()) {
      computeMassDistribution(node->second, mass); 
      node->mass += node->second->mass; 
      node->centreOfMass += node->second->centreOfMass; 
    }
    if (node->third != nullptr && !node->third->noParticles()) {
      computeMassDistribution(node->third, mass); 
      node->mass += node->third->mass; 
      node->centreOfMass += node->third->centreOfMass; 
    }
    if (node->fourth != nullptr && !node->fourth->noParticles()) {
      computeMassDistribution(node->fourth, mass); 
      node->mass += node->fourth->mass; 
      node->centreOfMass += node->fourth->centreOfMass;  
    }
    node->centreOfMass /=  node->mass; 
  }
} 

MathVector calculateForceBarnesHutPerVertex(shared_ptr<Node>& node, shared_ptr<Vertex>& targetParticle, double k, double theta_const){
  double distance, height, theta;
  MathVector force = {0, 0}; 

  // if it is a leaf
  if (node->numChild() == 0 && node->n != nullptr) {
    MathVector diff = node->centreOfMass-targetParticle->pos;
    distance = diff.abs(); 
    
    if (distance == 0) {
      return {0, 0};
    }
    height = node->box.c2.x - node->box.c1.x ; 
    force = (diff/distance)*(node->mass*rf(k, distance));
  } else {
    MathVector diff = node->centreOfMass-targetParticle->pos; 
    distance = diff.abs(); 
    height = node->box.c2.x - node->box.c1.x ; 
    theta = height/distance; 

    if (distance == 0) {
      return {0, 0};
    }
       
    if (distance != distance) {
      std::cerr << "[GraphVisualisation::Warning] W & L must be a bigger number" << endl; 
      return {0, 0}; 
    }
    if (theta < theta_const) {
      auto temp = diff/distance; 
      force = (diff/distance)*(node->mass*rf(k, distance));
    } else {
      if (node->first != nullptr) {
        force += calculateForceBarnesHutPerVertex(node->first, targetParticle, k, theta_const); 
      }
      if (node->second != nullptr) {
        force += calculateForceBarnesHutPerVertex(node->second, targetParticle, k, theta_const); 
      }
      if (node->third != nullptr) {
        force += calculateForceBarnesHutPerVertex(node->third, targetParticle, k, theta_const); 
      }
      if (node->fourth != nullptr) {
        force += calculateForceBarnesHutPerVertex(node->fourth, targetParticle, k, theta_const); 
      }
    }
  }
  return force; 
}

void calculateRepulsiveForce_barnesHutAlgo(
  vector<shared_ptr<Vertex>>& vertices,
  vector<vector<double>>& adjMax,
  double k,
  double width,
  double length,
  double mass,
  bool dynamic,
  double theta) {
  MathVector diff, force; 
  double diffABS, abs; 
  int numVertices = vertices.size(); 
  // generate tree
  shared_ptr<Node> tree = make_shared<Node>();
  generateTree(vertices, width, length, tree, dynamic); 
  computeMassDistribution(tree, mass);
  for (int i = 0; i < numVertices; i++) {
    force = calculateForceBarnesHutPerVertex(tree, vertices[i], k, theta); 
    vertices[i]->disp = force;
  }
}

void calculateForceBarnesHut(
  vector<shared_ptr<Vertex>>& vertices,
  vector<vector<double>>& adjMax,
  double k,
  double width,
  double length,
  double mass,
  bool dynamic,
  double theta) { 

  calculateRepulsiveForce_barnesHutAlgo(
    vertices,
    adjMax,
    k,
    width,
    length, 
    mass,
    dynamic, 
    theta);
  calculateAttrativeForce(vertices, adjMax, k); 
}

void directedForceAlgorithm(
  vector<shared_ptr<Vertex>>& vertices,
  vector<vector<double>>& adjMax,
  int L,
  int W,
  int iterations,
  int algoType,
  double theta,
  double mass,
  bool dynamic) {

  int numVertices = vertices.size(); 
  int area = W*L;  
  double k = sqrt(area/numVertices); 
  double coeff, t; 

  MathVector diff; 
  double diffABS, abs; 

  //ProgressBar bar;
  //bar.set_bar_width(50);
  //bar.fill_bar_progress_with("■");
  //bar.fill_bar_remainder_with(" ");
  //float progress;
  // in each iterations
  t = 1; 
  for (int iter = 1; iter <= iterations; iter++) {
    if (algoType == 1) {
      calculateForceBruteForce(vertices, adjMax, k); 
    } else {
      // by default
      calculateForceBarnesHut(vertices, adjMax, k, W, L, mass, dynamic, theta); 
    }

    // for both different algorithm, you will need this
    for (int i = 0; i < vertices.size(); i++) {      
      abs = vertices[i]->disp.abs(); 
      vertices[i]->pos += vertices[i]->disp/abs * min(abs, t); 
    }

    t = cool(t); 
    // progress bar
    //progress = (static_cast<double>(iter/iterations))* 100;
    //bar.update(progress);
  }
  std::cerr << std::endl;
}

void duplicateTxtFile(std::string file, std::string copyPath) {
  std::ifstream inFile(file,       std::ios::binary);
  std::ofstream outFile(copyPath, std::ios::binary);
  outFile << inFile.rdbuf();
  outFile.close();
}

std::unordered_map<std::string, int> parseTxtFile(
  std::string path,
  std::vector<std::shared_ptr<Vertex>>& vertices,
  std::vector<std::vector<double>>& edges, 
  std::string outputPath,
  bool with_coloring) {

  std::string text, n1, n2;
  std::unordered_map<std::string, int> table; 
  std::ifstream infile(path);
  std::ofstream outfile;
  double indexN1 = 0, indexN2 = 0; 
  int iN1 = 0, iN2 = 0, end = 0, wl = 0;
  double weight = 0.0;
  bool prev = false; 
  std::vector<bool> temp; 
  std::vector<Vertex> vtemp; 
  bool has_weight = false;
  bool prev_color = false;
  int new_node = 0;
  std::string weight_str;
  outfile.open(outputPath);
  while (getline(infile, text)) {
    prev = false;
    
    
    for (int i = 0; i < text.size(); i++) {
      if (text[i] == '-' && text[i+1] == '-') {
        iN1 = i;
        iN2 = i+2;
        prev = false; 
      } else if (text[i] == ',') {
        has_weight = true;
        wl = i;
        break;
      }
    }
    end = text.size();
    //std::cerr << "in1: " << iN1 << " ,in2: " << iN2 << std::endl;
    n1 = text.substr(0, iN1);
    if (!has_weight) {
      int len = end - iN2; 
      n2 = text.substr(iN2, len);
    } else {
      int len = wl - iN2;
      n2 = text.substr(iN2, len);
      weight_str = text.substr(wl+1, end);
      //std:: cerr << "weight_str: " << weight_str <<"|" << std::endl;
      //std::cerr << "n1: " << n1 << " ,n2: " << n2 << " wl: " << wl << " end: " << end << std::endl;
      if (weight_str.find('e') < weight_str.length()) {
        // exponential value
        std::istringstream os(weight_str);
        os >> weight;
      } else  {
        weight = std::stoi(weight_str);
      }
      //std::cerr << "weight" << weight << std::endl;
      /*
      if (weight<50)
        continue;*/
    }

    outfile << text + "\n";
    // Output the text from the file
    auto it1 = table.find(n1); 
    if (it1 == table.end()) {
      // add n1 to table & assign index
      indexN1 = table.size(); 
      table.insert(make_pair(n1, table.size())); 
      vtemp.push_back({{0, 0}, {0, 0}}); 
      new_node++;
    } else {
      // get index
      indexN1 = it1->second; 
    }

    auto it2 = table.find(n2); 
    if (it2 == table.end()) {
      // add n1 to table & assign index
      indexN2 = table.size(); 
      table.insert(make_pair(n2, table.size())); 
      vtemp.push_back({{0, 0}, {0, 0}});
      new_node++;
    } else {
      // get index
      indexN2 = it2->second; 
    }
    if (new_node > 0) {
      int num_edges = edges.size();
      for (auto& col : edges) {
        for (int i = 0; i < new_node; i++) {
          col.push_back(0);
        }
      }
      for (int i = 0; i < new_node; i++) {
        std::vector<double> temp(table.size(), 0);
        edges.push_back(temp);
      }

    }
    
    if (has_weight) {
      edges[indexN1][indexN2] = weight; 
      edges[indexN2][indexN1] = weight;  
    } else {
      edges[indexN1][indexN2] = 1;
      edges[indexN2][indexN1] = 1;  
    }
    new_node = 0;
  }
  outfile.close();

  for (int i = 0; i < vtemp.size(); i++) {
    vertices.push_back(std::make_shared<Vertex>(vtemp[i])); 
  }

  //std::cerr << "[GraphVisualisation] number of nodes.." << vertices.size() << std::endl;
  return table;
}

void generateOutputFile(
  std::string inputPath,
  std::string path,
  std::vector<std::shared_ptr<Vertex>>& vertices,
  std::unordered_map<std::string, int> mapTable) {

  std::ifstream infile(inputPath); 
  std::ofstream outfile(path, std::ios_base::app); 

  outfile << "-" << std::endl; 
  for (int i = 0; i < vertices.size(); i++) {
      std::string vertex;
      for (auto map : mapTable) {
        if (map.second == i) {
          vertex = map.first;
          break;
        }
      }
      outfile << vertex << "|(" << vertices[i]->pos.x 
              << "," << vertices[i]->pos.y << ")" << std::endl;
  }
  outfile << "^" << std::endl; 
  outfile << "50,50"<< std::endl; 
}

void GenerateGraphFromDirectedForceAlgorithm(
  std::string input,
  std::string output) {
  // params
  int iterations = 1000, width = 1000, length = 1000;  
  int algoType = 1;  // default is barnes hut
  bool dynamic = true, random = false, color = false;
  double mass = 1, theta = 0.01; 
  std::vector<std::vector<double>> edges;
  std::unordered_map<std::string, int> map_table;
  std::vector<std::shared_ptr<Vertex>> vertices;
  //std::cerr << "\n[GraphVisualisation] " << input << std::endl;
  //std::cerr << "[GraphVisualisation] Reading vertices" << std::endl;

  if (vertices.size() < 1000) {
    width = vertices.size();
    length = vertices.size();
  }
  
  map_table = parseTxtFile(input, vertices, edges, output, color);
  if (vertices.size() == 0) {
    return;
  }
  
  initVerticesPosition(vertices, width, length, random); 

  //std::cerr << "[GraphVisualisation] calculating, iterations: " 
  //          << iterations << std::endl;
  directedForceAlgorithm(
    vertices, 
    edges, 
    width, 
    length, 
    iterations, 
    algoType, 
    theta, 
    mass, 
    dynamic);

  //std::cerr << "[GraphVisualisation] Generating output" << std::endl;
  generateOutputFile(input, output, vertices, map_table);
    
  //std::cerr << "[GraphVisualisation] txt file generated" << std::endl;
}

void GenerateGraphFromDirectedForceAlgorithm(
  std::string input,
  std::string output,
  std::vector<std::shared_ptr<Vertex>>& vertices,
  std::vector<std::vector<double>>& edges,
  std::unordered_map<std::string, int>& map_table) {
  // params
  int iterations = 1000, width = 1000, length = 1000;  
  int algoType = 1;  // default is barnes hut
  bool dynamic = true, random = false, color = false;
  double mass = 1, theta = 0.01; 

  //std::cerr << "\n[GraphVisualisation] " << input << std::endl;
  //std::cerr << "[GraphVisualisation] Reading vertices" << std::endl;

  map_table = parseTxtFile(input, vertices, edges, output, color);
  if (vertices.size() == 0) {
    return;
  }
  
  if (vertices.size() < 1000) {
    width = vertices.size();
    length = vertices.size();
  }
  

  initVerticesPosition(vertices, width, length, random); 

  //std::cerr << "[GraphVisualisation] calculating, iterations: " 
  //          << iterations << std::endl;
  directedForceAlgorithm(
    vertices, 
    edges, 
    width, 
    length, 
    iterations, 
    algoType, 
    theta, 
    mass, 
    dynamic);

  //std::cerr << "[GraphVisualisation] Generating output" << std::endl;
  generateOutputFile(input, output, vertices, map_table);
    
  //std::cerr << "[GraphVisualisation] txt file generated" << std::endl;
}

}  // namespace directedforce
