#include <chrono>
#include <random>
#include <memory>
#include <math.h>  
#include <cmath>  
#include <unordered_map>
#include <fstream>
#include <iostream>

#include "algorithm.h"
#include "progressBar.h"

#include <iostream>

namespace directedforce {


double af(double k, double x, double weight){
  return x*x/k*weight;
}

//k is const
double rf(double k, double z){
  double k2 = k; 
  return -k2*k2/z; 
}

double cool(double t){
  if (t > 0.001)
    return t*0.99; 
  else 
    return 0.001; 
}

double fRand(double fMin, double fMax){
  double f = (double)rand() / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

void initVerticesPosition(std::vector<std::shared_ptr<Vertex>>& vertices){
  int numV = vertices.size(); 
  int iter=1; 
  double angle;

  angle = 2.0 * M_PI / numV;
  for (int i=0;i<numV;i++){
    vertices[i]->pos.x = cos(angle*i); 
    vertices[i]->pos.y = sin(angle*i); 
  }
}

void calculateAttrativeForce(std::vector<std::shared_ptr<Vertex>>& vertices, std::vector<std::vector<double>>& adjMax, double k){
  MathVector diff; 
  double diffABS; 
  int numVertices = vertices.size(); 
  for (int i=0; i<numVertices; i++){
    for (int r=0; r<i; r++){
      if (adjMax[i][r]>0){
        diff = vertices[i]->pos - vertices[r]->pos; 
        diffABS = diff.abs(); 
        if (diffABS!=0){
          vertices[i]->disp -= diff/diffABS*af(k,diffABS,adjMax[i][r]); 
          vertices[r]->disp += diff/diffABS*af(k,diffABS,adjMax[i][r]); 
        }
      }
    }
  }
}

bool insert(std::shared_ptr<Node>& node, std::shared_ptr<Vertex>& particle){
  if (node->box.in(particle->pos)){
    if (node->noParticles()){
      node->n = particle;  
    }
    else{
      //node have children
      //get quadrant
      if (node->n!=nullptr){ 
        if (node->n->pos.x == particle->pos.x && node->n->pos.y == particle->pos.y){
          return false; 
          particle->pos += (node->box.c2 - node->box.c1)*2; 
        }
        std::shared_ptr<Node>& nQuadrant = node->getQuadrant(node->n->pos); 
        if (nQuadrant->noParticles()){
          nQuadrant->n = node->n; 
        }
        else{
          insert(nQuadrant, node->n); 
        }  
        node->n = nullptr;  
      }
      std::shared_ptr<Node>& quadrant = node->getQuadrant(particle->pos);
      if (quadrant->noParticles()){
        quadrant->n = particle; 
      }
      else{
        insert(quadrant, particle); 
      }
    }
  }
  else{
    std::cerr << "[GraphVisualisation::WARNING] Increase width/length" << std::endl; 
  }
  return true; 
}

Box getBoundingBox(std::vector<std::shared_ptr<Vertex>>& vertices){
  double xMin=0, yMin=0, xMax=0, yMax=0; 
  for (int i=0;i<vertices.size();i++){
    if (vertices[i]->pos.x < xMin){
      xMin = vertices[i]->pos.x; 
    }
    if (vertices[i]->pos.x > xMax){
      xMax = vertices[i]->pos.x; 
    }

    if (vertices[i]->pos.y < yMin){
      yMin = vertices[i]->pos.y; 
    }
    
    if (vertices[i]->pos.y > yMax){
      yMax = vertices[i]->pos.y; 
    }
  }
  Box box = {{xMin,yMin}, {xMax,yMin},{xMax,yMax}, {xMin,yMax}}; 
  return box; 
}

void initTree(std::shared_ptr<Node>& root, std::vector<std::shared_ptr<Vertex>>& vertices) {
  root->first = nullptr; 
  root->second = nullptr;
  root->third = nullptr; 
  root->fourth = nullptr; 
  root->n = nullptr;
  
  root->box = getBoundingBox(vertices);     
}

void generateTree(std::vector<std::shared_ptr<Vertex>>& vertices, std::shared_ptr<Node>& root){
 
  int numVertices = vertices.size(); 
  //init tree
  initTree(root, vertices); 
  for (int i=0; i<numVertices; i++){
    auto result = insert(root, vertices[i]); 
    while (!result){
      vertices[i]->pos.x = fRand(0,1000); 
      vertices[i]->pos.y = fRand(0,1000); 
      result = insert(root, vertices[i]); 
    }
  } 
}

void computeMassDistribution(std::shared_ptr<Node>& node, double mass=1){
  if (node->numChild()==0 && node->n!=nullptr){
    node->mass = mass;
    node->centreOfMass.x = node->n->pos.x; 
    node->centreOfMass.y = node->n->pos.y; 
  }
  else{
    if (node->first!=nullptr&&!node->first->noParticles()){
      computeMassDistribution(node->first, mass); 
      node->mass += node->first->mass; 
      node->centreOfMass += node->first->centreOfMass; 
    }
    if (node->second!=nullptr&&!node->second->noParticles()){
      computeMassDistribution(node->second, mass); 
      node->mass += node->second->mass; 
      node->centreOfMass += node->second->centreOfMass; 
    }
    if (node->third!=nullptr&&!node->third->noParticles()){
      computeMassDistribution(node->third, mass); 
      node->mass += node->third->mass; 
      node->centreOfMass += node->third->centreOfMass; 
    }
    if (node->fourth!=nullptr&&!node->fourth->noParticles()){
      computeMassDistribution(node->fourth, mass); 
      node->mass += node->fourth->mass; 
      node->centreOfMass += node->fourth->centreOfMass;  
    }
    node->centreOfMass /=  node->mass; 
  }
} 

MathVector calculateForceBarnesHutPerVertex(std::shared_ptr<Node>& node, std::shared_ptr<Vertex>& targetParticle, double k, double theta_const){
  double distance, height, theta;
  MathVector force = {0,0}; 

  //if it is a leaf
  if (node->numChild()==0&&node->n!=nullptr){
    MathVector diff = node->centreOfMass-targetParticle->pos;
    distance = diff.abs(); 
    
    if (distance==0){
      return {0,0};
    }
    height = node->box.c2.x - node->box.c1.x ; 
    force = (diff/distance)*(node->mass*rf(k,distance));
  }
  else{
    MathVector diff = node->centreOfMass-targetParticle->pos; 
    distance = diff.abs(); 
    height = node->box.c2.x - node->box.c1.x ; 
    theta = height/distance; 

    if (distance==0){
      return {0,0};
    }
       
    if (distance!=distance){
      std::cerr << "[GraphVisualisation::Warning] W & L must be a bigger number" << std::endl; 
      return {0,0}; 
    }
    if (theta < theta_const){
      auto temp = diff/distance; 
      force = (diff/distance)*(node->mass*rf(k,distance));
    }
    else{
      if (node->first!=nullptr){
        force += calculateForceBarnesHutPerVertex(node->first, targetParticle, k, theta_const); 
      }
      if (node->second!=nullptr){
        force += calculateForceBarnesHutPerVertex(node->second, targetParticle, k, theta_const); 
      }
      if (node->third!=nullptr){
        force += calculateForceBarnesHutPerVertex(node->third, targetParticle, k, theta_const); 
      }
      if (node->fourth!=nullptr){
        force += calculateForceBarnesHutPerVertex(node->fourth, targetParticle, k, theta_const); 
      }
    }
  }
  return force; 
}

void calculateRepulsiveForce_barnesHutAlgo(std::vector<std::shared_ptr<Vertex>>& vertices, std::vector<std::vector<double>>& adjMax, double k, double theta){
  MathVector diff, force; 
  double diffABS, abs; 
  int numVertices = vertices.size(); 
  //generate tree
  std::shared_ptr<Node> tree = std::make_shared<Node>();
  generateTree(vertices, tree); 
  computeMassDistribution(tree);
  for (int i=0; i<numVertices; i++){
    force = calculateForceBarnesHutPerVertex(tree, vertices[i], k, theta); 
    vertices[i]->disp = force; 
  }
}

void calculateForceBarnesHut(std::vector<std::shared_ptr<Vertex>>& vertices, std::vector<std::vector<double>>& adjMax, double k, double theta){ 
  calculateRepulsiveForce_barnesHutAlgo(vertices,adjMax, k, theta); 
  calculateAttrativeForce(vertices,adjMax,k); 
}

void directedForceAlgorithm(std::vector<std::shared_ptr<Vertex>>& vertices, std::vector<std::vector<double>>& adjMax){
  int numVertices = vertices.size(); 
  double k = sqrt(1000*1000/numVertices); 
  int iterations = 1000; 
  double theta = 0.01;
  double coeff, t; 

  MathVector diff; 
  double diffABS, abs; 

  //ProgressBar bar;
  //bar.set_bar_width(50);
  //bar.fill_bar_progress_with("■");
  //bar.fill_bar_remainder_with(" ");
  //float progress;
  //in each iterations
  t = 1; 
  for (int iter=1; iter<=iterations; iter++){
    calculateForceBarnesHut(vertices, adjMax, k, theta); 

    //for both different algorithm, you will need this
    for (int i=0; i<vertices.size(); i++){      
      abs = vertices[i]->disp.abs(); 
      vertices[i]->pos += vertices[i]->disp/abs * std::min(abs, t); 
    }

    t = cool(t); 
    // progress bar
    //progress = ((double) iter/iterations )* 100;
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
  double indexN1, indexN2; 
  int iN1, iN2, end, wl, weight; 
  bool prev = false; 
  std::vector<bool> temp; 
  std::vector<Vertex> vtemp; 
  bool has_weight = false;
  outfile.open(outputPath);
  bool prev_color = false;
  int new_node = 0;
  while (getline (infile, text)) {
    prev = false; 
    //std::cout << text << endl;
    
    if (with_coloring && !prev_color) {
      if (text[0] != '%') {
        outfile << text + "\n";
        continue;
      } else if (text[0] == '%') {
        outfile << text + "\n";
        prev_color = true;
        continue;
      }
    }
    for (int i = 0; i < text.size(); i++) {
      if (text[i] == '-' && !prev) {
        prev = true; 
        iN1 = i;
      } else if (text[i] != '-' && prev) {
        iN2 = i; 
        prev = false; 
      } else if (text[i] == ',') {
        has_weight = true;
        wl = i;
      } else if (text[i] == '\n') {
        end = i; 
      }
    }
    n1 = text.substr(0, iN1);
    if (!has_weight) {
      int len = end - iN2; 
      n2 = text.substr(iN2, len);
    } else {
      int len = wl - iN2;
      n2 = text.substr(iN2, len);
      weight = std::stoi(text.substr(wl+1, end-wl));
    }

    outfile << text + "\n";
    // Output the text from the file
    auto it1 = table.find(n1); 
    if (it1 == table.end()) {
      //add n1 to table & assign index
      indexN1 = table.size(); 
      table.insert(make_pair(n1, table.size())); 
      vtemp.push_back({{0, 0}, {0, 0}}); 
      new_node++;
    } else {
      //get index
      indexN1 = it1->second; 
    }

    auto it2 = table.find(n2); 
    if (it2 == table.end()) {
      //add n1 to table & assign index
      indexN2 = table.size(); 
      table.insert(make_pair(n2, table.size())); 
      vtemp.push_back({{0, 0}, {0, 0}});
      new_node++;
    } else {
      //get index
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
    //edges.push_back({}); 
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

  std::cerr << "[GraphVisualisation] number of nodes.." << vertices.size() << std::endl;
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

void GenerateGraphFromDirectedForceAlgorithm(std::string input, std::string output) {
  std::vector<std::shared_ptr<Vertex>> vertices;
  std::vector<std::vector<double>> edges;
  std::cerr << "[GraphVisualisation] Reading vertices" << std::endl;
  std::unordered_map<std::string, int> map_table;
  
  map_table = parseTxtFile(input, vertices, edges, output, false);
  initVerticesPosition(vertices); 

  //ModerateEdges(edges, vertices.size());

  /*
  for (auto r:edges){
    for (auto c:r){
      cerr << c << " ";
    }
    cerr<<endl;
  }*/
  
  std::cerr << "[GraphVisualisation] calculating, iterations: 1000" << std::endl;
  // for debug only
  // generateOutputFile(argv[optind], initpath, vertices, width, length);
  directedForceAlgorithm(vertices, edges);
  std::cerr << "[GraphVisualisation] Generating output" << std::endl;
  generateOutputFile(input, output, vertices, map_table);
}

}