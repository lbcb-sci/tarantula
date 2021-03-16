#include "pileogram.hpp"

namespace tarantula{

Pileogram::Pileogram(std::uint32_t id, std::uint32_t contigLen){
  this->id = id; 
  this->contigLen = contigLen; 
  this->begin = 0; 
  this->end = contigLen; 
  initData(); 
}

Pileogram::Pileogram(){
  
}

void Pileogram::addLayer(std::vector<biosoup::Overlap> overlap){
  for (const auto& ol: overlap){
    for (u_int32_t i=ol.rhs_begin;i<=ol.rhs_end;i++){
      data[i] += 1; 
    }
  }
}

void Pileogram::addLayer(uint32_t begin, uint32_t end){
  for (uint32_t i=begin;i<=end;i++){
    data[i] += 1;
  }
}

void Pileogram::initData(){
  for (uint32_t i=0;i<contigLen;i++){
    data.push_back(0); 
  }
}
}