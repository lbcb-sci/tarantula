// Copyright (c) 2020 Robert Vaser

#include <unordered_map>
#include <tuple>

#include "biosoup/overlap.hpp"
#include "ram/minimizer_engine.hpp"

#include "graph.hpp"

//to remove 
#include <iostream>

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace tarantula {

Graph::Graph(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool)
    : thread_pool_(thread_pool ?
        thread_pool :
        std::make_shared<thread_pool::ThreadPool>(1)) {
}

struct Sortkey{
  inline bool operator() (const std::unique_ptr<biosoup::NucleicAcid>& struct1, const std::unique_ptr<biosoup::NucleicAcid>& struct2){
    return (struct1->name.compare(struct2->name)<0) ? true : false; 
  }
};

void Graph::Construct(std::vector<std::unique_ptr<biosoup::NucleicAcid>>& targets ,std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences) {
  //paramters for RAM 
  uint32_t k=21, w=11, bandwidth = 100, chain = 2, matches = 25, gap =100; 
  bool is_ava = true ; 
  double frequency = 0.01; 
  uint32_t LONG_READ = 0.8*2500; //0.8*read_length 

  //sort the sequence first 

  std::cerr << "[tarantula::Construct] Sorting: " << sequences.size() << std::endl;
  std::sort(sequences.begin(), sequences.end(), Sortkey());

  // map hic_reads to assembly
  ram::MinimizerEngine minimizer_engine{
    thread_pool_,
    k,
    w,
    bandwidth,
    chain,
    matches,
    gap};
  

  minimizer_engine.Minimize(targets.begin(), targets.end());
  minimizer_engine.Filter(frequency);

  std::uint64_t num_targets = biosoup::NucleicAcid::num_objects;
  biosoup::NucleicAcid::num_objects = 0; 
  int numPair=0; 

  if (sequences.empty()) {
    return;
  }

  std::vector<std::future<std::vector<biosoup::Overlap>>> futures;
  for (long unsigned int i=0; i< sequences.size(); i++) {
    if (i==0 && sequences[i]->name.compare(sequences[i+1]->name)!=0){
      continue; 
    }
    else if (i==sequences.size()-1 && sequences[i]->name.compare(sequences[i-1]->name)!=0){
      continue;
    }
    else if ( i!=0 && i!= sequences.size()-1 && sequences[i]->name.compare(sequences[i-1]->name)!=0 && sequences[i]->name.compare(sequences[i+1]->name)!=0){
      continue; 
    }

    numPair+=1;
    //do assembly here
    if (is_ava && sequences[i]->id >= num_targets) {
      continue;
    }

    futures.emplace_back(thread_pool_->Submit(
      [&] (const std::unique_ptr<biosoup::NucleicAcid>& sequence)
          -> std::vector<biosoup::Overlap> {
        return minimizer_engine.Map(sequence, false, false);
      },
      std::ref(sequences[i])));
  }
  std::cerr << "[tarantula::Construct] Number of pair: " << futures.size() << std::endl;

  std::unordered_map<uint32_t, std::vector<std::vector<biosoup::Overlap>>> readPairs;
  //match the pair and put it into unordered_map
  for (auto& it: futures){
    auto pair = it.get(); 
    if (pair.size()<1){
      //ignore as there are no overlaps
      continue; 
    }
    auto search = readPairs.find(pair[0].rhs_id);
    if (search == readPairs.end()){
      //not found
      std::vector<std::vector<biosoup::Overlap>> temp;
      temp.emplace_back(pair); 
      readPairs.insert(std::make_pair(pair[0].rhs_id, temp));
    }
    else{
      search->second.emplace_back(pair); 
    }
  }

  //go thru the pair to get rid of multiple overlap / no overlap
  std::vector<uint32_t> getRid; 
  for (auto const& rp : readPairs){
    if (rp.second.size()!=2){
      getRid.emplace_back(rp.first);
      //discard..., but shouldnt happen at all
    }
    
    for (auto const& ol: rp.second){
      if (ol.size()<1){
        //discard
        getRid.emplace_back(rp.first);
        break; 
      }
      if (ol.size()>1){
         auto numLR=0; 
         for (auto overlap: ol){
           if ((overlap.rhs_end - overlap.rhs_begin) > LONG_READ){
             numLR++;
           }
           if (numLR > 1){
             //if there are more than 1 long read discard
             getRid.emplace_back(rp.first);
             break;
           }
         }
         if (numLR==0){
           //dicard too because all are short reads
           getRid.emplace_back(rp.first);
         }
      }
    }
  }

  //for now is discard
  for (auto key: getRid){
    readPairs.erase(key); 
  }

  std::cerr << "[tarantula::Construct] Number of good read pair: " << readPairs.size() << std::endl;

  //then check if the first & second map to the same contig

  //then create piles
   

}

}  // namespace tarantula
