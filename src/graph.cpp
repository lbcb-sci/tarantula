// Copyright (c) 2020 Robert Vaser

#include <iostream>

#include "graph.hpp"




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
  double frequency = 0.01; 

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
  
  int numPair=0; 

  if (sequences.empty()) {
    return;
  }

  std::vector<std::future<void>> void_futures; 
  std::unordered_map<std::string, std::vector<std::vector<biosoup::Overlap>>> readPairs;
  for (long unsigned int i=0; i< sequences.size()-1; i++) {
    if (sequences[i]->name.compare(sequences[i+1]->name)==0){
      numPair+=1;
      //pair
      Process(void_futures, minimizer_engine, sequences[i], sequences[i+1], readPairs); 
      i++; 
    }
    else{
      continue; 
    }
  }

  for (const auto& it : void_futures) {
    it.wait();
  }

  std::cerr << "[tarantula::Construct] Number of good read pair: " << readPairs.size() << std::endl;

  for (auto const& rp: readPairs){
    if (rp.second[0].size()!=1){
      std::cerr << "error" << std::endl; 
    }
    if (rp.second[0].size()!=1){
      std::cerr << "error" << std::endl; 
    }
  }

  //then create pilo-o-gram per contig, qn is how to get all the contigs?
  CreateGraph(targets); 
  std::cerr << "[tarantula::Construct] Graph created, number of nodes: " << contigs.size() << std::endl;
  FillPileogram(readPairs); 
  std::cerr << "[tarantula::Construct] Pile-o-gram created, number of nodes: " << contigs.size() << std::endl;
  
  return; 
}

//technically this also can be multi-thread - just remove the for loop
void Graph::CreateGraph(std::vector<std::unique_ptr<biosoup::NucleicAcid>>& targets){
  //go thru all the read pair and then create node?
  for (auto const& target: targets){ 
    Node node = Node(target->id, target->inflated_len); 
    contigs.insert({target->id,node}); 
    //everything init except for the pilogram data part --> that one need the readpair to update
  }
}

void Graph::FillPileogram(std::unordered_map<std::string, std::vector<std::vector<biosoup::Overlap>>>& readPairs){
  //go through all the read pairs
  std::unordered_map<std::uint32_t,Node>::iterator found;
  std::tuple<std::uint32_t, std::uint32_t> overlap;
  //find the contig, then add layer to pileogram
  for (const auto& rp: readPairs){
    found = contigs.find(rp.second[0][0].rhs_id); 
    if (found==contigs.end()){
      //not found 
      std::cerr << "ERROR contig not found" << std::endl;
    }
    else{
      overlap = GetOverlap(rp.second[0][0],rp.second[1][0]); 
      found->second.pileogram.addLayer(std::get<0>(overlap), std::get<1>(overlap)); 
      std::cerr << "[tarantula::Construct] Fill pile-o-gram - contig: " << rp.second[0][0].rhs_id << ", length: " << found->second.pileogram.contigLen << "| begin: " << std::get<0>(overlap) << " end: " << std::get<1>(overlap) << std::endl;
    }
  }
}

std::tuple<std::uint32_t, std::uint32_t> Graph::GetOverlap(biosoup::Overlap ol1, biosoup::Overlap ol2){
  std::tuple<std::uint32_t, std::uint32_t> overlap;
  if (ol1.rhs_end > ol2.rhs_begin){
    overlap = std::make_tuple(ol1.rhs_begin, ol2.rhs_end); 
  }
  else{
    overlap = std::make_tuple(ol2.rhs_begin, ol1.rhs_end); 
  }
  return overlap; 
}

void Graph::Process(std::vector<std::future<void>>& futures, ram::MinimizerEngine& minimizer_engine,std::unique_ptr<biosoup::NucleicAcid>& sequence1, std::unique_ptr<biosoup::NucleicAcid>& sequence2, std::unordered_map<std::string, std::vector<std::vector<biosoup::Overlap>>>& readPairs){
  futures.emplace_back(thread_pool_->Submit(
    [&] (const std::unique_ptr<biosoup::NucleicAcid>& sequence1,const std::unique_ptr<biosoup::NucleicAcid>& sequence2)
        -> void {
          std::vector<std::vector<biosoup::Overlap>> minimizerResult; 
          minimizerResult.emplace_back(minimizer_engine.Map(sequence1, false, false)); 
          minimizerResult.emplace_back(minimizer_engine.Map(sequence2, false, false)); 
          auto longReadLen = sequence1->inflated_len*0.8; 

          if (minimizerResult[0].size()<1 || minimizerResult[1].size()<1 ){
            //no overlap, discard 
            return; 
          }

          for (auto& rp: minimizerResult){
            if (rp.size()==1){
              continue;
            }
            std::vector<int> discard; 
            biosoup::Overlap temp; 
            auto numLR=0; 
            for (auto const& ol: rp){
              if (ol.rhs_end-ol.rhs_begin > longReadLen){
                numLR++;
                temp = ol; 
              }
              if (numLR > 1){
                //if there are more than 1 long read discard
                return; 
              }
            }
            if (numLR==0){
              //dicard too because all are short reads
              return; 
            }
            rp.clear(); 
            rp.emplace_back(temp); 
          }

          if (minimizerResult[0][0].rhs_id!=minimizerResult[1][0].rhs_id){
            return; 
          }

          //add into the map
          readPairs.insert(std::make_pair(sequence1->name, minimizerResult)); 
          return; 
    },
    std::ref(sequence1), std::ref(sequence2)));
}

}  // namespace tarantula
