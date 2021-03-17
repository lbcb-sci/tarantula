// Copyright (c) 2021 Cecilia Lee, Robert Vaser

#include "pileogram.hpp"

namespace tarantula {

Pileogram::Pileogram(std::uint32_t id, std::uint32_t contig_len)
  : id(id),
    contig_len(contig_len),
    begin(0),
    end(contig_len),
    data(contig_len, 0) {
}

void Pileogram::AddLayer(std::vector<biosoup::Overlap> overlap) {
  for (const auto& ol : overlap) {
    for (std::uint32_t i = ol.rhs_begin; i <= ol.rhs_end; i++) {
      data[i] += 1;
    }
  }
}

void Pileogram::AddLayer(std::uint32_t begin, std::uint32_t end) {
  for (std::uint32_t i = begin; i <= end; i++) {
    data[i] += 1;
  }
}

}  // namespace tarantula
