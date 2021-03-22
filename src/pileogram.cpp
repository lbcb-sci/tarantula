// Copyright (c) 2021 Cecilia Lee, Robert Vaser

#include <algorithm>

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

void Pileogram::AddLayer(std::vector<std::pair<std::uint32_t, std::uint32_t>>& overlaps) {
  std::vector<std::uint32_t> boundaries;
  for (const auto& it : overlaps) {
    boundaries.emplace_back((it.first + 1)  << 1);
    boundaries.emplace_back((it.second - 1) << 1 | 1);
  }
  std::sort(boundaries.begin(), boundaries.end());

  std::uint32_t coverage = 0;
  std::uint32_t last_boundary = 0;
  for (const auto& it : boundaries) {
    if (coverage > 0) {
      for (std::uint32_t i = last_boundary; i < (it >> 1); ++i) {
        data[i] += coverage;
      }
    }
    last_boundary = it >> 1;
    coverage += it & 1 ? -1 : 1;
  }
}


}  // namespace tarantula
