#include "pileogram.hpp"

namespace tarantula {

Pileogram::Pileogram(std::uint32_t id, std::uint32_t contig_len)
  : id(id),
    contig_len(contig_len),
    begin(0),
    end(contig_len),
    data(contig_len, 0) {
}

Pileogram::Pileogram() {
}

void Pileogram::addLayer(std::vector<biosoup::Overlap> overlap) {
  for (const auto& ol : overlap) {
    for (u_int32_t i = ol.rhs_begin; i <= ol.rhs_end; i++) {
      data[i] += 1;
    }
  }
}

void Pileogram::addLayer(uint32_t begin, uint32_t end) {
  for (uint32_t i = begin; i <= end; i++) {
    data[i] += 1;
  }
}
}  // namespace tarantula
