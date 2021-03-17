
#ifndef TARANTULA_PILEOGRAM_HPP_
#define TARANTULA_PILEOGRAM_HPP_

#include <cstdint>
#include <vector>

#include "biosoup/overlap.hpp"
#include "biosoup/nucleic_acid.hpp"

namespace tarantula {
class Pileogram{
  public:
    Pileogram(std::uint32_t id, std::uint32_t contig_len);
    Pileogram();

    Pileogram(const Pileogram&) = delete;
    Pileogram& operator=(const Pileogram&) = delete;

    Pileogram(Pileogram&&) = default;
    Pileogram& operator=(Pileogram&&) = default;

    ~Pileogram() = default;
    void addLayer(std::vector<biosoup::Overlap> overlap);
    void addLayer(uint32_t begin, uint32_t end);
    void initData();

    std::uint32_t id;
    std::uint32_t contig_len;
    std::uint32_t begin;
    std::uint32_t end;
    std::vector<std::uint64_t> data;
};

}  // namespace tarantula

#endif  // TARANTULA_PILEOGRAM_HPP_
