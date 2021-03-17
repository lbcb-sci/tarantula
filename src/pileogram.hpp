// Copyright (c) 2021 Cecilia Lee, Robert Vaser

#ifndef TARANTULA_PILEOGRAM_HPP_
#define TARANTULA_PILEOGRAM_HPP_

#include <cstdint>
#include <vector>

#include "biosoup/overlap.hpp"
#include "cereal/types/vector.hpp"

namespace tarantula {

class Pileogram {
 public:
  Pileogram(std::uint32_t id, std::uint32_t contig_len);
  Pileogram() = default;

  Pileogram(const Pileogram&) = default;
  Pileogram& operator=(const Pileogram&) = default;

  Pileogram(Pileogram&&) = default;
  Pileogram& operator=(Pileogram&&) = default;

  ~Pileogram() = default;

  void AddLayer(std::vector<biosoup::Overlap> overlap);
  void AddLayer(std::uint32_t begin, std::uint32_t end);

  template<class Archive>
  void serialize(Archive& archive) {  // NOLINT
    archive(
        CEREAL_NVP(id),
        CEREAL_NVP(begin),
        CEREAL_NVP(end),
        CEREAL_NVP(data));
  }

  std::uint32_t id;
  std::uint32_t contig_len;
  std::uint32_t begin;
  std::uint32_t end;
  std::vector<std::uint64_t> data;
};

}  // namespace tarantula

#endif  // TARANTULA_PILEOGRAM_HPP_
