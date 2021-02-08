// Copyright (c) 2020 Robert Vaser

#include "graph.hpp"

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "gtest/gtest.h"

namespace tarantula {
namespace test {

class TarantulaTest: public ::testing::Test {
 public:
  void SetUp() {
    // auto p = bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastaParser>(TEST_DATA);  NOLINT
    g = Graph{nullptr};
  }

  Graph g;
};

TEST_F(TarantulaTest, Test) {
}

}  // namespace test
}  // namespace tarantula
