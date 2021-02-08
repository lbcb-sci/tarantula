// Copyright (c) 2021 Robert Vaser

#include <getopt.h>

#include <cstdint>
#include <iostream>
#include <vector>

#include "graph.hpp"

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/timer.hpp"

std::atomic<std::uint32_t> biosoup::Sequence::num_objects{0};

namespace {

static struct option options[] = {
  {"threads", required_argument, nullptr, 't'},
  {"version", no_argument, nullptr, 'v'},
  {"help", no_argument, nullptr, 'h'},
  {nullptr, 0, nullptr, 0}
};

std::unique_ptr<bioparser::Parser<biosoup::Sequence>> CreateParser(
    const std::string& path) {
  auto is_suffix = [] (const std::string& s, const std::string& suff) {
    return s.size() < suff.size() ? false :
        s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
  };

  if (is_suffix(path, ".fasta") || is_suffix(path, ".fasta.gz") ||
      is_suffix(path, ".fna")   || is_suffix(path, ".fna.gz") ||
      is_suffix(path, ".fa")    || is_suffix(path, ".fa.gz")) {
    try {
      return bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastaParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }
  if (is_suffix(path, ".fastq") || is_suffix(path, ".fastq.gz") ||
      is_suffix(path, ".fq")    || is_suffix(path, ".fq.gz")) {
    try {
      return bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastqParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }

  std::cerr << "[tarantula::CreateParser] error: file " << path
            << " has unsupported format extension (valid extensions: .fasta, "
            << ".fasta.gz, .fna, .fna.gz, .fa, .fa.gz, .fastq, .fastq.gz, "
            << ".fq, .fq.gz)"
            << std::endl;
  return nullptr;
}

void Help() {
  std::cout <<
      "usage: tarantula [options ...] <target> <sequences> [<sequences> ...]\n"
      "\n"
      "  # default output is to stdout\n"
      "  <target>/<sequences>\n"
      "    input file in FASTA/FASTQ format (can be compressed with gzip)\n"
      "\n"
      "  options:\n"
      "    -t, --threads <int>\n"
      "      default: 1\n"
      "      number of threads\n"
      "    --version\n"
      "      prints the version number\n"
      "    -h, --help\n"
      "      prints the usage\n";
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> input_paths;

  std::uint32_t num_threads = 1;

  std::string optstr = "t:h";
  int arg;
  while ((arg = getopt_long(argc, argv, optstr.c_str(), options, nullptr)) != -1) {  // NOLINT
    switch (arg) {
      case 't': num_threads = atoi(optarg); break;
      case 'v': std::cout << VERSION << std::endl; return 0;
      case 'h': Help(); return 0;
      default: return 1;
    }
  }

  if (argc == 1) {
    Help();
    return 0;
  }

  for (std::int32_t i = optind; i < argc; ++i) {
    input_paths.emplace_back(argv[i]);
  }

  if (input_paths.size() < 2) {
    std::cerr << "[tarantula::] error: missing input file!" << std::endl;
    return 1;
  }

  biosoup::Timer timer{};
  timer.Start();

  auto tparser = CreateParser(input_paths[0]);
  if (tparser == nullptr) {
    return 1;
  }
  std::vector<std::unique_ptr<biosoup::Sequence>> targets;
  try {
    targets = tparser->Parse(-1);
  } catch (std::invalid_argument& exception) {
    std::cerr << exception.what() << std::endl;
    return 1;
  }
  std::cerr << "[tarantula::] num targets = " << targets.size() << std::endl;

  std::vector<std::unique_ptr<biosoup::Sequence>> sequences;
  for (std::uint32_t i = 1; i < input_paths.size(); ++i) {
    auto sparser = CreateParser(input_paths[i]);
    if (sparser == nullptr) {
      return 1;
    }
    std::vector<std::unique_ptr<biosoup::Sequence>> part;
    try {
      part = sparser->Parse(-1);
    } catch (std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return 1;
    }
    sequences.insert(
        sequences.end(),
        std::make_move_iterator(part.begin()),
        std::make_move_iterator(part.end()));
  }
  std::cerr << "[tarantula::] num sequences = " << sequences.size() << std::endl;  // NOLINT

  auto thread_pool = std::make_shared<thread_pool::ThreadPool>(num_threads);

  tarantula::Graph graph{thread_pool};

  timer.Stop();
  std::cerr << "[tarantula::] " << std::fixed << timer.elapsed_time() << "s"
            << std::endl;

  return 0;
}
