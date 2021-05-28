// Copyright (c) 2020 Robert Vaser

#include "graph.hpp"

#include <fstream>

#include "biosoup/timer.hpp"
#include "cereal/archives/binary.hpp"
#include "ram/minimizer_engine.hpp"

namespace tarantula {

Graph::Graph(
    bool checkpoints,
    std::shared_ptr<thread_pool::ThreadPool> thread_pool)
    : thread_pool_(thread_pool ?
        thread_pool :
        std::make_shared<thread_pool::ThreadPool>(1)),
      stage_(-1),
      checkpoints_(checkpoints),
      links_(),
      support_() {
}

void Graph::Construct(
    const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& targets,
    const std::string& path,
    const std::string& path_pair) {
  if (targets.empty() || path.empty() || path_pair.empty()) {
    return;
  }

  if (stage_ == -1 && checkpoints_) {  // checkpoint test
    Store();
  }

  biosoup::Timer timer{};

  if (stage_ == -1) {  // find links
    for (const auto& it : targets) {
      targets_.emplace_back(it->inflated_len);
    }

    timer.Start();

    ram::MinimizerEngine minimizer_engine{thread_pool_, 21, 11, 100, 2, 25, 100};  // NOLINT
    minimizer_engine.Minimize(targets.begin(), targets.end());
    minimizer_engine.Filter(0.0001);

    std::cerr << "[tarantula::Graph::Construct] minimized target sequences "
              << std::fixed << timer.Stop() << "s"
              << std::endl;
  }

  if (stage_ == -1) {  // checkpoint
    ++stage_;
    if (checkpoints_) {
      timer.Start();
      Store();
      std::cerr << "[tarantula::Graph::Construct] reached checkpoint "
                << std::fixed << timer.Stop() << "s"
                << std::endl;
    }
  }
}

void Graph::Store() const {
  std::ofstream os("tarantula.cereal");
  try {
    cereal::BinaryOutputArchive archive(os);
    archive(*this);
  } catch (std::exception&) {
    throw std::logic_error(
        "[tarantula::Graph::Store] error: unable to store archive");
  }
}

void Graph::Load() {
  std::ifstream is("tarantula.cereal");
  try {
    cereal::BinaryInputArchive archive(is);
    archive(*this);
  } catch (std::exception&) {
    throw std::logic_error(
        "[tarantula::Graph::Load] error: unable to load archive");
  }
}

}  // namespace tarantula
