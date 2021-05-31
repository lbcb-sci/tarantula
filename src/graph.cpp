// Copyright (c) 2020 Robert Vaser

#include "graph.hpp"

#include <algorithm>
#include <fstream>
#include <stdexcept>

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
      unique_(),
      ambiguous_() {
}

void Graph::Construct(
    const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& targets,
    const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
    const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& pairs) {
  if (targets.empty() || sequences.empty() || sequences.size() != pairs.size()) {  // NOLINT
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

    timer.Start();

    std::vector<std::future<std::vector<Link>>> futures;
    for (std::size_t i = 0, j = 0, bytes = 0; i < sequences.size(); ++i) {
      bytes += sequences[i]->inflated_len + pairs[i]->inflated_len;
      futures.emplace_back(thread_pool_->Submit(
          [&] (std::uint32_t i) -> std::vector<Link> {
            std::vector<Link> dst;

            auto fo = minimizer_engine.Map(sequences[i], false, false);
            std::sort(fo.begin(), fo.end(),
                [] (const biosoup::Overlap& lhs,
                    const biosoup::Overlap& rhs) -> bool {
                  return lhs.rhs_end - lhs.rhs_begin > rhs.rhs_end - rhs.rhs_begin;  // NOLINT
                });
            for (std::size_t j = 0; j < fo.size(); ++j) {
              if (fo[j].rhs_end - fo[j].rhs_begin < 0.42 * sequences[i]->inflated_len) {  // NOLINT
                fo.resize(j);
                break;
              }
            }
            if (fo.empty()) {
              return dst;
            }

            auto so = minimizer_engine.Map(pairs[i], false, false);
            std::sort(so.begin(), so.end(),
                [] (const biosoup::Overlap& lhs,
                    const biosoup::Overlap& rhs) -> bool {
                  return lhs.rhs_end - lhs.rhs_begin > rhs.rhs_end - rhs.rhs_begin;  // NOLINT
                });
            for (std::size_t j = 0; j < so.size(); ++j) {
              if (so[j].rhs_end - so[j].rhs_begin < 0.42 * pairs[i]->inflated_len) {  // NOLINT
                so.resize(j);
                break;
              }
            }
            if (so.empty()) {
              return dst;
            }

            for (const auto& it : fo) {
              for (const auto& jt : so) {
                dst.emplace_back(
                    it.rhs_id, it.rhs_begin,
                    jt.rhs_id, jt.rhs_begin);
              }
            }
            return dst;
          },
          i));

      if (i != sequences.size() - 1 && bytes < (1ULL << 30)) {
        continue;
      }
      bytes = 0;

      for (auto& it : futures) {
        auto links = it.get();
        if (links.empty()) {
          continue;
        } else if (links.size() == 1) {
          if (unique_.capacity() == unique_.size()) {
            unique_.reserve(unique_.capacity() * 1.5);
          }
          unique_.emplace_back(links.front());
        } else {
          if (ambiguous_.capacity() <= ambiguous_.size() + links.size()) {
            ambiguous_.reserve(ambiguous_.capacity() * 1.5);
          }
          ambiguous_.insert(ambiguous_.end(), links.begin(), links.end());
        }
      }
      futures.clear();

      std::cerr << "[tarantula::Graph::Construct] mapped " << i - j + 1 << " sequence pairs "  // NOLINT
                << std::fixed << timer.Stop() << "s"
                << std::endl;
      j = i;

      timer.Start();
    }

    std::cerr << "[tarantula::Graph::Construct] stored " << unique_.size()
              << " unique links" << std::endl;
    std::cerr << "[tarantula::Graph::Construct] stored " << ambiguous_.size()
              << " ambiguous links" << std::endl;
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
