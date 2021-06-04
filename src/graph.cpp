// Copyright (c) 2020 Robert Vaser

#include "graph.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <fstream>
#include <random>
#include <stdexcept>
#include <unordered_set>

#include "biosoup/timer.hpp"
#include "cereal/archives/binary.hpp"
#include "ram/minimizer_engine.hpp"

namespace tarantula {

std::uint32_t Graph::Node::num_objects{0};
std::uint32_t Graph::Edge::num_objects{0};

Graph::Graph(
    bool checkpoints,
    std::shared_ptr<thread_pool::ThreadPool> thread_pool)
    : thread_pool_(thread_pool ?
        thread_pool :
        std::make_shared<thread_pool::ThreadPool>(1)),
      stage_(-1),
      checkpoints_(checkpoints),
      targets_(),
      unique_(),
      ambiguous_(),
      nodes_(),
      edges_() {
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
    biosoup::NucleicAcid::num_objects = 0;
    for (const auto& it : targets) {
      targets_.emplace_back(*it);
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

            auto map = [&] (const std::unique_ptr<biosoup::NucleicAcid>& na)
                -> std::vector<biosoup::Overlap> {
              auto ovl = minimizer_engine.Map(na, false, false);
              std::sort(ovl.begin(), ovl.end(),
                  [] (const biosoup::Overlap& lhs,
                      const biosoup::Overlap& rhs) -> bool {
                    return lhs.rhs_end - lhs.rhs_begin >
                           rhs.rhs_end - rhs.rhs_begin;
                  });
              for (std::size_t j = 0; j < ovl.size(); ++j) {
                if (ovl[j].rhs_end - ovl[j].rhs_begin < 0.42 * na->inflated_len) {  // NOLINT
                  ovl.resize(j);
                  break;
                }
              }
              return ovl;
            };

            auto first = map(sequences[i]);
            if (first.empty()) {
              return dst;
            }

            auto second = map(pairs[i]);
            if (second.empty()) {
              return dst;
            }

            for (const auto& it : first) {
              for (const auto& jt : second) {
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

void Graph::Scaffold() {
  for (std::size_t i = 0; i < targets_.size(); ++i) {
    CreateSubgraph(i, 25000);
    CreateForceDirectedLayout(std::to_string(i) + ".json");
  }
}

void Graph::CreateSubgraph(std::size_t i, std::size_t resolution) {
  nodes_.clear();
  Node::num_objects = 0;

  edges_.clear();
  Edge::num_objects = 0;

  if (i >= targets_.size() || targets_[i].inflated_len < resolution) {
    return;
  }

  auto occurence = [] (const std::string& str, const std::string& pattern)
      -> std::size_t {
    std::size_t dst = 0;
    std::string::size_type pos = 0;
    while ((pos = str.find(pattern, pos)) != std::string::npos) {
      ++dst;
      pos += pattern.size();
    }
    return dst;
  };

  std::vector<std::string> enzymes = {"GATC", "GAATC", "GACTC", "GAGTC", "GATTC"};  // NOLINT

  for (std::size_t j = 0; j < targets_[i].inflated_len; j += resolution) {
    nodes_.emplace_back(std::unique_ptr<Node>(new Node()));

    std::string window = targets_[i].InflateData(j, resolution);
    for (const auto& it : enzymes) {
      nodes_.back()->cuts += occurence(window, it);
    }
  }

  for (std::size_t j = 0; j < nodes_.size() - 1; ++j) {
    edges_.emplace_back(std::unique_ptr<Edge>(new Edge()));
    edges_.back()->tail = nodes_[j].get();
    edges_.back()->head = nodes_[j + 1].get();
    nodes_[j]->edges.emplace_back(edges_.back().get());

    edges_.emplace_back(std::unique_ptr<Edge>(new Edge()));
    edges_.back()->tail = nodes_[j + 1].get();
    edges_.back()->head = nodes_[j].get();
    nodes_[j + 1]->edges.emplace_back(edges_.back().get());

    edges_[edges_.size() - 2]->pair = edges_[edges_.size() - 1].get();
    edges_[edges_.size() - 1]->pair = edges_[edges_.size() - 2].get();
  }
  for (auto& it : edges_) {
    it->strength = 1;
  }

  auto add_links = [&] (const std::vector<Link>& links) -> void {
    for (const auto& it : links) {
      if (it.lhs_id == i && it.rhs_id == i) {
        auto tail_id = it.lhs_pos / resolution;
        auto head_id = it.rhs_pos / resolution;
        if (tail_id == head_id) {
          continue;
        }

        bool is_found = false;
        for (const auto& jt : nodes_[tail_id]->edges) {
          if (jt->head->id == head_id) {
            jt->weight++;
            jt->pair->weight++;
            is_found = true;
            break;
          }
        }

        if (is_found) {
          continue;
        }

        edges_.emplace_back(std::unique_ptr<Edge>(new Edge()));
        edges_.back()->tail = nodes_[tail_id].get();
        edges_.back()->head = nodes_[head_id].get();
        nodes_[tail_id]->edges.emplace_back(edges_.back().get());

        edges_.emplace_back(std::unique_ptr<Edge>(new Edge()));
        edges_.back()->tail = nodes_[head_id].get();
        edges_.back()->head = nodes_[tail_id].get();
        nodes_[head_id]->edges.emplace_back(edges_.back().get());

        edges_[edges_.size() - 2]->pair = edges_[edges_.size() - 1].get();
        edges_[edges_.size() - 1]->pair = edges_[edges_.size() - 2].get();
      }
    }
  };

  add_links(unique_);
  for (auto& it : edges_) {
    if (it->weight != 0) {
      it->strength += std::log(it->weight);
      it->weight = 0;
    }
  }
  add_links(ambiguous_);
  for (auto& it : edges_) {
    if (it->weight != 0) {
      it->strength += std::log(it->weight);
      it->weight = 0;
    }
  }

  std::vector<std::vector<double>> strengths(10);
  std::size_t bad = 0;
  for (std::size_t j = 0; j < edges_.size(); j += 2) {
    if (edges_[j]->strength != edges_[j]->pair->strength) {
      ++bad;
    }
    std::size_t distance = std::min(9, std::abs(
        static_cast<std::int32_t>(edges_[j]->head->id) -
        static_cast<std::int32_t>(edges_[j]->tail->id)));
    strengths[distance].emplace_back(edges_[j]->strength);
  }
  std::cerr << "Component " << i << " [" << nodes_.size() << "] " << std::endl;
  std::cerr << "Bad edges = " << bad << std::endl;
  for (std::size_t j = 0; j < strengths.size(); ++j) {
    std::cerr << j << ": [" << strengths[j].size() << "] ";
    if (strengths[j].empty()) {
      std::cerr << std::endl;
      continue;
    }
    std::sort(strengths[j].begin(), strengths[j].end());
    std::cerr << strengths[j].front() << " "
              << strengths[j][1 * strengths[j].size() / 4] << " "
              << strengths[j][2 * strengths[j].size() / 4] << " "
              << strengths[j][3 * strengths[j].size() / 4] << " "
              << strengths[j].back() << std::endl;
  }
  std::cerr << std::endl;
}

void Graph::CreateForceDirectedLayout(const std::string& path) {
  if (nodes_.empty() || edges_.empty()) {
    return;
  }

  std::ofstream os;
  if (!path.empty()) {
    os.open(path);
    os << "{" << std::endl;
  }

  static std::uint64_t seed = 21;
  seed <<= 1;

  std::mt19937 generator(seed);
  std::uniform_real_distribution<> distribution(0., 1.);

  struct Point {
   public:
    Point() = default;
    Point(double x, double y)
        : x(x),
          y(y) {}

    bool operator==(const Point& other) const {
      return x == other.x && y == other.y;
    }
    Point operator+(const Point& other) const {
      return Point(x + other.x, y + other.y);
    }
    Point& operator+=(const Point& other) {
      x += other.x;
      y += other.y;
      return *this;
    }
    Point operator-(const Point& other) const {
      return Point(x - other.x, y - other.y);
    }
    Point operator*(double c) const {
      return Point(x * c, y * c);
    }
    Point& operator/=(double c) {
      x /= c;
      y /= c;
      return *this;
    }
    double Norm() const {
      return std::sqrt(x * x + y * y);
    }

    double x;
    double y;
  };

  struct Quadtree {
   public:
    Quadtree(Point nucleus, double width)
        : nucleus(nucleus),
          width(width),
          center(0, 0),
          mass(0),
          subtrees() {
    }

    bool Add(const Point& p) {
      if (nucleus.x - width > p.x || p.x > nucleus.x + width ||
          nucleus.y - width > p.y || p.y > nucleus.y + width) {
        return false;
      }
      ++mass;
      if (mass == 1) {
        center = p;
      } else if (subtrees.empty()) {
        if (center == p) {
          return true;
        }
        double w = width / 2;
        subtrees.emplace_back(Point(nucleus.x + w, nucleus.y + w), w);
        subtrees.emplace_back(Point(nucleus.x - w, nucleus.y + w), w);
        subtrees.emplace_back(Point(nucleus.x - w, nucleus.y - w), w);
        subtrees.emplace_back(Point(nucleus.x + w, nucleus.y - w), w);
        for (auto& it : subtrees) {
          if (it.Add(center)) {
            break;
          }
        }
      }
      for (auto& it : subtrees) {
        if (it.Add(p)) {
          break;
        }
      }
      return true;
    }

    void Centre() {
      if (subtrees.empty()) {
        return;
      }
      center = Point(0, 0);
      for (auto& it : subtrees) {
        it.Centre();
        center += it.center * it.mass;
      }
      center /= mass;
    }

    Point Force(const Point& p, double k) const {
      auto delta = p - center;
      auto distance = delta.Norm();
      if (width * 2 / distance < 1) {
        return delta * (mass * (k * k) / (distance * distance));
      }
      delta = Point(0, 0);
      for (const auto& it : subtrees) {
        delta += it.Force(p, k);
      }
      return delta;
    }

    Point nucleus;
    double width;
    Point center;
    std::uint32_t mass;
    std::vector<Quadtree> subtrees;
  };

  std::vector<std::unordered_set<std::uint32_t>> components;
  std::vector<bool> is_visited(nodes_.size(), false);
  for (std::uint32_t i = 0; i < nodes_.size(); ++i) {
    if (nodes_[i] == nullptr || is_visited[i]) {
      continue;
    }

    components.resize(components.size() + 1);

    std::deque<std::uint32_t> que = { i };
    while (!que.empty()) {
      std::uint32_t j = que.front();
      que.pop_front();

      if (is_visited[j]) {
        continue;
      }
      const auto& node = nodes_[j];
      is_visited[node->id] = true;
      components.back().emplace(node->id);

      for (auto it : node->edges) {
        que.emplace_back(it->head->id);
      }
    }
  }
  std::vector<bool>().swap(is_visited);

  std::uint32_t c = 0;
  for (const auto& component : components) {
    if (component.size() < 6) {
      continue;
    }

    std::uint32_t num_iterations = 1000;
    double k = std::sqrt(1. / static_cast<double>(component.size()));
    double t = 0.1;
    double dt = t / static_cast<double>(num_iterations + 1);

    std::vector<Point> points(nodes_.size());
    for (const auto& it : component) {
      points[it].x = distribution(generator);
      points[it].y = distribution(generator);
    }

    for (std::uint32_t i = 0; i < num_iterations; ++i) {
      Point x = {0, 0}, y = {0, 0};
      for (const auto& n : component) {
        x.x = std::min(x.x, points[n].x);
        x.y = std::max(x.y, points[n].x);
        y.x = std::min(y.x, points[n].y);
        y.y = std::max(y.y, points[n].y);
      }
      double w = (x.y - x.x) / 2, h = (y.y - y.x) / 2;

      Quadtree tree(Point(x.x + w, y.x + h), std::max(w, h) + 0.01);
      for (const auto& n : component) {
        tree.Add(points[n]);
      }
      tree.Centre();

      std::vector<std::future<void>> thread_futures;
      std::vector<Point> displacements(nodes_.size(), Point(0, 0));

      auto thread_task = [&](std::uint32_t n) -> void {
        auto displacement = tree.Force(points[n], k);
        for (auto e : nodes_[n]->edges) {
          auto m = e->head->id;
          auto delta = points[n] - points[m];
          auto distance = delta.Norm();
          if (distance < 0.01) {
            distance = 0.01;
          }
          displacement += delta * (-1. * distance * e->strength / k);
        }
        auto length = displacement.Norm();
        if (length < 0.01) {
          length = 0.1;
        }
        displacements[n] = displacement * (t / length);
        return;
      };

      for (const auto& n : component) {
        thread_futures.emplace_back(thread_pool_->Submit(thread_task, n));
      }
      for (const auto& it : thread_futures) {
        it.wait();
      }
      for (const auto& n : component) {
        points[n] += displacements[n];
      }

      t -= dt;
    }

    for (const auto& it : edges_) {
      if (it == nullptr || (it->tail->id ^ it->head->id) != 1) {
        continue;
      }
      auto n = it->tail->id;
      auto m = it->head->id;

      if (component.find(n) != component.end() &&
          component.find(m) != component.end()) {
        it->weight = (points[n] - points[m]).Norm();
      }
    }

    if (!path.empty()) {
      os << "    \"component_" << c++ << "\": {" << std::endl;

      bool is_first_node = true;
      os << "      \"nodes\": {" << std::endl;
      for (const auto& it : component) {
        if (!is_first_node) {
          os << "," << std::endl;
        }
        is_first_node = false;
        os << "        \"" << it << "\": [";
        os << points[it].x << ", ";
        os << points[it].y << "]";
      }
      os << std::endl << "      }," << std::endl;

      bool is_first_edge = true;
      os << "      \"edges\": [" << std::endl;
      std::uint32_t nmin = -1, nmax = 0;
      for (const auto& it : component) {
        nmin = std::min(nmin, it);
        nmax = std::max(nmax, it);
      }
      for (std::uint32_t i = nmin; i < nmax; ++i) {
        if (!is_first_edge) {
          os << "," << std::endl;
        }
        is_first_edge = false;
        os << "        [\"" << i << "\", \"" << i + 1 << "\"]";
      }

      os << std::endl << "      ]" << std::endl;
      os << "    }";
    }
  }

  if (!path.empty()) {
    os << std::endl << "}";
    os << std::endl;
    os.close();
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
