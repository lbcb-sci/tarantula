// Copyright (c) 2020 Robert Vaser

#ifndef TARANTULA_GRAPH_HPP_
#define TARANTULA_GRAPH_HPP_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "cereal/access.hpp"
#include "cereal/types/vector.hpp"
#include "thread_pool/thread_pool.hpp"

namespace tarantula {

class Graph {
 public:
  explicit Graph(
      bool checkpoints,
      std::shared_ptr<thread_pool::ThreadPool> thread_pool = nullptr);

  Graph(const Graph&) = delete;
  Graph& operator=(const Graph&) = delete;

  Graph(Graph&&) = default;
  Graph& operator=(Graph&&) = default;

  ~Graph() = default;

  int stage() const {
    return stage_;
  }

  void Construct(
      const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& targets,
      const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
      const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& pairs);

  void Scaffold();

  void Load();  // cereal load wrapper

  void Store() const;  // cereal store wrapper

 private:
  friend cereal::access;

  Graph() = default;  // needed for cereal

  template<class Archive>
  void serialize(Archive& archive) {  // NOLINT
    archive(stage_, checkpoints_, targets_, unique_, ambiguous_);
  }

  void CreateSubgraph(std::size_t i, std::size_t w);

  void CreateForceDirectedLayout(const std::string& path = "");

  struct Link {
   public:
    Link() = default;  // needed for cereal

    Link(
        std::uint16_t lhs_id, std::uint32_t lhs_pos,
        std::uint16_t rhs_id, std::uint32_t rhs_pos)
        : lhs_id(lhs_id),
          rhs_id(rhs_id),
          lhs_pos(lhs_pos),
          rhs_pos(rhs_pos) {}

    Link(const Link&) = default;
    Link& operator=(const Link&) = default;

    Link(Link&&) = default;
    Link& operator=(Link&&) = default;

    ~Link() = default;

    template<class Archive>
    void serialize(Archive& archive) {  // NOLINT
      archive(lhs_id, rhs_id, lhs_pos, rhs_pos);
    }

    std::uint16_t lhs_id;
    std::uint16_t rhs_id;
    std::uint32_t lhs_pos;
    std::uint32_t rhs_pos;
  };

  struct Node;
  struct Edge;

  struct Node {
   public:
    Node()
        : id(num_objects++) {}

    static std::uint32_t num_objects;

    std::uint32_t id;
    std::vector<Edge*> edges;
  };

  struct Edge {
   public:
    Edge()
        : id(num_objects++) {}

    static std::uint32_t num_objects;

    std::uint32_t id;
    std::uint32_t strength;
    double weight;
    Node* tail;
    Node* head;
  };

  std::shared_ptr<thread_pool::ThreadPool> thread_pool_;

  int stage_;
  bool checkpoints_;
  std::vector<std::uint32_t> targets_;
  std::vector<Link> unique_;
  std::vector<Link> ambiguous_;
  std::vector<std::unique_ptr<Node>> nodes_;
  std::vector<std::unique_ptr<Edge>> edges_;
};

}  // namespace tarantula

#endif  // TARANTULA_GRAPH_HPP_
