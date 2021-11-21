#ifndef _GRAPH_H_
#define _GRAPH_H_
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

// When you refuse to write a splay tree yourself
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>


template<typename T>
class Graph {
public:
  using sparse_vector = typename __gnu_pbds::tree<size_t, T, std::less<size_t>, __gnu_pbds::rb_tree_tag, __gnu_pbds::tree_order_statistics_node_update>;
  using const_sparse_vector = const __gnu_pbds::tree<size_t, T, std::less<size_t>, __gnu_pbds::rb_tree_tag, __gnu_pbds::tree_order_statistics_node_update>;
  using adj_matrix = typename std::map<size_t, sparse_vector>;
private:
  adj_matrix adj;

public:

  struct edge {
    size_t u, v;
    T w;

    edge(size_t u, size_t v, const T &w) : u{u}, v{v}, w{w} {}
  };

  class edge_container {
    typename adj_matrix::const_iterator _begin, _end;
  public:
    edge_container(const adj_matrix &adj): _begin{adj.begin()}, _end{adj.end()} {}
    class iterator {
      typename adj_matrix::const_iterator it_u;
      bool reset_v;
      typename sparse_vector::const_iterator it_v;

      public:
        iterator(typename adj_matrix::const_iterator it):
          it_u{it}, reset_v{true}, it_v{} {}

        iterator operator++() {
          if (reset_v) {
            it_v = it_u->second.begin();
            reset_v = false;
          }
          if (++it_v == it_u->second.end()) {
            ++it_u;
            reset_v = true;
          }
          return *this;
        }

        edge operator*() const {
          typename sparse_vector::const_iterator it;
          if (reset_v) {
            it = it_u->second.begin();
          } else {
            it = it_v;
          }
          return edge(it_u->first, it->first, it->second);
        }

        bool operator == (iterator other) const {
          if (it_u == other.it_u) {
            if (reset_v == other.reset_v || it_v == other.it_v) {
              return true;
            }
          }
          return false;
        }

        bool operator != (iterator other) const {
          return !(*this == other);
        }
    };

    iterator begin() const {
      return iterator(_begin);
    }

    iterator end() const {
      return iterator(_end);
    }
  };

  Graph(): adj{} {}

  size_t size() const {
    return adj.count();
  }

  void add_edge(size_t u, size_t v, const T &w) {
    adj[u][v] += w;
  }

  void add_edge(size_t u, size_t v, T &&w) {
    adj[u][v] += std::move(w);
  }

  void set_edge_weight(size_t u, size_t v, const T &w) {
    adj[u][v] = w;
  }

  void set_edge_weight(size_t u, size_t v, T &&w) {
    adj[u][v] = std::move(w);
  }

  void remove_edge(size_t u, size_t v) {
    if (has_edge(u, v)) {
      adj[u].erase(v);
      if (adj[u].empty()) {
        adj.erase(u);
      }
    } else {
      std::stringstream ss;
      ss << "edge (" << u << ", " << v << ") does not exist";
      throw std::domain_error(ss.str());
    }
  }

  bool has_edge(size_t u, size_t v) const {
    return adj.count(u) && adj[u].count(v);
  }

  const T &get_weight(size_t u, size_t v) const {
    if (has_edge(u, v)) {
      return adj[u][v];
    } else {
      std::stringstream ss;
      ss << "edge (" << u << ", " << v << ") does not exist";
      throw std::domain_error(ss.str());
    }
  }

  edge_container edges() const {
    return edge_container(adj);
  }

  size_t deg(size_t u) const {
    auto it = adj.find(u);
    if (it != adj.end()) {
      return it->second.size();
    } else {
      return 0;
    }
  }

  const_sparse_vector &neighbours(size_t u) const {
    return adj.at(u);
  }
};
#endif // _GRAPH_H_
