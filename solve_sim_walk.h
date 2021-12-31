#ifndef _SOLVE_SIM_WALK_H_
#define _SOLVE_SIM_WALK_H_

#include <algorithm>
#include <climits>
#include <cmath>
#include <map>
#include <queue>
#include <random>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include "graph.h"
#include "solve_min_deg_elim.h"

/**
 * Laplacian solver using simulated walks, removing 1/10 the vertices each time
 *
 * Probably very very slow
 */

namespace solve_sim_walk {
  const size_t GE_LIMIT = 1e4;
  const unsigned RANDOM_WALK_ATTEMPTS = 10;
  const unsigned RANDOM_WALK_ATTEMPT_LENGTH = 30;
  const double RHO_ADJUSTMENT_FACTOR = 0.01;
  const double ITERATIVE_ADJUSTMENT_FACTOR = 0.001;

  template<typename T>
  using adj_list = std::vector<std::vector<std::pair<size_t, T>>>;

  template<typename T>
  std::tuple<bool, size_t, T> sim_random_walk(
      size_t, const adj_list<T>&, size_t, const std::vector<bool>&, unsigned=UINT_MAX);

  template<typename T>
  std::vector<std::tuple<bool, size_t, T>> sim_random_walks(
      size_t n, const Graph<T> &G, size_t u, const std::vector<bool> &C, unsigned walks, unsigned max_len=UINT_MAX);

  template<typename T>
  std::vector<T> solve_iterative(size_t n, const Graph<T> &G, std::vector<T> &b, const T &eps, unsigned iterations);

  template<typename T>
  std::vector<T> solve(size_t n, const Graph<T> &G, std::vector<T> &b, const T &eps, bool normalize=true);

  template<typename T>
  std::vector<T> solve_vertex_elimination(
      size_t n, const Graph<T> &G, std::vector<T> &b, const T &eps);

  template<typename T>
  std::vector<T> solve_gaussian(size_t n, const Graph<T> &G, const std::vector<T> &b, const T &eps);

  template<typename T>
  std::vector<std::vector<std::pair<size_t, T>>> graph_to_adj_list(const Graph<T> &G);

  template<typename T>
  adj_list<T> graph_to_adj_list(size_t n, const Graph<T> &G) {
    adj_list<T> adj(n);
    for (size_t i = 0; i < n; ++i) {
      if (G.deg(i) > 0) {
        for (auto e : G.neighbours(i)) {
          adj[i].push_back(e);
        }
      }
    }
    return adj;
  }

  template<typename T>
  std::tuple<bool, size_t, T> sim_random_walk(
      size_t n, const adj_list<T> &adj, size_t u, const std::vector<bool> &C, unsigned max_len) {
    static std::mt19937 gen{0xdeadbeef};

    T dist = T(0);
    size_t len = 0;
    while (!C[u] && len < max_len) {
      size_t idx = std::uniform_int_distribution<size_t>(0, adj[u].size() - 1)(gen);

      T weight;
      std::tie(u, weight) = adj[u][idx];

      dist += 1 / weight;
      len += 1;
    }

    if (!C[u]) {
      return std::make_tuple(false, 0, T(0));
    } else {
      return std::make_tuple(true, u, dist);
    }
  }

  template<typename T>
  std::vector<std::tuple<bool, size_t, T>> sim_random_walks(
      size_t n, const adj_list<T> &adj, size_t u, const std::vector<bool> &C, unsigned walks, unsigned max_len) {
    static std::mt19937 gen{0xdeadbeef};

    // cache degrees; lookup is surprisingly expensive
    std::vector<size_t> deg(n);
    for (size_t i = 0; i < n; ++i) {
      deg[i] = adj[i].size();
    }

    std::vector<std::tuple<bool, size_t, T>> ret;
    // cur, cnt, len, dist
    std::queue<std::tuple<size_t, unsigned, unsigned, T>> queue;
    queue.emplace(u, walks, 0, 0);
    while (!queue.empty()) {
      size_t cur;
      unsigned long cnt, len;
      T dist;
      std::tie(cur, cnt, len, dist) = queue.front();
      queue.pop();
      if (C[cur] || len >= max_len) {
        for (unsigned i = 0; i < cnt; ++i) {
          ret.emplace_back(C[cur], cur, dist);
        }
        continue;
      }
      std::vector<size_t> dest_idx;
      for (size_t i = 0; i < cnt; ++i) {
        size_t idx = std::uniform_int_distribution<size_t>(0, deg[cur] - 1)(gen);
        dest_idx.push_back(idx);
      }

      std::sort(dest_idx.begin(), dest_idx.end());
      const size_t deg_cur = deg[cur];
      size_t last = deg_cur; // sentinel value
      unsigned long last_cnt = 0;
      for (size_t idx : dest_idx) {
        if (idx != last) {
          if (last != deg_cur) {
            size_t v; T weight;
            std::tie(v, weight) = adj[cur][last];
            queue.emplace(v, last_cnt, len + 1, dist + T(1) / weight);
          }
          last = idx;
          last_cnt = 0;
        }
      }
      size_t v; T weight;
      std::tie(v, weight) = adj[cur][last];
      queue.emplace(v, last_cnt, len + 1, dist + weight);
    }
    return ret;
  }

  template<typename T>
  std::vector<T> solve_iterative(size_t n, const Graph<T> &G, std::vector<T> &b,
      const T &eps, unsigned iterations) {
    std::vector<T> x = solve(n, G, b, eps, 1);
    for (size_t i = 0; i < iterations; ++i) {
      std::vector<T> x_diff(n);
      T max_diff = T(0);
      for (size_t u = 0; u < n; ++u) {
        x_diff[u] = b[u];
        T weighted_degree{0};
        for (auto e : G.neighbours(u)) {
          x_diff[u] += b[e.first] * e.second;
          weighted_degree += e.second;
        }
        x_diff[u] -= b[u] * weighted_degree;
        max_diff = std::max(std::abs(x_diff[u]), max_diff);
      }
      if (max_diff < eps) {
        break;
      }
      std::vector<T> x_delta = solve(n, G, x_diff, eps, 1);
      for (size_t j = 0; j < n; ++j) {
        x[j] += T(ITERATIVE_ADJUSTMENT_FACTOR) * x_delta[j];
      }
    }
    return x;
  }

  template<typename T>
  std::vector<T> solve(size_t n, const Graph<T> &G, std::vector<T> &b, const T &eps, bool normalize) {
    std::vector<T> x;
    if (n <= GE_LIMIT) {
      x = solve_min_deg_elim::solve_graphs(n, G, b, eps);
    } else {
      x = solve_vertex_elimination(n, G, b, eps);
    }
    if (normalize) {
      T avg = T(0);
      for (size_t u = 0; u < n; ++u) {
        avg += x[u];
      }
      avg /= T(n);
      for (size_t u = 0; u < n; ++u) {
        x[u] -= avg;
      }
    }
    return x;
  }

  std::vector<size_t> sample(size_t n, size_t k) {
    static std::mt19937 gen{0xdeadbeef};
    std::vector<size_t> vals(n);
    for (size_t i = 0; i < n; ++i) {
      vals[i] = i;
    }
    std::shuffle(vals.begin(), vals.end(), gen);
    return std::vector<size_t>(vals.begin(), vals.begin() + k);
  }

  template<typename T>
  std::vector<size_t> get_largest_degree(size_t n, const Graph<T> &G, size_t k) {
    std::vector<size_t> vals(n), deg(n);
    for (size_t i = 0; i < n; ++i) {
      vals[i] = i;
      deg[i] = G.deg(i);
    }
    std::sort(vals.begin(), vals.end(), [=](size_t a, size_t b) {
        return deg[a] > deg[b];
    });
    return std::vector<size_t>(vals.begin(), vals.begin() + k);
  }

  template<typename T>
  std::vector<T> solve_vertex_elimination(
      size_t n, const Graph<T> &G, std::vector<T> &b, const T &eps) {

    const unsigned rho = std::round(RHO_ADJUSTMENT_FACTOR / (eps * eps)) + 1;

    Graph<T> H;
    std::vector<T> b2;
    std::vector<bool> in_C(n, false);
    // let C be 10% of vertices with largest degree
    std::vector<size_t> C = get_largest_degree(n, G, n / 10);
    for (size_t c : C) {
      in_C[c] = 1;
    }
    adj_list<T> adj_list = graph_to_adj_list(G);
    for (auto e : G.edges()) {
      if (in_C[e.u]) {
        continue;
      }
      unsigned cnt = 0;
      for (auto a : sim_random_walks(n, adj_list, e.u, in_C, RANDOM_WALK_ATTEMPTS, RANDOM_WALK_ATTEMPT_LENGTH)) {
        if (!std::get<0>(a)) {
          ++cnt;
        }
      }
      if (cnt >= 0.5 * RANDOM_WALK_ATTEMPTS) {
        in_C[e.u] = true;
        C.push_back(e.u);
      }
    }
    std::sort(C.begin(), C.end());
    // std::cerr << "|C|/|V| = " << 1.0 * C.size() / n << std::endl;
    // F is the remainder
    std::vector<size_t> F;
    for (size_t i = 0; i < n; ++i) {
      if (in_C[i]) {
        continue;
      }
      F.push_back(i);
    }
    std::vector<size_t> compression_map(n, n + 1);
    for (size_t i = 0; i < C.size(); ++i) {
      compression_map[C[i]] = i;
    }
    for (size_t i = 0; i < F.size(); ++i) {
      compression_map[F[i]] = i;
    }
    std::vector<std::map<size_t, T>> H_adj(n);
    for (auto e : G.edges()) {
      if (e.u > e.v) { // assume undirected for now
        continue;
      }
      auto u_walks = sim_random_walks(n, adj_list, e.u, in_C, rho);
      auto v_walks = sim_random_walks(n, adj_list, e.v, in_C, rho);
      for (size_t i = 0; i < u_walks.size(); ++i) {
        T H_r = T(1) / e.w;
        bool _;
        size_t u_walk_last; T u_walk_dist;
        size_t v_walk_last; T v_walk_dist;
        std::tie(_, u_walk_last, u_walk_dist) = u_walks[i];
        std::tie(_, v_walk_last, v_walk_dist) = v_walks[i];

        H_r += u_walk_dist;
        H_r += v_walk_dist;

        size_t nu = compression_map[u_walk_last];
        size_t nv = compression_map[v_walk_last];
        assert(std::binary_search(C.begin(), C.end(), u_walk_last));
        assert(std::binary_search(C.begin(), C.end(), v_walk_last));
        T nw = T(rho) / H_r;
        H_adj[nu][nv] += nw;
        H_adj[nv][nu] += nw;
      }
    }
    for (size_t u = 0; u < n; ++u) {
      for (auto e : H_adj[u]) {
        H.add_edge(u, e.first, T(1) / e.second);
      }
    }
    // populate nb for C
    std::vector<T> C_b(C.size());
    for (size_t f : F) {
      for (auto a : sim_random_walks(n, G, f, in_C, rho)) {
        size_t dest = std::get<1>(a);
        C_b[dest] += b[dest] / rho;
      }
    }
    std::vector<T> solved_x(n);
    // solve C recursively
    auto solved_C_x = solve(C.size(), H, C_b, eps, false);
    for (size_t i = 0; i < C.size(); ++i) {
      solved_x[C[i]] = solved_C_x[i];
    }
    // solve the remainder
    std::vector<T> F_b(F.size());
    for (size_t f : F) {
      for (auto a : sim_random_walks(n, adj_list, f, in_C, rho)) {
        solved_x[f] += solved_x[std::get<1>(a)] / rho;
      }
      T weighted_degree(0);
      for (auto e : G.neighbours(f)) {
        weighted_degree += e.second;
      }
      solved_x[f] += b[f] / weighted_degree;
    }
    return solved_x;
  }

  template<typename T>
  std::vector<T> solve_gaussian(size_t n, const Graph<T> &G,
      const std::vector<T> &b, const T &eps) {
    std::vector<std::vector<T>> L(n);
    std::vector<T> nb(b);
    for (auto & v : L) {
      v.resize(n);
    }
    for (auto e : G.edges()) {
      L[e.u][e.v] = -e.w;
      L[e.u][e.u] += e.w;
    }
    // Gaussian time
    for (size_t i = 0; i < n - 1; ++i) {
      size_t idx = i;
      for (; idx < n; ++idx) {
        if (L[idx][i] != T(0)) {
          break;
        }
      }
      if (idx == n) {
        throw std::runtime_error("matrix has rank less than n - 1");
      } else if (idx > i) {
        std::swap(L[i], L[idx]);
        std::swap(nb[i], nb[idx]);
      }
      for (size_t j = i + 1; j < n; ++j) {
        L[i][j] /= L[i][i];
      }
      nb[i] /= L[i][i];
      L[i][i] = T(1);
      for (size_t j = idx + 1; j < n; ++j) {
        if (L[j][i] == T(0)) {
          continue;
        }
        for (size_t k = i + 1; k < n; ++k) {
          L[j][k] -= L[i][k] * L[j][i] / L[i][i];
        }
        nb[j] -= nb[i] * L[j][i] / L[i][i];
        L[j][i] = T(0);
      }
    }
    // last row is all zeros, set to 1
    for (size_t i = n - 2; ; --i) {
      for (size_t j = i + 1; j < n; ++j) {
        nb[i] -= L[i][j] * nb[j];
      }
      if (i == 0) break;
    }
    return nb;
  }

};

#endif // _SOLVE_SIM_WALK_H_
