#ifndef _SOLVE_SIM_WALK_H_
#define _SOLVE_SIM_WALK_H_

#include <algorithm>
#include <cassert>
#include <climits>
#include <iostream>
#include <random>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include "graph.h"

/**
 * Laplacian solver using simulated walks, removing half the vertices at a time.
 *
 * Probably very very slow
 */

namespace solve_sim_walk {
  const size_t GE_LIMIT = 1500;
  const size_t RANDOM_WALK_ATTEMPTS = 10;
  const size_t RANDOM_WALK_ATTEMPT_LENGTH = 30;
  const double RHO_ADJUSTMENT_FACTOR = 0.1;
  const double ITERATIVE_ADJUSTMENT_FACTOR = 0.001;

  template<typename T>
  std::tuple<bool, size_t, T> sim_random_walk(
      const Graph<T> &G, size_t u, std::vector<bool> &C, size_t max_len=UINT_MAX);

  template<typename T>
  std::vector<T> solve_iterative(size_t n, const Graph<T> &G, std::vector<T> &b, const T &eps, size_t iterations);

  template<typename T>
  std::vector<T> solve(size_t n, const Graph<T> &G, std::vector<T> &b, const T &eps, bool normalize=true);

  template<typename T>
  std::vector<T> solve_vertex_elimination(
      size_t n, const Graph<T> &G, std::vector<T> &b, const T &eps);

  template<typename T>
  std::vector<T> solve_gaussian(size_t n, const Graph<T> &G, const std::vector<T> &b, const T &eps);

  template<typename T>
  std::tuple<bool, size_t, T> sim_random_walk(
      const Graph<T> &G, size_t u, std::vector<bool> &C, size_t max_len) {
    static std::mt19937 gen{0xdeadbeef};

    assert(G.deg(u) > 0);
    T dist = T(0);
    size_t len = 0;
    while (!C[u] && len < max_len) {
      size_t idx = std::uniform_int_distribution<size_t>(0, G.deg(u) - 1)(gen);

      T weight;
      std::tie(u, weight) = *G.neighbours(u).find_by_order(idx);

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
  std::vector<T> solve_iterative(size_t n, const Graph<T> &G, std::vector<T> &b,
      const T &eps, size_t iterations) {
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
    if (n < GE_LIMIT) {
      x = solve_gaussian(n, G, b, eps);
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
    std::vector<size_t> vals(n);
    for (size_t i = 0; i < n; ++i) {
      vals[i] = i;
    }
    std::sort(vals.begin(), vals.end(), [=](size_t a, size_t b) {
        return G.deg(a) > G.deg(b);
    });
    for (size_t i = 0; i < vals.size(); ++i) {
      if (i < k) {
        std::cerr << "Keeping vertex of degree " << G.deg(vals[i]) << std::endl;
      } else {
        std::cerr << "Eliminating vertex of degree " << G.deg(vals[i]) << std::endl;
      }
    }
    return std::vector<size_t>(vals.begin(), vals.begin() + k);
  }

  template<typename T>
  std::vector<T> solve_vertex_elimination(
      size_t n, const Graph<T> &G, std::vector<T> &b, const T &eps) {

    T rho = RHO_ADJUSTMENT_FACTOR / (eps * eps);
    Graph<T> H1;
    std::vector<T> b2;
    std::vector<bool> in_C(n, false);
    // let C be 10% of vertices with largest degree
    std::vector<size_t> C = get_largest_degree(n, G, n / 10);
    for (size_t c : C) {
      in_C[c] = 1;
    }
    for (auto e : G.edges()) {
      if (in_C[e.u]) {
        continue;
      }
      size_t cnt = 0;
      for (size_t i = 0; i < RANDOM_WALK_ATTEMPTS; ++i) {
        if (!std::get<0>(sim_random_walk(G, e.u, in_C, RANDOM_WALK_ATTEMPT_LENGTH))) {
          cnt += 1;
        }
      }
      if (cnt / RANDOM_WALK_ATTEMPTS >= 0.9) {
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
    for (auto e : G.edges()) {
      if (e.u > e.v) { // assume undirected for now
        continue;
      }
      for (size_t i = 0; i < rho; ++i) {
        T H1_r = e.w;
        bool dummy;
        size_t u_walk_last; T u_walk_dist;
        size_t v_walk_last; T v_walk_dist;
        std::tie(dummy, u_walk_last, u_walk_dist) = sim_random_walk(G, e.u, in_C);
        std::tie(dummy, v_walk_last, v_walk_dist) = sim_random_walk(G, e.v, in_C);

        H1_r += u_walk_dist;
        H1_r += v_walk_dist;

        size_t nu = compression_map[u_walk_last];
        size_t nv = compression_map[v_walk_last];
        T nw = T(1) / H1_r;
        H1.add_edge(nu, nv, nw);
        H1.add_edge(nv, nu, nw);
      }
    }
    // populate nb for C
    std::vector<T> C_b(C.size());
    for (size_t i = 0; i < rho; ++i) {
      T H2_b = T(0);
      for (size_t f : F) {
        // walk 1 / rho fraction into C
        size_t dest = compression_map[std::get<1>(sim_random_walk(G, f, in_C))];
        C_b[dest] += b[dest] / rho;
      }
    }
    std::vector<T> solved_x(n);
    // solve C recursively
    auto solved_C_x = solve(C.size(), H1, C_b, eps, false);
    for (size_t i = 0; i < C.size(); ++i) {
      solved_x[C[i]] = solved_C_x[i];
    }
    // solve the remainder
    std::vector<T> F_b(F.size());
    for (size_t f : F) {
      for (size_t i = 0; i < rho; ++i) {
        solved_x[f] += solved_x[std::get<1>(sim_random_walk(G, f, in_C))] / rho;
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
