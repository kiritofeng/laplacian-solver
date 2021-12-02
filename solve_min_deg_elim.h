#ifndef __SOLVE_MIN_DEG_ELIM_H__
#define __SOLVE_MIN_DEG_ELIM_H__

#include <algorithm>
#include <map>
#include <set>
#include <tuple>
#include <vector>
#include <utility>

#include "graph.h"

namespace solve_min_deg_elim {
  template<typename T>
  void normalize(size_t n, std::vector<T> &x) {
    T avg = T(0);
    for (size_t u = 0; u < n; ++u) {
      avg += x[u];
    }
    if (avg != T(0)) { // is not orthogonal to all-ones vector, normalize
      avg /= T(n);
      for (size_t u = 0; u < n; ++u) {
        x[u] -= avg;
      }
    }
  }

  template<typename T>
  std::vector<T> solve_graphs(size_t n, const Graph<T> &G, const std::vector<T> &b, const T &eps) {
    std::vector<std::map<size_t, T>> adj(n);
    std::set<std::pair<size_t, size_t>> S;
    std::vector<T> nb(b), weighted_deg(n);
    std::vector<size_t> deg(n);
    for (size_t u = 0; u < n; ++u) {
      deg[u] = G.deg(u);
      for (auto e : G.neighbours(u)) {
        adj[u][e.first] = e.second;
        weighted_deg[u] += e.second;
      }
      S.insert({deg[u], u});
    }
    std::vector<std::tuple<size_t, size_t, T>> edges;
    while (!S.empty()) {
      size_t u = S.begin()->second;
      S.erase(S.begin());
      if (adj[u].empty()) {
        nb[u] = T(0);
      }
      for (auto e1 : adj[u]) {
        size_t v = e1.first;
        S.erase({deg[v], v});
        // push b to the other vertices, and establish clique
        T w = adj[v][u];
        for (auto e2 : adj[u]) {
          if (e2.first == v) {
            continue;
          }
          auto it = adj[v].find(e2.first);
          if (it == adj[v].end()) {
            deg[v] += 1;
          }
          adj[v][e2.first] += e2.second * w / weighted_deg[u];
          weighted_deg[v] += e2.second * w / weighted_deg[u];
        }
        nb[v] += nb[u] * w / weighted_deg[u];

        weighted_deg[v] -= w;
        deg[v] -= 1;
        adj[v].erase(u);

        edges.emplace_back(u, v, w);
        S.insert({deg[v], v});
      }
    }
    size_t last = n + 1;
    std::reverse(edges.begin(), edges.end());
    for (auto e : edges) {
      size_t u, v; T w;
      std::tie(u, v, w) = e;
      if (u != last) {
        if (last != n + 1) {
          nb[last] /= weighted_deg[last];
        }
        last = u;
      }
      nb[u] += nb[v] * w;
    }
    if (last != n + 1) {
      nb[last] /= weighted_deg[last];
    }
    normalize(n, nb);
    return nb;
  }
};

#endif // __SOLVE_MIN_DEG_ELIM_H__
