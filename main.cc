#include <iostream>
#include <vector>

#include "graph.h"
#include "solve_sim_walk.h"

static const double EPS = 1e-6;

int main() {
  std::cin.tie(0)->sync_with_stdio(0);
  size_t n, m;
  std::cin >> n >> m;
  Graph<double> G;
  std::vector<double> b(n);
  for (size_t i = 0; i < m; ++i) {
    size_t u, v;
    double w;
    std::cin >> u >> v >> w;
    G.add_edge(u, v, w);
    G.add_edge(v, u, w);
  }
  for (size_t i = 0; i < n; ++i) {
    std::cin >> b[i];
  }
  auto x = solve_sim_walk::solve(n, G, b, EPS);
  for (auto xi : x) {
    std::cout << xi << std::endl;
  }
}
