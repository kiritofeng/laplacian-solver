#include <iostream>
#include <vector>

#include "graph.h"
#include "solve_min_deg_elim.h"

const double EPS = 1e-25;

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
    G.add_edge(u - 1, v - 1, w);
    G.add_edge(v - 1, u - 1, w);
  }
  for (size_t i = 0; i < n; ++i) {
    std::cin >> b[i];
  }
  auto x = solve_min_deg_elim::solve_graphs(n, G, b, EPS);
  for (auto xi : x) {
    std::cout << xi << std::endl;
  }
}
