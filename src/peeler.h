#ifndef PEELER_H
#define PEELER_H

#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <mpi.h>

struct Edge {
    int u, v;
    int cycle_support;
    int flow_support;
};

struct Graph {
    std::vector<Edge> edges;
};

void compute_supports(std::vector<Graph>& partitions, int rank);
std::tuple<int, int> distributedMaxTrussNumbers(std::vector<Graph>& partitions, int rank);
void update_cycle_support(Edge& edge, int new_support);
void update_flow_support(Edge& edge, int new_support);
void UpSupc(Edge& e_star, std::unordered_map<int, std::unordered_set<int>>& out_neighbors, std::unordered_map<int, std::unordered_set<int>>& in_neighbors, std::unordered_map<int, int>& LinkSta);
void UpSupf(Edge& e_star, std::unordered_map<int, std::unordered_set<int>>& out_neighbors, std::unordered_map<int, std::unordered_set<int>>& in_neighbors, std::unordered_map<int, int>& LinkSta);
void perform_flow_decomposition(std::unordered_set<Edge>& edges, std::unordered_map<int, std::unordered_set<int>>& out_neighbors, std::unordered_map<int, std::unordered_set<int>>& in_neighbors, std::unordered_map<int, int>& LinkSta);
void stratified_local_peel_processing(std::unordered_set<Edge>& F_i, const std::unordered_map<Edge, int>& cycle_supports, const std::unordered_map<Edge, int>& flow_supports, int k_c, int k_f, std::unordered_map<int, std::unordered_set<int>>& out_neighbors, std::unordered_map<int, std::unordered_set<int>>& in_neighbors, std::unordered_map<int, int>& LinkSta);

#endif // PEELER_H