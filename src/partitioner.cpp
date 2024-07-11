#include "partitioner.h"
#include <iostream>
#include <unordered_map>
#include <queue>

using namespace std;

int score_cycle_triangles(const Edge& e, const unordered_set<Edge>& Fi, const unordered_map<int, unordered_set<int>>& out_neighbors, const unordered_map<int, unordered_set<int>>& in_neighbors) {
    int score = 0;
    for (int w : out_neighbors.at(e.v)) {
        if (in_neighbors.at(e.u).find(w) != in_neighbors.at(e.u).end()) {
            Edge e1 = {e.v, w};
            Edge e2 = {w, e.u};
            if (Fi.find(e1) != Fi.end() && Fi.find(e2) != Fi.end()) {
                score++;
            }
        }
    }
    return score;
}

vector<unordered_set<Edge>> type_aware_balanced_partitioner(const unordered_map<int, unordered_set<int>>& out_neighbors, const unordered_map<int, unordered_set<int>>& in_neighbors, int p, double epsilon, const string& triangle_type) {
    vector<unordered_set<Edge>> partitions(p);
    vector<int> C(p, 0);
    vector<unordered_set<Edge>> D(p);
    unordered_set<Edge> unassigned_edges;

    // Initialize unassigned_edges
    for (const auto& kv : out_neighbors) {
        int u = kv.first;
        for (int v : kv.second) {
            unassigned_edges.insert({u, v});
        }
    }

    // Calculate initial support values
    unordered_map<Edge, int> support;
    for (const auto& e : unassigned_edges) {
        support[e] = score_cycle_triangles(e, {}, out_neighbors, in_neighbors); // Use the appropriate score function based on triangle_type
    }

    int C_max = 0;
    for (const auto& kv : support) {
        C_max += kv.second;
    }
    C_max /= p;

    while (!unassigned_edges.empty()) {
        for (int i = 0; i < p; ++i) {
            if (C[i] >= C_max || C[i] >= epsilon * *min_element(C.begin(), C.end())) {
                continue;
            }

            Edge e;
            if (!D[i].empty()) {
                auto it = max_element(D[i].begin(), D[i].end(), [&](const Edge& a, const Edge& b) {
                    return score_cycle_triangles(a, partitions[i], out_neighbors, in_neighbors) < score_cycle_triangles(b, partitions[i], out_neighbors, in_neighbors);
                });
                e = *it;
                D[i].erase(it);
            } else {
                e = *unassigned_edges.begin();
            }

            partitions[i].insert(e);
            unassigned_edges.erase(e);
            C[i] += support[e];

            // Update D and support values
            for (int j = 0; j < p; ++j) {
                if (j != i) {
                    D[j].erase(e);
                }
            }
        }
    }

    return partitions;
}

std::vector<std::unordered_set<Edge>> stratified_balanced_partitioner(const Graph& G, int p) {
    std::vector<std::unordered_set<Edge>> partitions(p);
    std::unordered_map<Edge, TrussnessPair> trussness_map = compute_trussness_pairs(G);

    // Obtain k_c(e) and k_f(e) for each edge
    std::vector<Edge> edges = G.get_edges();
    std::vector<int> k_c(edges.size()), k_f(edges.size());
    for (size_t i = 0; i < edges.size(); ++i) {
        k_c[i] = trussness_map[edges[i]].k_c;
        k_f[i] = trussness_map[edges[i]].k_f;
    }

    int k_cmax = *std::max_element(k_c.begin(), k_c.end());
    int k_fmax = *std::max_element(k_f.begin(), k_f.end());
    int B_max = std::accumulate(k_f.begin(), k_f.end(), 0) / p;

    int i = 0;
    for (int k_j = 1; k_j <= k_cmax && i < p - 1; ++k_j) {
        int k_l = 1;
        while (k_l <= k_fmax) {
            std::unordered_set<Edge> H;
            for (size_t e = 0; e < edges.size(); ++e) {
                if (k_f[e] <= k_l && k_c[e] >= k_j) {
                    H.insert(edges[e]);
                }
            }
            if (std::accumulate(H.begin(), H.end(), 0, [&](int sum, const Edge& e) {
                return sum + trussness_map[e].k_f;
            }) <= B_max) {
                partitions[i] = H;
                ++i;
                break;
            }
            ++k_l;
        }
    }

    // Assign remaining edges to the last partition
    std::unordered_set<Edge> remaining_edges;
    for (const auto& edge : edges) {
        if (std::find(partitions.begin(), partitions.end(), edge) == partitions.end()) {
            remaining_edges.insert(edge);
        }
    }
    partitions[p - 1] = remaining_edges;

    return partitions;
}