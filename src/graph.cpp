#include "graph.h"
#include <fstream>
#include <iostream>

void read_graph_data(const string& filename, int rank, int size, unordered_map<int, unordered_set<int>>& out_neighbors, unordered_map<int, unordered_set<int>>& in_neighbors, vector<Edge>& edges) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Failed to open file: " << filename << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int u, v;
    int edge_count = 0;
    while (file >> u >> v) {
        if (edge_count % size == rank) {
            out_neighbors[u].insert(v);
            in_neighbors[v].insert(u);
            edges.push_back({u, v, 0, 0});
        }
        edge_count++;
    }

    file.close();
}

void compute_supports(unordered_map<int, unordered_set<int>>& out_neighbors, unordered_map<int, unordered_set<int>>& in_neighbors, vector<Edge>& edges) {
    // Compute initial cycle and flow supports for edges
    for (auto& edge : edges) {
        int u = edge.u, v = edge.v;
        unordered_set<int> cycle_set, flow_set;

        for (int w : out_neighbors[v]) {
            if (out_neighbors[u].find(w) != out_neighbors[u].end()) {
                cycle_set.insert(w);
            }
        }

        for (int w : in_neighbors[u]) {
            if (in_neighbors[v].find(w) != in_neighbors[v].end()) {
                flow_set.insert(w);
            }
        }

        edge.cycle_support = cycle_set.size();
        edge.flow_support = flow_set.size();
    }
}

void update_supports(vector<Edge>& edges, vector<Edge>& received_edges) {
    // Update cycle and flow supports based on received messages
    for (auto& received_edge : received_edges) {
        for (auto& edge : edges) {
            if (edge.u == received_edge.u && edge.v == received_edge.v) {
                edge.cycle_support = received_edge.cycle_support;
                edge.flow_support = received_edge.flow_support;
            }
        }
    }
}

void peel_edges(vector<Edge>& edges, int k_c, int k_f) {
    // Peel edges that do not meet the support criteria
    edges.erase(remove_if(edges.begin(), edges.end(), [k_c, k_f](const Edge& edge) {
        return edge.cycle_support < k_c || edge.flow_support < k_f;
    }), edges.end());
}

void prepare_messages(vector<Edge>& edges, vector<Edge>& messages) {
    // Prepare messages to be sent in the next superstep
    messages = edges;
}