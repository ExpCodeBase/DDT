#ifndef GRAPH_H
#define GRAPH_H

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>

using namespace std;

struct Edge {
    int u, v;
    int cycle_support, flow_support;
};

void read_graph_data(const string& filename, int rank, int size, unordered_map<int, unordered_set<int>>& out_neighbors, unordered_map<int, unordered_set<int>>& in_neighbors, vector<Edge>& edges);
void compute_supports(unordered_map<int, unordered_set<int>>& out_neighbors, unordered_map<int, unordered_set<int>>& in_neighbors, vector<Edge>& edges);
void update_supports(vector<Edge>& edges, vector<Edge>& received_edges);
void peel_edges(vector<Edge>& edges, int k_c, int k_f);
void prepare_messages(vector<Edge>& edges, vector<Edge>& messages);

#endif // GRAPH_H