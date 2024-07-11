#include "peeler.h"
#include <unordered_map>
#include <mpi.h>

void compute_supports(std::vector<Graph>& partitions, int rank) {
    Graph& F = partitions[rank];
    std::unordered_map<int, std::unordered_set<int>> out_neighbors;
    std::unordered_map<int, std::unordered_set<int>> in_neighbors;

    // Initialize supports to zero
    for (auto& edge : F.edges) {
        edge.cycle_support = 0;
        edge.flow_support = 0;
    }

    // Build adjacency lists
    for (const auto& edge : F.edges) {
        out_neighbors[edge.u].insert(edge.v);
        in_neighbors[edge.v].insert(edge.u);
    }

    // Compute supports
    for (auto& edge : F.edges) {
        int u = edge.u;
        int v = edge.v;

        // Cycle support
        for (int w : in_neighbors[u]) {
            if (out_neighbors[v].find(w) != out_neighbors[v].end()) {
                edge.cycle_support++;
            }
        }

        // Flow support
        for (int w : out_neighbors[u]) {
            if (in_neighbors[v].find(w) != in_neighbors[v].end()) {
                edge.flow_support++;
            }
        }
        for (int w : in_neighbors[v]) {
            if (out_neighbors[u].find(w) != out_neighbors[u].end()) {
                edge.flow_support++;
            }
        }
        for (int w : out_neighbors[v]) {
            if (in_neighbors[u].find(w) != in_neighbors[u].end()) {
                edge.flow_support++;
            }
        }
    }

    // Synchronize supports across all workers
    int edges_count = F.edges.size();
    MPI_Allreduce(MPI_IN_PLACE, &edges_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, F.edges.data(), edges_count * sizeof(Edge), MPI_BYTE, MPI_SUM, MPI_COMM_WORLD);
}

std::tuple<int, int> distributedMaxTrussNumbers(std::vector<Graph>& partitions, int rank) {
    int local_max_k_c = 0, local_max_k_f = 0;
    for (const auto& edge : G.edges) {
        for (const auto& pair : edge.pairs) {
            local_max_k_c = std::max(local_max_k_c, pair.k_c);
            local_max_k_f = std::max(local_max_k_f, pair.k_f);
        }
    }

    int global_max_k_c, global_max_k_f;
    MPI_Allreduce(&local_max_k_c, &global_max_k_c, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&local_max_k_f, &global_max_k_f, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    return std::make_tuple(global_max_k_c, global_max_k_f);
}

void update_cycle_support(Edge& edge, int new_support) {
    edge.cycle_support = new_support;
}

void update_flow_support(Edge& edge, int new_support) {
    edge.flow_support = new_support;
}

void UpSupc(Edge& e_star, std::unordered_map<int, std::unordered_set<int>>& out_neighbors, std::unordered_map<int, std::unordered_set<int>>& in_neighbors, std::unordered_map<int, int>& LinkSta) {
    int u = e_star.u;
    int v = e_star.v;

    // (1) Zero the cycle support for e_star itself
    update_cycle_support(e_star, 0);

    // (2) Traverse edges spanned between {u, v} and T_i^io(e_star)
    for (int h : out_neighbors[u]) {
        if (out_neighbors[v].find(h) != out_neighbors[v].end()) {
            // Decrement cycle support for edges incident to h
            for (auto& edge : edges_incident_to(h)) {
                update_cycle_support(edge, edge.cycle_support - 1);
            }
            // Disqualify u and v in the local indexes of h
            disqualify_in_local_index(h, u, v);
        }
    }

    // Update LinkSta because of the removed edge e_star
    LinkSta[u * num_vertices + v] = 0;
    LinkSta[v * num_vertices + u] = 0;
}

void UpSupf(Edge& e_star, std::unordered_map<int, std::unordered_set<int>>& out_neighbors, std::unordered_map<int, std::unordered_set<int>>& in_neighbors, std::unordered_map<int, int>& LinkSta) {
    int u = e_star.u;
    int v = e_star.v;

    // (1) Zero the flow support for e_star itself
    update_flow_support(e_star, 0);

    if (LinkSta[u * num_vertices + v] == 1 || LinkSta[v * num_vertices + u] == 1) {
        // (2) If u and v are uni-linked
        for (int h : in_neighbors[u]) {
            if (in_neighbors[v].find(h) != in_neighbors[v].end()) {
                // Decrement flow support for edges incident to h
                for (auto& edge : edges_incident_to(h)) {
                    update_flow_support(edge, edge.flow_support - 1);
                }
                // Disqualify u and v in the local indexes of h
                disqualify_in_local_index(h, u, v);
            }
        }
    } else if (LinkSta[u * num_vertices + v] == 2 && LinkSta[v * num_vertices + u] == 2) {
        // (3) If u and v are bi-linked
        for (int h : out_neighbors[u]) {
            if (out_neighbors[v].find(h) != out_neighbors[v].end()) {
                // Decrement flow support for edges incident to h
                for (auto& edge : edges_incident_to(h)) {
                    update_flow_support(edge, edge.flow_support - 1);
                }
                // Disqualify u and v in the local indexes of h
                disqualify_in_local_index(h, u, v);
            }
        }
    }

    // Update LinkSta because of the removed edge e_star
    LinkSta[u * num_vertices + v] = 0;
    LinkSta[v * num_vertices + u] = 0;
}

void perform_flow_decomposition(std::unordered_set<Edge>& edges, std::unordered_map<int, std::unordered_set<int>>& out_neighbors, std::unordered_map<int, std::unordered_set<int>>& in_neighbors, std::unordered_map<int, int>& LinkSta) {
    // Initialize k_f to 0
    int k_f = 0;

    // Continue until all edges are peeled
    while (!edges.empty()) {
        // Peel edges with unqualified flow support
        for (auto it = edges.begin(); it != edges.end(); ) {
            if (it->flow_support < k_f) {
                UpSupf(*it, out_neighbors, in_neighbors, LinkSta);
                it = edges.erase(it);
            } else {
                ++it;
            }
        }

        // Increase k_f
        k_f++;
    }
}

void stratified_local_peel_processing(std::unordered_set<Edge>& F_i, const std::unordered_map<Edge, int>& cycle_supports, const std::unordered_map<Edge, int>& flow_supports, int k_c, int k_f, std::unordered_map<int, std::unordered_set<int>>& out_neighbors, std::unordered_map<int, std::unordered_set<int>>& in_neighbors, std::unordered_map<int, int>& LinkSta) {
    // Implementation of the stratified local-peel processing logic
    // Ensure that all edges which form a flow triangle within the interval on each fragment are included
    std::unordered_set<Edge> F_i_alpha;

    for (const auto& edge : F_i) {
        // Check if the edge forms a flow triangle within the interval
        if (is_flow_triangle_within_interval(edge, F_i)) {
            F_i_alpha.insert(edge);
        }
    }

    // Update cycle support for edges in F_i_alpha
    for (const auto& edge : F_i_alpha) {
        if (cycle_supports.find(edge) != cycle_supports.end()) {
            // Update cycle support for the edge
            update_cycle_support(edge, cycle_supports.at(edge));
        }
    }

    // Update flow support for edges in F_i_alpha
    for (const auto& edge : F_i_alpha) {
        if (flow_supports.find(edge) != flow_supports.end()) {
            // Update flow support for the edge
            update_flow_support(edge, flow_supports.at(edge));
        }
    }

    // Perform flow decomposition by lines 1-6 of Alg.~\ref{algo:disbatpeel}
    perform_flow_decomposition(F_i_alpha, out_neighbors, in_neighbors, LinkSta);

    // Peel edges based on the updated cycle and flow supports
    for (auto it = F_i.begin(); it != F_i.end(); ) {
        if (it->cycle_support < k_c || it->flow_support < k_f) {
            it = F_i.erase(it);
        } else {
            ++it;
        }
    }
}