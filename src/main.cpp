#include <mpi.h>
#include "graph.h"
#include "partitioner.h"
#include "peeler.h"
#include "mpi_utils.h"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Load graph data
    unordered_map<int, unordered_set<int>> out_neighbors;
    unordered_map<int, unordered_set<int>> in_neighbors;
    vector<Edge> local_edges;
    read_graph_data("graph.e", rank, size, out_neighbors, in_neighbors, local_edges);

    // Select partitioner type
    PartitionerType partitioner_type = STRATIFIED_BALANCED; // or STRATIFIED_BALANCED
    std::vector<std::unordered_set<Edge>> partitions;

    // Perform partitioning
    if (partitioner_type == TYPE_AWARE_BALANCED) {
        partitions = type_aware_balanced_partitioner(G, size, 1.1, "c");
    } else if (partitioner_type == STRATIFIED_BALANCED) {
        partitions = stratified_balanced_partitioner(G, size);
    }

    // Compute supports
    compute_supports(partitions, rank);

    // Compute maximal truss numbers for cycle and flow
    int k_c_max, k_f_max;
    std::make_tuple(k_c_max, k_f_max)  = distributedMaxTrussNumbers(partitions, rank);

    int superstep = 1;
    bool terminate = false;

    // Perform local-peel processing on each partition
    if (rank < size) {
        std::unordered_set<Edge> F_i = partitions[rank];
        if (partitioner_type == STRATIFIED_BALANCED) {
            stratified_local_peel_processing(F_i, cycle_supports, flow_supports, k_c_max, k_f_max, out_neighbors, in_neighbors, LinkSta);
        }
    }


    MPI_Finalize();
    return 0;
}