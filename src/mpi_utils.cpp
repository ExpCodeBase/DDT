#include "mpi_utils.h"
#include <mpi.h>
#include <iostream>

vector<Edge> mpi_receive_edges(int rank, int size) {
    vector<Edge> received_edges;
    int received_count;
    MPI_Status status;
    for (int i = 0; i < size; ++i) {
        if (i != rank) {
            MPI_Probe(i, 0, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_BYTE, &received_count);
            received_edges.resize(received_count / sizeof(Edge));
            MPI_Recv(received_edges.data(), received_count, MPI_BYTE, i, 0, MPI_COMM_WORLD, &status);
        }
    }
    return received_edges;
}

void mpi_send_edges(const vector<Edge>& messages, int rank, int size) {
    for (int i = 0; i < size; ++i) {
        if (i != rank) {
            MPI_Send(messages.data(), messages.size() * sizeof(Edge), MPI_BYTE, i, 0, MPI_COMM_WORLD);
        }
    }
}

void collect_results(const vector<Edge>& local_edges, int rank, int size) {
    if (rank == 0) {
        vector<Edge> results;
        results.insert(results.end(), local_edges.begin(), local_edges.end());

        for (int i = 1; i < size; ++i) {
            int received_count;
            MPI_Status status;
            MPI_Probe(i, 0, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_BYTE, &received_count);
            vector<Edge> received_edges(received_count / sizeof(Edge));
            MPI_Recv(received_edges.data(), received_count, MPI_BYTE, i, 0, MPI_COMM_WORLD, &status);
            results.insert(results.end(), received_edges.begin(), received_edges.end());
        }

        // Print results
        for (const auto& edge : results) {
            cout << "Edge (" << edge.u << ", " << edge.v << "): Cycle Support = " << edge.cycle_support << ", Flow Support = " << edge.flow_support << endl;
        }
    } else {
        MPI_Send(local_edges.data(), local_edges.size() * sizeof(Edge), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
    }
}