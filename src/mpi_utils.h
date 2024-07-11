#ifndef MPI_UTILS_H
#define MPI_UTILS_H

#include <vector>
#include "graph.h"

using namespace std;

vector<Edge> mpi_receive_edges(int rank, int size);
void mpi_send_edges(const vector<Edge>& messages, int rank, int size);
void collect_results(const vector<Edge>& local_edges, int rank, int size);

#endif // MPI_UTILS_H