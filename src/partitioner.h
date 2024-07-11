#ifndef PARTITIONER_H
#define PARTITIONER_H

#include <vector>
#include <unordered_set>
#include "graph.h"

struct TrussnessPair {
    int k_c, k_f;
};

enum PartitionerType {
    TYPE_AWARE_BALANCED,
    STRATIFIED_BALANCED
};

std::vector<std::unordered_set<Edge>> type_aware_balanced_partitioner(const Graph& G, int p, double epsilon, const std::string& triangle_type);
std::vector<std::unordered_set<Edge>> stratified_balanced_partitioner(const Graph& G, int p);

#endif // PARTITIONER_H