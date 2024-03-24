#ifndef EXAMPLES_ANALYTICAL_APPS_FLAGS_H_
#define EXAMPLES_ANALYTICAL_APPS_FLAGS_H_

#include <gflags/gflags_declare.h>

DECLARE_bool(directed);
DECLARE_string(application);
DECLARE_string(efile);
DECLARE_string(vfile);
DECLARE_string(out_prefix);
DECLARE_string(jobid);

DECLARE_int64(bfs_source);
DECLARE_int64(sssp_source);
DECLARE_double(pr_d);
DECLARE_int32(pr_mr);
DECLARE_int32(cdlp_mr);

DECLARE_int32(degree_threshold);

DECLARE_bool(segmented_partition);
DECLARE_bool(rebalance);
DECLARE_int32(rebalance_vertex_factor);

DECLARE_bool(serialize);
DECLARE_bool(deserialize);
DECLARE_string(serialization_prefix);

DECLARE_int32(app_concurrency);

DECLARE_string(lb);
#endif  // EXAMPLES_ANALYTICAL_APPS_FLAGS_H_
