#include "flags.h"

#include <gflags/gflags.h>

#include <limits>

/* flags related to the job. */
DEFINE_string(application, "", "application name");
DEFINE_string(efile, "", "edge file");
DEFINE_string(vfile, "", "vertex file");
DEFINE_string(out_prefix, "", "output directory of results");
DEFINE_bool(directed, false, "input graph is directed or not.");

/* flags related to specific applications. */
DEFINE_int64(bfs_source, 0, "source vertex of bfs.");
DEFINE_int32(cdlp_mr, 10, "max rounds of cdlp.");
DEFINE_int64(sssp_source, 0, "source vertex of sssp.");
DEFINE_double(pr_d, 0.85, "damping_factor of pagerank");
DEFINE_int32(pr_mr, 10, "max rounds of pagerank");
DEFINE_int32(degree_threshold, std::numeric_limits<int>::max(),
             "Filtering threshold for some algorithms");

DEFINE_bool(segmented_partition, true,
            "whether to use segmented partitioning.");
DEFINE_bool(rebalance, true, "whether to rebalance graph after loading.");
DEFINE_int32(rebalance_vertex_factor, 0, "vertex factor of rebalancing.");

DEFINE_bool(serialize, false, "whether to serialize loaded graph.");
DEFINE_bool(deserialize, false, "whether to deserialize graph while loading.");
DEFINE_string(serialization_prefix, "",
              "where to load/store the serialization files");

DEFINE_int32(app_concurrency, -1, "concurrency of application");

DEFINE_string(lb, "cm",
              "Load balancing policy, these options can be used: auto, "
              "auto_static, none, cta, cm, wm, strict");
