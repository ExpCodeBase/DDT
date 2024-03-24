
#ifndef EXAMPLES_ANALYTICAL_APPS_DDT_DDT_H_
#define EXAMPLES_ANALYTICAL_APPS_DDT_DDT_H_

#include <grape/grape.h>

#include "ddt/ddt_context.h"

namespace grape {

template <typename FRAG_T>
class DDT : public ParallelAppBase<FRAG_T, DDTContext<FRAG_T>>,
            public ParallelEngine {
 public:
  INSTALL_PARALLEL_WORKER(DDT<FRAG_T>, DDTContext<FRAG_T>, FRAG_T)
  using vertex_t = typename fragment_t::vertex_t;

  static constexpr bool need_split_edges = true;

  void PEval(const fragment_t& frag, context_t& ctx,
             message_manager_t& messages) {
    using depth_type = typename context_t::depth_type;

    messages.InitChannels(thread_num(), 2 * 1023 * 64, 2 * 1024 * 64);

    ctx.current_depth = 1;

    vertex_t source;
    bool native_source = frag.GetInnerVertex(ctx.source_id, source);

    auto inner_vertices = frag.InnerVertices();

    // init double buffer which contains updated vertices using bitmap
    ctx.curr_inner_updated.Init(inner_vertices, GetThreadPool());
    ctx.next_inner_updated.Init(inner_vertices, GetThreadPool());

#ifdef PROFILING
    ctx.exec_time -= GetCurrentTime();
#endif

    auto& channel_0 = messages.Channels()[0];

    // run first round DDT, update unreached vertices
    if (native_source) {
      ctx.partial_result[source] = 0;
      auto oes = frag.GetOutgoingAdjList(source);
      for (auto& e : oes) {
        auto u = e.get_neighbor();
        if (ctx.partial_result[u] == std::numeric_limits<depth_type>::max()) {
          ctx.partial_result[u] = 1;
          if (frag.IsOuterVertex(u)) {
            channel_0.SyncStateOnOuterVertex<fragment_t>(frag, u);
          } else {
            ctx.curr_inner_updated.Insert(u);
          }
        }
      }
    }

#ifdef PROFILING
    ctx.exec_time += GetCurrentTime();
    ctx.postprocess_time -= GetCurrentTime();
#endif

    messages.ForceContinue();

#ifdef PROFILING
    ctx.postprocess_time += GetCurrentTime();
#endif
  }

  void IncEval(const fragment_t& frag, context_t& ctx,
               message_manager_t& messages) {
    using depth_type = typename context_t::depth_type;

    auto& channels = messages.Channels();

    depth_type next_depth = ctx.current_depth + 1;
    int thrd_num = thread_num();
    ctx.next_inner_updated.ParallelClear(GetThreadPool());

#ifdef PROFILING
    ctx.preprocess_time -= GetCurrentTime();
#endif

    messages.ParallelProcess<fragment_t, EmptyType>(
        thrd_num, frag, [&ctx](int tid, vertex_t v, EmptyType) {
          if (ctx.partial_result[v] == std::numeric_limits<depth_type>::max()) {
            ctx.partial_result[v] = ctx.current_depth;
            ctx.curr_inner_updated.Insert(v);
          }
        });

#ifdef PROFILING
    ctx.preprocess_time += GetCurrentTime();
    ctx.exec_time -= GetCurrentTime();
#endif

    double rate = 0;
    if (ctx.avg_degree > 10) {
      auto ivnum = frag.GetInnerVerticesNum();
      rate = static_cast<double>(
                 ctx.curr_inner_updated.ParallelCount(GetThreadPool())) /
             static_cast<double>(ivnum);
      if (rate > 0.1) {
        auto inner_vertices = frag.InnerVertices();
        auto outer_vertices = frag.OuterVertices();
        ForEach(outer_vertices, [next_depth, &frag, &ctx, &channels](
                                    int tid, vertex_t v) {
          if (ctx.partial_result[v] == std::numeric_limits<depth_type>::max()) {
            auto ies = frag.GetIncomingAdjList(v);
            for (auto& e : ies) {
              auto u = e.get_neighbor();
              if (ctx.curr_inner_updated.Exist(u)) {
                ctx.partial_result[v] = next_depth;
                channels[tid].SyncStateOnOuterVertex<fragment_t>(frag, v);
                break;
              }
            }
          }
        });
        ForEach(inner_vertices, [next_depth, &frag, &ctx](int tid, vertex_t v) {
          if (ctx.partial_result[v] == std::numeric_limits<depth_type>::max()) {
            auto oes = frag.GetOutgoingInnerVertexAdjList(v);
            for (auto& e : oes) {
              auto u = e.get_neighbor();
              if (ctx.curr_inner_updated.Exist(u)) {
                ctx.partial_result[v] = next_depth;
                ctx.next_inner_updated.Insert(v);
                break;
              }
            }
          }
        });
      } else {
        ForEach(ctx.curr_inner_updated, [next_depth, &frag, &ctx, &channels](
                                            int tid, vertex_t v) {
          auto oes = frag.GetOutgoingAdjList(v);
          for (auto& e : oes) {
            auto u = e.get_neighbor();
            if (ctx.partial_result[u] ==
                std::numeric_limits<depth_type>::max()) {
              ctx.partial_result[u] = next_depth;
              if (frag.IsOuterVertex(u)) {
                channels[tid].SyncStateOnOuterVertex<fragment_t>(frag, u);
              } else {
                ctx.next_inner_updated.Insert(u);
              }
            }
          }
        });
      }
    } else {
      ForEach(ctx.curr_inner_updated, [next_depth, &frag, &ctx, &channels](
                                          int tid, vertex_t v) {
        auto oes = frag.GetOutgoingAdjList(v);
        for (auto& e : oes) {
          auto u = e.get_neighbor();
          if (ctx.partial_result[u] == std::numeric_limits<depth_type>::max()) {
            ctx.partial_result[u] = next_depth;
            if (frag.IsOuterVertex(u)) {
              channels[tid].SyncStateOnOuterVertex<fragment_t>(frag, u);
            } else {
              ctx.next_inner_updated.Insert(u);
            }
          }
        }
      });
    }

#ifdef PROFILING
    ctx.exec_time += GetCurrentTime();
    ctx.postprocess_time -= GetCurrentTime();
#endif

    ctx.current_depth = next_depth;
    if (!ctx.next_inner_updated.Empty()) {
      messages.ForceContinue();
    }

    ctx.next_inner_updated.Swap(ctx.curr_inner_updated);

#ifdef PROFILING
    ctx.postprocess_time += GetCurrentTime();
#endif
  }
};

}  // namespace grape

#endif  // EXAMPLES_ANALYTICAL_APPS_DDT_DDT_H_
