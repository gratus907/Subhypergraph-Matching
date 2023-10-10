#pragma once
#include <boost/dynamic_bitset.hpp>
#include "SubgraphMatching/BipartiteConstraint.h"
#include "SubgraphMatching/CandidateSpace.h"


const int INVALID = -2;
namespace GraphLib {
    static std::string space(int idx) {
        std::string s;
        for (int i = 0; i < 2*idx; i++) {s.push_back(' ');}
        return s;
    }

    void print(boost::dynamic_bitset<> const & b, int m = -1, int idx = -1) {
        if (m == -1) m = b.size();
        if (idx >= 0) std::cout << space(idx);
        for (int i = 0; i < m; ++i) {
            std::cout << b[i];
        }
        std::cout << std::endl;
    }

namespace SubgraphMatching {
    class BacktrackEngine {
    public:
        unsigned long long num_embeddings = 0, traversed_nodes = 0, pruned_nodes = 0;
    private:
        BipartiteConstraint *GlobalBipartiteConstraint;
        DataGraph *data_;
        PatternGraph *query_;
        CandidateSpace *CS;
        SubgraphMatchingOption opt_;
        int *seen, *isolated_vertex_candidates;
        std::vector<std::vector<std::vector<int>>> local_candidates;
        std::vector<int> which_local_candidate;
        int root = -1;
        std::vector<int> M;
        unsigned long long cnt = 0, conflicts = 0, dead_end = 0, bp_failure = 0;


        void print_cands(std::vector<int> &v, int who) {
            for (auto &x : v) {
                fprintf(stderr, "%d(%d) ", x, CS->GetCandidate(who, x));
            }
            fprintf(stderr, "\n");
        }

        int ChooseExtendableVertex(int idx) {
            int u = -1, cmu_size = 1e9;
            for (int i = 0; i < query_->GetNumVertices(); i++) {
                if (which_local_candidate[i] == -1) continue;
                if (M[i] != -1) continue;
                int local_candidate_size = local_candidates[i][which_local_candidate[i]].size();
                if (local_candidate_size == 0)
                    return INVALID;
                if (local_candidate_size == 1) {
                    return i;
                }
                if (which_local_candidate[i] + 1 == query_->GetDegree(i)) {
                    continue;
                }
                if (local_candidate_size < cmu_size) {
                    cmu_size = local_candidate_size;
                    u = i;
                }
            }
            return u;
        }

        bool PropagateExtendableVertex(int u, int v_idx, int print_idx=0) {
            auto &q_nbrs = query_->GetNeighbors(u);
            for (int i = 0; i < q_nbrs.size(); i++) {
                int q_nbr = q_nbrs[i];
                auto &cand_nbr = CS->GetCandidateNeighbors(u, v_idx, q_nbr);
//                fprintf(stderr, "!PropagateExtendableVertex %d->%d to sz = %lu : ", u, q_nbr, cand_nbr.size());

                if (which_local_candidate[q_nbr] == -1) {
                    which_local_candidate[q_nbr] = 0;
                    std::copy(cand_nbr.begin(), cand_nbr.end(),
                              std::back_inserter(local_candidates[q_nbr][which_local_candidate[q_nbr]]));
                }
                else {
                    auto &prev_cand = local_candidates[q_nbr][which_local_candidate[q_nbr]];
                    std::set_intersection(
                            cand_nbr.begin(), cand_nbr.end(),
                            prev_cand.begin(), prev_cand.end(),
                            std::back_inserter(local_candidates[q_nbr][which_local_candidate[q_nbr]+1]));
                    which_local_candidate[q_nbr]++;
                }
                if (local_candidates[q_nbr][which_local_candidate[q_nbr]].size() == 0) {
                    for (int j = i; j >= 0; j--) {
                        int processed_q_nbr = q_nbrs[j];
                        if (which_local_candidate[processed_q_nbr] >= 0) {
                            local_candidates[processed_q_nbr][which_local_candidate[processed_q_nbr]].clear();
                            which_local_candidate[processed_q_nbr]--;
                        }
                    }
                    return false;
                }
//                fprintf(stderr, "Set local_candidates[%d][%d] to size %lu vector : ", q_nbr, which_local_candidate[q_nbr],
//                        local_candidates[q_nbr][which_local_candidate[q_nbr]].size());
//                print_cands(local_candidates[q_nbr][which_local_candidate[q_nbr]], q_nbr);
            }
            return true;
        }

        void ReleaseNeighbors(int u) {
            for (int q_nbr : query_->GetNeighbors(u)) {
                if (which_local_candidate[q_nbr] >= 0) {
                    local_candidates[q_nbr][which_local_candidate[q_nbr]].clear();
                    which_local_candidate[q_nbr]--;
                }
            }
        }

        UnionFind isolated_vertex_groups;
        std::vector <std::vector<int>> isolated_bipartite_graph;
        std::vector <int> bp_cand_idx, isolated_vertices, isolated_candidates;
        long long bipartite_count_dp[1<<15][2];
        std::vector<int> distinct_candidates;
        long long CountMaximumMatchings(int l, int r) {
            if (l + 1 == r) {
                long long a, b, common;
                a = b = common = 0;
                int u1 = isolated_vertices[l];
                auto &u1_cands = local_candidates[u1][which_local_candidate[u1]];
                int u2 = isolated_vertices[r];
                auto &u2_cands = local_candidates[u2][which_local_candidate[u2]];
                for (int &uc : u1_cands) {
                    int v = CS->GetCandidate(u1, uc);
                    if (seen[v] != -1) continue;
                    bp_cand_idx[v] = 1; a++;
                }
                for (int &uc : u2_cands) {
                    int v = CS->GetCandidate(u2, uc);
                    if (seen[v] != -1) continue;
                    if (bp_cand_idx[v] == 1) common++;
                    b++;
                }
                for (int &uc : u1_cands){
                    int v = CS->GetCandidate(u1, uc);
                    if (seen[v] != -1) continue;
                    bp_cand_idx[v] = -1;
                }
                return (a - common) * b + common * (b - 1);
            }
            else {
                cnt++;
                int num_distinct_candidates = 0;
                for (int i = l; i <= r; i++) {
                    int u = isolated_vertices[i];
                    auto &u_cands = local_candidates[u][which_local_candidate[u]];
                    for (int &uc : u_cands) {
                        int v = CS->GetCandidate(u, uc);
                        if (seen[v] != -1) continue;
                        if (bp_cand_idx[v] == -1) {
                            distinct_candidates.push_back(v);
                            bp_cand_idx[v] = ++num_distinct_candidates;
                        }
                        isolated_bipartite_graph[bp_cand_idx[v]].push_back(i-l);
                    }
                }
                int max_bitmask = (1 << (r - l + 1));

                for (int i = 0; i < max_bitmask; i++) {
                    bipartite_count_dp[i][0] = 0;
                }
                int cur = 1, bef = 0;
                bipartite_count_dp[0][0] = bipartite_count_dp[0][1] = 1;
                for (int i = 1; i <= num_distinct_candidates; i++) {
                    for (int b = 1; b < max_bitmask; b++) {
                        bipartite_count_dp[b][cur] = bipartite_count_dp[b][bef];
                    }
                    for (int left_idx : isolated_bipartite_graph[i]) {
                        for (int b = 1; b < max_bitmask; b++) {
                            if (b & (1 << left_idx)) {
                                int b_prev = b - (1<<left_idx);
                                bipartite_count_dp[b][cur] += bipartite_count_dp[b_prev][bef];
                            }
                        }
                    }
                    cur = 1 - cur; bef = 1 - bef;
                }
                for (int i = l; i <= r; i++) {
                    int u = isolated_vertices[i];
                    auto &u_cands = local_candidates[u][which_local_candidate[u]];
                    for (int &uc : u_cands) {
                        int v = CS->GetCandidate(u, uc);
                        if (seen[v] != -1) continue;
                        if (bp_cand_idx[v] != -1) {
                            isolated_bipartite_graph[bp_cand_idx[v]].clear();
                            bp_cand_idx[v] = -1;
                        }
                    }
                }
                distinct_candidates.clear();
                return bipartite_count_dp[max_bitmask-1][num_distinct_candidates%2];
            }
        }
        void RevertIsolatedCandidates() {
            for (auto &v : isolated_candidates) {
                isolated_vertex_candidates[v] = -1;
            }
            isolated_candidates.clear();
        }
        unsigned long long MatchIsolatedVertices() {
            long long num_answers = 1;
            isolated_vertices.clear();
            isolated_vertex_groups.init();
            std::vector <long long> temp_candidates(query_->GetNumVertices(), 0);
            for (int i = 0; i < query_->GetNumVertices(); i++) {
                if (M[i] == -1 and which_local_candidate[i] + 1 == query_->GetDegree(i)) {
                    isolated_vertices.push_back(i);
                    auto &candidates = local_candidates[i][which_local_candidate[i]];
                    for (int idx = 0; idx < candidates.size(); idx++) {
                        int v = CS->GetCandidate(i, candidates[idx]);
                        if (seen[v] != -1) {
                            continue;
                        }
                        if (isolated_vertex_candidates[v] == -1) {
                            isolated_vertex_candidates[v] = i;
                            isolated_candidates.push_back(v);
                        }
                        else {
                            int grouped_query_vertex = isolated_vertex_candidates[v];
                            isolated_vertex_groups.unite(grouped_query_vertex, i);
                        }
                        temp_candidates[i]++;
                    }
                    if (temp_candidates[i] == 0) {
                        bp_failure++;
                        return 0;
                    }
                }
            }
            std::sort(isolated_vertices.begin(), isolated_vertices.end(),
                      [this](int a, int b)->bool {
                    return isolated_vertex_groups.find(a) < isolated_vertex_groups.find(b);
            });
            int st = 0, ed = 0;
            int k = isolated_vertices.size();
            while (st < k and ed < k) {
                if ((ed + 1 < k) and
                    (isolated_vertex_groups.find(isolated_vertices[ed+1]) == isolated_vertex_groups.find(isolated_vertices[ed]))) {
                    ed++;
                }
                else {
                    if (st == ed)
                        num_answers *= temp_candidates[isolated_vertices[st]];
                    else
                        num_answers *= CountMaximumMatchings(st, ed);
                    st = ed = ed+1;
                }
            }
            return num_answers;
        }


//        bool CheckForEmbedding() {
//            for (int u = 0; u < query_->GetNumVertices(); u++) {
//                if (M[u] == -1){
//                    fprintf(stderr, "??? Incomplete Embedding was detected: ");
//                    for (int i = 0; i < query_->GetNumVertices(); i++) {
//                        fprintf(stderr, "%d(%d) ", M[i],CS->GetCandidate(i, M[i]));
//                    }
//                    fprintf(stderr, "\n#NUM_MATCHED_NBR/NUM_NBR: ");
//                    for (int i = 0; i < query_->GetNumVertices(); i++) {
//                        fprintf(stderr, "%d/%d ", which_local_candidate[i]+1,query_->GetDegree(i));
//                    }
//                    return false;
//                }
//                int v = CS->GetCandidate(u, M[u]);
//                seen[v] = 0;
//            }
//
//            for (int u = 0; u < query_->GetNumVertices(); u++) {
//                int v = CS->GetCandidate(u, M[u]);
//                if (seen[v] == 1) {
//                    fprintf(stderr, "??? Non-Injective Embedding was detected: ");
//                    for (int i = 0; i < query_->GetNumVertices(); i++) {
//                        fprintf(stderr, "%d(%d) ", M[i],CS->GetCandidate(i, M[i]));
//                    }
//                    return false;
//                }
//                seen[v] = 1;
//            }
//
//            for (int u = 0; u < query_->GetNumVertices(); u++) {
//                int v = CS->GetCandidate(u, M[u]);
//                seen[v] = 0;
//            }
//
//            for (int u = 0; u < query_->GetNumVertices(); u++) {
//                for (int uc : query_->GetNeighbors(u)) {
//                    int v = CS->GetCandidate(u, M[u]);
//                    int vc = CS->GetCandidate(uc, M[uc]);
//                    if (data_->GetEdgeIndex(v, vc) == -1) {
//                        fprintf(stderr, "??? Edge Missing Embedding was detected: ");
//                        for (int i = 0; i < query_->GetNumVertices(); i++) {
//                            fprintf(stderr, "%d(%d) ", M[i],CS->GetCandidate(i, M[i]));
//                        }
//                        fprintf(stderr, "%d-%d should be an edge, but it is not\n",v,vc);
//                        return false;
//                    }
//                }
//            }
//            return true;
//        }

        bool FindEmbeddings(int idx) {
            traversed_nodes++;
            if (traversed_nodes % 5'000'000 == 0) {
                fprintf(stderr, "Traversed nodes: %llu Embeddings: %llu\n", traversed_nodes, num_embeddings);
                fflush(stderr);
            }
            int u = ChooseExtendableVertex(idx);
            if (u == -1) {
                unsigned long long num_found = MatchIsolatedVertices();
                num_embeddings += num_found;
                RevertIsolatedCandidates();
                return (num_found > 0);
            }
            bool found = false;
            std::vector<int>& loc_cands = local_candidates[u][which_local_candidate[u]];
            for (int v_idx : loc_cands) {
                int v = CS->GetCandidate(u, v_idx);
                if (seen[v] != -1) {
                    conflicts++;
                    continue;
                }
                M[u] = v_idx;
                seen[v] = u;
                bool extendable = PropagateExtendableVertex(u, M[u], idx);
                if (extendable) {
                    found |= FindEmbeddings(idx+1);
                    ReleaseNeighbors(u);
                }
                else dead_end++;
                M[u] = -1;
                seen[v] = -1;
            }
            return found;
        }

    public:
        BacktrackEngine(DataGraph *data, SubgraphMatchingOption opt) {
            data_ = data;
            opt_ = opt;
            CS = new CandidateSpace(data, opt);
            seen = new int[data->GetNumVertices()];
            isolated_vertex_candidates = new int[data->GetNumVertices()];
            GlobalBipartiteConstraint = new BipartiteConstraint(opt.MAX_QUERY_VERTEX, data->GetNumVertices());
            memset(seen, -1, sizeof(int) * data->GetNumVertices());
            memset(isolated_vertex_candidates, -1, sizeof(int) * data->GetNumVertices());
            which_local_candidate.resize(opt.MAX_QUERY_VERTEX);
            local_candidates.resize(opt.MAX_QUERY_VERTEX);
            for (int i = 0; i < opt.MAX_QUERY_VERTEX; i++) {
                local_candidates[i].resize(opt.MAX_QUERY_VERTEX);
            }
            isolated_bipartite_graph.resize(data_->GetNumVertices());
            for (int i = 0; i < data_->GetNumVertices(); i++) {
                isolated_bipartite_graph.reserve(64);
            }
            bp_cand_idx.resize(data_->GetNumVertices(), -1);
            isolated_vertices.reserve(opt_.MAX_QUERY_VERTEX);
            isolated_candidates.reserve(data_->GetNumVertices());
            distinct_candidates.reserve(64);
        };
        ~BacktrackEngine(){
            delete[] seen;
        };


        void ResetOptions() {
            if (!data_->FourCycleEnumerated()) {
                CS->opt.structure_filter = std::min(opt_.structure_filter, TRIANGLE_SAFETY);
            }
            CS->opt.use_cs_index = false;
        }
        void Match(PatternGraph *query) {
            ResetOptions();
            query_ = query;
            num_embeddings = traversed_nodes = pruned_nodes = bp_failure = dead_end = conflicts = 0;
            memset(seen, -1, data_->GetNumVertices());
            memset(isolated_vertex_candidates, -1, data_->GetNumVertices());
            std::fill(which_local_candidate.begin(), which_local_candidate.end(), -1);
            isolated_vertex_groups = UnionFind(query_->GetNumVertices());
            M.resize(query_->GetNumVertices(), -1);
            CS->BuildCS(query_);
            if (opt_.max_num_matches == 0) {
                return;
            }
            std::vector <int> num_cands(query_->GetNumVertices());
            for (int i = 0; i < query_->GetNumVertices(); i++) {
                num_cands[i] = CS->GetCandidateSetSize(i);
            }
            root = std::min_element(num_cands.begin(), num_cands.end()) - num_cands.begin();
            for (int i = 0; i < CS->GetCandidateSetSize(root); i++) {
                int v = CS->GetCandidate(root, i);
                M[root] = i;
                seen[v] = root;
                if (PropagateExtendableVertex(root, i, 0))
                    FindEmbeddings(1);
                ReleaseNeighbors(root);
                seen[v] = -1;
                M[root] = -1;
            }
            printf("Total BP Count : %llu\n",cnt);
            printf("Dead end nodes : %llu\n",dead_end);
            printf("Conflicts : %llu\n",conflicts);
            printf("BP-Failure : %llu\n",bp_failure);
        };

        unsigned long long GetNumEmbeddings() {
            return num_embeddings;
        }

        unsigned long long GetNumRecursiveCalls() {
            return traversed_nodes;
        }
    };
}
}