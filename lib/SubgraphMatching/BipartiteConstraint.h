#include "SubgraphMatching/DataGraph.h"
#include "SubgraphMatching/PatternGraph.h"
#include "DataStructure/Graph.h"
#include "Base/Base.h"
#include "Base/BasicAlgorithms.h"

namespace GraphLib::SubgraphMatching {
    struct BipartiteConstraint {
        int left_size = 0, right_size = 0;
        // adjacency list
        int **adj, *adj_size, **adj_index;
        bool *used;
        // store matching
        int *left, *right;
        // true if there is no matching, or matching might be broken
        bool matched = true;

        std::pair<int, int> last_removed_edge;

        BipartiteConstraint(int left_size_, int right_size_) {
            left_size = left_size_;
            right_size = right_size_;
            adj = new int *[left_size];
            adj_size = new int[left_size];
            adj_index = new int *[left_size];
            left = new int[left_size];
            right = new int[right_size];
            used = new bool[left_size];
            memset(left,  -1, sizeof(int) * left_size);
            memset(right, -1, sizeof(int) * right_size);
            for (int i = 0; i < left_size; i++) {
                adj[i] = new int[right_size];
                adj_size[i] = 0;
                adj_index[i] = new int[right_size];
            }
            matched = false;
        }

        void Reset(int left_size_, int right_size_) {
            memset(left, -1, sizeof(int) * left_size);
            memset(right, -1, sizeof(int) * right_size);
            memset(used, false, sizeof(bool) * left_size);
            for (int i = 0; i < left_size; i++) {
                adj_size[i] = 0;
                memset(adj[i], 0, sizeof(int) * right_size);
                memset(adj_index[i], 0, sizeof(int) * right_size);
            }
            left_size = left_size_;
            right_size = right_size_;
        }

        void AddEdge(int u, int v) {
            adj[u][adj_size[u]] = v;
            adj_index[u][v] = adj_size[u];
            adj_size[u]++;
        }

        void RemoveEdge(int u, int v, bool try_to_fix = false) {
            int removed_edge_idx = adj_index[u][v];
            if (adj_size[u] > 1) {
                adj_index[u][adj[u][adj_size[u] - 1]] = removed_edge_idx;
                std::swap(adj[u][adj_size[u] - 1], adj[u][removed_edge_idx]);
            }
            adj_size[u]--;
            if (try_to_fix) {
                if (left[u] == v) {
                    matched = false;
                    left[u] = right[v] = -1;
                    last_removed_edge = {u, v};
                    if (adj_size[u] == 0) return;
                    // Try remaining edges first, they are likely to be matched
                    for (int i = removed_edge_idx; i < adj_size[u]; i++) {
                        int v_new = adj[u][i];
                        if (right[v_new] == -1) {
                            right[v_new] = u;
                            left[u] = v_new;
                            matched = true;
                            return;
                        }
                    }
                    std::memset(used, false, sizeof(bool) * left_size);
                    matched = FindAugmentingPath(u);
                }
            }
        }

        void Recover() {
            int u, v; std::tie(u, v) = last_removed_edge;
            left[u] = v; right[v] = u;
            last_removed_edge = {-1, -1};
            matched = true;
        }

        bool FindAugmentingPath(int u) {
            if (used[u]) return false;
            used[u] = true;
            for (int i = 0; i < adj_size[u]; i++) {
                int v = adj[u][i];
                int c = right[v];
                if (c == -1 or FindAugmentingPath(c)) {
                    left[u] = v;
                    right[v] = u;
                    return true;
                }
            }
            return false;
        }

        int FindMaximumMatching() {
            memset(left, -1, sizeof(int) * left_size);
            memset(right, -1, sizeof(int) * right_size);
            int ans = 0;
            for (int u = 0; u < left_size; u++) {
                for (int i = 0; i < adj_size[u]; i++) {
                    int v = adj[u][i];
                    if (right[v] == -1) {
                        left[u] = v;
                        right[v] = u;
                        ans++;
                        break;
                    }
                }
            }
            for (int u = 0; u < left_size; u++) {
                if (left[u] == -1) {
                    std::memset(used, false, sizeof(bool) * left_size);
                    if (FindAugmentingPath(u)) {
                        ans++;
                    }
                }
            }
            matched = (ans == left_size);
            return ans;
        }
    };

}
