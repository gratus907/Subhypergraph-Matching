#include "DataStructure/Graph.h"

namespace GraphLib {
/***
 * @brief Enumerate all small cycles
 */
    void Graph::EnumerateLocalTriangles() {
        local_triangles.resize(GetNumEdges());
        int num_triangles = 0;
        std::vector<int> nbr_edge_id(GetNumVertices(), -1);
        for (int i = 0; i < GetNumEdges(); i++) {
            auto &[u, v] = edge_list[i];
            for (int &fst_edge : GetAllIncidentEdges(u)) {
                nbr_edge_id[GetOppositePoint(fst_edge)] = fst_edge;
            }
            for (int &snd_edge : GetAllIncidentEdges(v)) {
                int opp = GetOppositePoint(snd_edge);
                if (nbr_edge_id[opp] != -1) {
                    local_triangles[i].emplace_back(std::tuple(opp, nbr_edge_id[opp], snd_edge));
                    num_triangles++;
                }
            }
            for (int &fst_edge : GetAllIncidentEdges(u)) {
                nbr_edge_id[GetOppositePoint(fst_edge)] = -1;
            }
        }
    }


    void Graph::EnumerateLocalFourCycles() {
        Timer timer; timer.Start();
        local_four_cycles.resize(GetNumEdges());
        double total_required = 0, done = 0;
        for (int i = 0; i < GetNumEdges(); i++) {
            auto &[u, v] = edge_list[i];
            total_required += GetDegree(u) * GetDegree(v);
        }
        if (total_required > 1e10) {
            local_four_cycles.resize(0);
            return;
        }

        long long num_four_cycles = 0;
        for (int i = 0; i < GetNumEdges(); i++) {
            auto &[u, v] = edge_list[i];
            if (GetDegree(u) < GetDegree(v)) {
                for (int &fourth_edge_opp : GetAllIncidentEdges(u)) {
                    int fourth_vertex = GetOppositePoint(fourth_edge_opp);
                    if (fourth_vertex == v) continue;
                    if (GetDegree(fourth_vertex) < GetDegree(v)) {
                        for (int &third_edge_opp : GetAllIncidentEdges(fourth_vertex)) {
                            int third_vertex = GetOppositePoint(third_edge_opp);
                            if (third_vertex == u) continue;
                            int snd_edge = GetEdgeIndex(v, third_vertex);
                            if (snd_edge != -1) {
                                int fourth_edge = fourth_edge_opp ^ 1;
                                int third_edge = third_edge_opp ^ 1;
                                int fst_diag = GetEdgeIndex(u, third_vertex);
                                int snd_diag = GetEdgeIndex(v, fourth_vertex);
                                local_four_cycles[i].emplace_back(FourMotif(
                                        {i, snd_edge, third_edge, fourth_edge},
                                        {fst_diag, snd_diag})
                                );
                                num_four_cycles++;
                            }
                        }
                    }
                    else {
                        for (int &snd_edge : GetAllIncidentEdges(v)) {
                            int third_vertex = GetOppositePoint(snd_edge);
                            if (third_vertex == u) continue;
                            int third_edge = GetEdgeIndex(third_vertex, fourth_vertex);
                            if (third_edge != -1) {
                                int fourth_edge = fourth_edge_opp ^ 1;
                                int fst_diag = GetEdgeIndex(u, third_vertex);
                                int snd_diag = GetEdgeIndex(v, fourth_vertex);
                                local_four_cycles[i].emplace_back(FourMotif(
                                        {i, snd_edge, third_edge, fourth_edge},
                                        {fst_diag, snd_diag})
                                );
                                num_four_cycles++;
                            }
                        }
                    }
                }
            }
            else {
                for (int &snd_edge : GetAllIncidentEdges(v)) {
                    int third_vertex = GetOppositePoint(snd_edge);
                    if (third_vertex == u) continue;
                    if (GetDegree(third_vertex) < GetDegree(u)) {
                        for (int &third_edge : GetAllIncidentEdges(third_vertex)) {
                            int fourth_vertex = GetOppositePoint(third_edge);
                            if (fourth_vertex == v) continue;
                            int fourth_edge = GetEdgeIndex(fourth_vertex, u);
                            if (fourth_edge != -1) {
                                int fst_diag = GetEdgeIndex(u, third_vertex);
                                int snd_diag = GetEdgeIndex(v, fourth_vertex);
                                local_four_cycles[i].emplace_back(FourMotif(
                                        {i, snd_edge, third_edge, fourth_edge},
                                        {fst_diag, snd_diag})
                                );
                                num_four_cycles++;
                            }
                        }
                    }
                    else {
                        for (int &fourth_edge_opp : GetAllIncidentEdges(u)) {
                            int fourth_vertex = GetOppositePoint(fourth_edge_opp);
                            if (fourth_vertex == v) continue;
                            int third_edge = GetEdgeIndex(third_vertex, fourth_vertex);
                            if (third_edge != -1) {
                                int fourth_edge = fourth_edge_opp ^ 1;
                                int fst_diag = GetEdgeIndex(u, third_vertex);
                                int snd_diag = GetEdgeIndex(v, fourth_vertex);
                                local_four_cycles[i].emplace_back(FourMotif(
                                        {i, snd_edge, third_edge, fourth_edge},
                                        {fst_diag, snd_diag})
                                );
                                num_four_cycles++;
                            }
                        }
                    }
                }
            }
            done += GetDegree(u) * GetDegree(v);
        }
        timer.Stop();
    }

    void Graph::ChibaNishizeki() {
        std::vector<int> vertices_by_degree(GetNumVertices(), 0);
        for (int i = 0; i < GetNumVertices(); i++) vertices_by_degree[i] = i;
        std::sort(vertices_by_degree.begin(), vertices_by_degree.end(),[this](int a, int b)->bool{
            return GetDegree(a) > GetDegree(b);
        });

    }
}