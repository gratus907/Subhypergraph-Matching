#pragma once
#include "SubhypergraphMatching/DataHyperGraph.h"
#include "SubhypergraphMatching/PatternHyperGraph.h"

namespace GraphLib::SubHyperGraphMatching {
    struct SubHyperGraphMatchingOption{

    };

    class HyperCandidateSpace {
    private:
        SubHyperGraphMatchingOption opt;
        DataHyperGraph *data;
        PatternHyperGraph *query;
        std::vector<std::vector<int>> candidate_vertex_set_;
        std::vector<std::vector<int>> candidate_hyperedge_set_;

        std::vector<std::vector<int>> required_vertex_label_nbrs;
        std::vector<std::vector<int>> required_hyperedge_label_nbrs;

        std::vector<std::vector<int>> in_queue;

        std::vector<std::vector<int>> vertex_cs_index;
        std::vector<std::vector<int>> edge_cs_index;

        //@def vertex_cand_nbr_count[u][v][i]: size(CS(e|u,v)), where e is the i-th incident hyperedge of u
        std::vector<std::vector<std::vector<int>>> vertex_cand_nbr_count;

        //@def hyperedge_cand_nbr_count[e][f][i]: size(CS(u|e,f)), where u is the i-th contained vertex of e
        std::vector<std::vector<std::vector<int>>> hyperedge_cand_nbr_count;

        //@def vertex_incidence_count[u][v][i]: number of f in NbrCS(u, v), where f is the i-th incident hyperedge of v
        std::vector<std::vector<std::vector<int>>> vertex_incidence_count;

        //@def hyperedge_contained_count[e][f][i]: number of v in NbrCS(e, f), where v is the i-th contained vertex of f
        std::vector<std::vector<std::vector<int>>> hyperedge_contained_count;

        //@def vertex_label_nbr_count[u][v][l]: number of l-labeled hyperedges in NbrCS(u, v).
        std::vector<std::vector<std::vector<int>>> vertex_label_nbr_count;

        //@def hyperedge_label_nbr_count[e][f][i]: number of l-labeled vertices in NbrCS(e, f).
        std::vector<std::vector<std::vector<int>>> hyperedge_label_nbr_count;


    public:
        HyperCandidateSpace(DataHyperGraph *data_, PatternHyperGraph *query_, SubHyperGraphMatchingOption filter_option);
        ~HyperCandidateSpace();

        inline bool isVertexCandidate(int u, int v) const { return vertex_cs_index[u][v] >= 0; }
        inline bool isHyperedgeCandidate(int e, int f) const { return edge_cs_index[e][f] >= 0; }
        void BuildInitialHCS();
        void RefineHCS();
        void Initialize();
        void BuildHyperCandidateSpace();

        void PrintCSStatistics(int level = 0);

        bool HyperEdgeSafety(int e, int f);
        bool VertexSafety(int u, int v);

        void RemoveHyperedge(int e, int f);
        void RemoveVertex(int u, int v);
    };

    HyperCandidateSpace::HyperCandidateSpace(DataHyperGraph *data_, PatternHyperGraph *query_, SubHyperGraphMatchingOption filter_option) {
        opt = filter_option;
        data = data_;
        query = query_;
        Initialize();
        BuildHyperCandidateSpace();
    }

    HyperCandidateSpace::~HyperCandidateSpace() {
    }


    void HyperCandidateSpace::PrintCSStatistics(int level) {
        fprintf(log_to,"\e\[0;31m[CS Statistics]\n");
        int num_vertices = 0, num_edges = 0;
        for (int i = 0; i < candidate_vertex_set_.size(); i++) {
            num_vertices += candidate_vertex_set_[i].size();
            if (level >= 1) {
                int label = query->GetVertexLabel(i);
                fprintf(log_to, "[LOG] Vert %d(label=%d) : Matching %lu vertices\n",i,label,candidate_vertex_set_[i].size());
            }
        }
        fprintf(log_to,"Number of vertices: %d\n", num_vertices);
        for (int i = 0; i < candidate_hyperedge_set_.size(); i++) {
            num_edges += candidate_hyperedge_set_[i].size();
            if (level >= 1) {
                int label = query->GetHyperedgeLabel(i);
                fprintf(log_to, "[LOG] Edge %d(label=%d) : Matching %lu hyperedges\n",i,label,candidate_hyperedge_set_[i].size());
            }
        }
        fprintf(log_to,"Number of hyperedges: %d \e[0m\n", num_edges);
    };


    void HyperCandidateSpace::Initialize() {
        vertex_cs_index.resize(query->GetNumVertices(), std::vector<int>(data->GetNumVertices(), -1));
        candidate_vertex_set_.resize(query->GetNumVertices());

        vertex_cand_nbr_count.resize(query->GetNumVertices());
        vertex_incidence_count.resize(query->GetNumVertices());
        vertex_label_nbr_count.resize(query->GetNumVertices());

        for (int i = 0; i < query->GetNumVertices(); i++) {
            vertex_cand_nbr_count[i].resize(data->GetNumVertices(), std::vector<int>(query->GetDegree(i)));
            vertex_incidence_count[i].resize(data->GetNumVertices());
            for (int j = 0; j < data->GetNumVertices(); j++) {
                vertex_incidence_count[i][j].resize(data->GetDegree(j));
            }
            vertex_label_nbr_count[i].resize(data->GetNumVertices(), std::vector<int>(query->GetNumHyperedgeLabels()));
        }

        edge_cs_index.resize(query->GetNumHyperedges(), std::vector<int>(data->GetNumHyperedges(), -1));
        candidate_hyperedge_set_.resize(query->GetNumHyperedges());
        hyperedge_cand_nbr_count.resize(query->GetNumHyperedges());
        hyperedge_contained_count.resize(query->GetNumHyperedges());
        hyperedge_label_nbr_count.resize(query->GetNumHyperedges());
        in_queue.resize(query->GetNumHyperedges(), std::vector<int>(data->GetNumHyperedges(), 0));

        for (int i = 0; i < query->GetNumHyperedges(); i++) {
            hyperedge_cand_nbr_count[i].resize(data->GetNumHyperedges(), std::vector<int>(query->GetArity(i)));
            hyperedge_contained_count[i].resize(data->GetNumHyperedges());
            for (int j = 0; j < data->GetNumHyperedges(); j++) {
                hyperedge_contained_count[i][j].resize(data->GetArity(j));
            }
            hyperedge_label_nbr_count[i].resize(data->GetNumHyperedges(), std::vector<int>(query->GetNumVertexLabels()));
        }

        required_vertex_label_nbrs.resize(query->GetNumVertices(), std::vector<int>(query->GetNumHyperedgeLabels(), 0));
        for (int u = 0; u < query->GetNumVertices(); u++) {
            for (int e : query->GetIncidentHyperedges(u)) {
                required_vertex_label_nbrs[u][query->GetHyperedgeLabel(e)]++;
            }
        }

        required_hyperedge_label_nbrs.resize(query->GetNumHyperedges(), std::vector<int>(query->GetNumVertexLabels(), 0));
        for (int e = 0; e < query->GetNumHyperedges(); e++) {
            for (int u : query->GetHyperedge(e)) {
                required_hyperedge_label_nbrs[e][query->GetVertexLabel(u)]++;
            }
        }
    }

    void HyperCandidateSpace::BuildHyperCandidateSpace() {
        BuildInitialHCS();
        RefineHCS();
    }

    void HyperCandidateSpace::BuildInitialHCS() {
        for (int e = 0; e < query->GetNumHyperedges(); e++) {
            int query_hyperedge_label = query->GetHyperedgeLabel(e);
            for (int &f : data->GetHyperedgesByLabel(query_hyperedge_label)) {
                edge_cs_index[e][f] = candidate_hyperedge_set_[e].size();
                candidate_hyperedge_set_[e].push_back(f);
            }
        }
        for (int u = 0; u < query->GetNumVertices(); u++) {
            int query_vertex_label = query->GetVertexLabel(u);
            for (int &v : data->GetVerticesByLabel(query_vertex_label)) {
                for (int &e : query->GetIncidentHyperedges(u)) {
                    bool ok = false;
                    for (int &f : data->GetIncidentHyperedges(v)) {
                        ok |= isHyperedgeCandidate(e, f);
                    }
                    if (!ok) goto nxt_candidate;
                }
                vertex_cs_index[u][v] = candidate_vertex_set_[u].size();
                candidate_vertex_set_[u].push_back(v);
                nxt_candidate:;
            }
        }

        PrintCSStatistics(0);
        fflush(stderr);
//        for (int u = 0; u < query->GetNumVertices(); u++) {
//            fprintf(stderr, "CandVert[%d]: ",u);
//            for (int v : candidate_vertex_set_[u]) {
//                fprintf(stderr, "%d ",v);
//            }
//            fprintf(stderr, "\n");
//        }
//
//        for (int u = 0; u < query->GetNumHyperedges(); u++) {
//            fprintf(stderr, "CandEdge[%d]: ",u);
//            for (int v : candidate_hyperedge_set_[u]) {
//                fprintf(stderr, "%d ",v);
//            }
//            fprintf(stderr, "\n");
//        }
        for (int u = 0; u < query->GetNumVertices(); u++) {
            for (int v : candidate_vertex_set_[u]) {
                for (int i = 0; i < query->GetDegree(u); i++) {
                    int e = query->GetIncidentHyperedge(u, i);
                    for (int j = 0; j < data->GetDegree(v); j++) {
                        int f = data->GetIncidentHyperedge(v, j);
                        if (isHyperedgeCandidate(e, f)) {
                            vertex_cand_nbr_count[u][v][i]++;
                            vertex_incidence_count[u][v][j]++;
//                            fprintf(stderr, "Found (%d, %d[%d]) as incidence pair of (%d, %d)->%d\n", e, f, j, u, v,vertex_incidence_count[u][v][j]);
                            if (vertex_incidence_count[u][v][j] == 1) {
//                                fprintf(stderr, "Add 1 to vertex_label_nbr_count[%d][%d][%d]\n", u,v,data->GetHyperedgeLabel(f));
                                vertex_label_nbr_count[u][v][data->GetHyperedgeLabel(f)]++;
                            }
                        }
                    }
//                    fprintf(stderr, "Vertex_cand_nbr_count[%d][%d][%d]: %d\n", u,v,i,vertex_cand_nbr_count[u][v][i]);
                }
            }
        }

        for (int e = 0; e < query->GetNumHyperedges(); e++) {
            for (int f : candidate_hyperedge_set_[e]) {
                for (int i = 0; i < query->GetArity(e); i++) {
                    int u = query->GetContainedVertex(e, i);
                    for (int j = 0; j < data->GetArity(f); j++) {
                        int v = data->GetContainedVertex(f, j);
                        if (isVertexCandidate(u, v)) {
                            hyperedge_cand_nbr_count[e][f][i]++;
                            hyperedge_contained_count[e][f][j]++;
                            if (hyperedge_contained_count[e][f][j] == 1) {
                                hyperedge_label_nbr_count[e][f][data->GetVertexLabel(v)]++;
                            }
                        }
                    }
                }
            }
        }
    }


    bool HyperCandidateSpace::HyperEdgeSafety(int e, int f) {
        for (int i = 0; i < query->GetArity(e); i++) {
            if (hyperedge_cand_nbr_count[e][f][i] == 0) {
//                fprintf(stderr, "Remove Hyperedge (%d, %d) as there are no %dth nbr of %d\n",e,f,i,e);
                return false;
            }
        }
        for (int l = 0; l < query->GetNumVertexLabels(); l++) {
            if (hyperedge_label_nbr_count[e][f][l] < required_hyperedge_label_nbrs[e][l]) {
//                fprintf(stderr, "Remove Hyperedge (%d, %d) as there are only %d %d-nbrs while %d required\n",e,f,hyperedge_label_nbr_count[e][f][l],l,required_hyperedge_label_nbrs[e][l]);
                return false;
            }
        }
        return true;
    }

    void HyperCandidateSpace::RemoveHyperedge(int e, int f) {
//        for (int u : query->GetHyperedge(e)) {
//            for (int v : data->GetHyperedge(f)) {
//                fprintf(stderr, "%d", isVertexCandidate(u, v));
//            }
//            fprintf(stderr, "\n");
//        }
//        fprintf(stderr, "Remove HyperedgePair (%d, %d)\n",e,f);
        if (candidate_hyperedge_set_[e].size() == 1) {
            fprintf(stderr, "[ERROR] ??? HyperEdgeEmpty");
            exit(-1);
        }
        int idx = edge_cs_index[e][f];
        int last_hyperedge = candidate_hyperedge_set_[e].back();
        edge_cs_index[e][last_hyperedge] = idx;
        candidate_hyperedge_set_[e][idx] = last_hyperedge;
        edge_cs_index[e][f] = -1;
        candidate_hyperedge_set_[e].pop_back();
        for (int i = 0; i < query->GetArity(e); i++) {
            int u = query->GetContainedVertex(e, i);
            for (int j = 0; j < data->GetArity(f); j++) {
                int v = data->GetContainedVertex(f, j);
                if (isVertexCandidate(u, v)) {
                    int e_idx_for_u = query->GetInverseHyperedgeIndex(e, i);
                    int f_idx_for_v = data->GetInverseHyperedgeIndex(f, j);
                    vertex_cand_nbr_count[u][v][e_idx_for_u]--;
                    vertex_incidence_count[u][v][f_idx_for_v]--;
                    if (vertex_incidence_count[u][v][f_idx_for_v] == 0) {
                        vertex_label_nbr_count[u][v][data->GetHyperedgeLabel(f)]--;
                    }
                }
            }
        }
    }


    bool HyperCandidateSpace::VertexSafety(int u, int v) {
        for (int i = 0; i < query->GetDegree(u); i++) {
            if (vertex_cand_nbr_count[u][v][i] == 0) {
//                fprintf(stderr, "Remove Vertex (%d, %d) as there are no %dth nbr of %d\n",u,v,i,u);
                return false;
            }
        }
        for (int l = 0; l < query->GetNumHyperedgeLabels(); l++) {
            if (vertex_label_nbr_count[u][v][l] < required_vertex_label_nbrs[u][l]) {
//                fprintf(stderr, "Remove Vertex (%d, %d) as there are only %d %d-nbrs while %d required\n",u,v,vertex_label_nbr_count[u][v][l],l,required_vertex_label_nbrs[u][l]);
                return false;
            }
        }
        return true;
    }

    void HyperCandidateSpace::RemoveVertex(int u, int v) {
//        fprintf(stderr, "Remove VertexPair (%d, %d)\n",u,v);
        if (candidate_vertex_set_[u].size() == 1) {
            fprintf(stderr, "[ERROR] ??? VertexEmpty (%d, %d)\n", u, v);
            exit(-1);
        }
        int idx = vertex_cs_index[u][v];
        int last_vertex = candidate_vertex_set_[u].back();
        vertex_cs_index[u][last_vertex] = idx;
        candidate_vertex_set_[u][idx] = last_vertex;
        vertex_cs_index[u][v] = -1;
        candidate_vertex_set_[u].pop_back();
        for (int i = 0; i < query->GetDegree(u); i++) {
            int e = query->GetIncidentHyperedge(u, i);
            for (int j = 0; j < data->GetDegree(v); j++) {
                int f = data->GetIncidentHyperedge(v, j);
                if (isHyperedgeCandidate(e, f)) {
                    int u_idx_for_e = query->GetInverseVertexIndex(u, i);
                    int v_idx_for_f = data->GetInverseVertexIndex(v, j);
                    hyperedge_cand_nbr_count[e][f][u_idx_for_e]--;
                    hyperedge_contained_count[e][f][v_idx_for_f]--;
                    if (hyperedge_contained_count[e][f][v_idx_for_f] == 0) {
                        hyperedge_label_nbr_count[e][f][data->GetVertexLabel(v)]--;
                    }
                }
            }
        }
    }

    void HyperCandidateSpace::RefineHCS() {
        std::queue<std::pair<int, int>> refinement_queue;
        for (int e = 0; e < query->GetNumHyperedges(); e++) {
            for (int f : candidate_hyperedge_set_[e]) {
                refinement_queue.emplace(e, f);
                in_queue[e][f] = true;
            }
        }
        while (!refinement_queue.empty()) {
            auto [e, f] = refinement_queue.front();
            refinement_queue.pop();
            in_queue[e][f] = false;
            if (!HyperEdgeSafety(e, f)) {
                RemoveHyperedge(e, f);
                for (auto u : query->GetHyperedge(e)) {
                    for (auto v : data->GetHyperedge(f)) {
                        if (isVertexCandidate(u, v)) {
                            if (!VertexSafety(u, v)) {
                                RemoveVertex(u, v);
                                for (auto ec : query->GetIncidentHyperedges(u)) {
                                    for (auto fc : data->GetIncidentHyperedges(v)) {
                                        if (isHyperedgeCandidate(ec, fc) and (in_queue[ec][fc] == 0)) {
                                            refinement_queue.emplace(ec, fc);
                                            in_queue[ec][fc] = true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        PrintCSStatistics(0);
    }
}