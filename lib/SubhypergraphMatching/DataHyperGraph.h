#pragma once
#include "DataStructure/HyperGraph/HyperGraph.h"
#include "SubhypergraphMatching/PatternHyperGraph.h"

namespace GraphLib {
    namespace SubHyperGraphMatching {
        class DataHyperGraph : public GraphLib::HyperGraph {
        protected:
            std::vector<std::vector<int>> vertex_by_labels;
            std::vector<std::vector<int>> hyperedges_by_label;
        public:
            std::vector<int>& GetHyperedgesByLabel(const int l) {return hyperedges_by_label[l];}
            std::vector<int>& GetVerticesByLabel(const int l) {return vertex_by_labels[l];}
            void LoadDataGraph(std::string dataset, std::string path, PatternHyperGraph &P);
        };

        void DataHyperGraph::LoadDataGraph(std::string dataset, std::string path, PatternHyperGraph &P) {
            std::string hyperedge_file = path + "/" + dataset + "/hyperedges-" + dataset + ".txt";
            std::string vertex_label_file = path + "/" + dataset + "/node-labels-" + dataset + ".txt";
            std::cerr << "Read " << fileSize(hyperedge_file.c_str()) << " bytes from " << hyperedge_file << endl;
            std::ifstream fin(vertex_label_file);
            std::string line;
            std::vector<int> tmp_vertex_label;
            while (getline(fin, line)) {
                int l = stoi(line)-1;
                tmp_vertex_label.push_back(P.GetMappedVertexLabel(l));
            }
            num_vertex = tmp_vertex_label.size();
            std::vector<int> vertex_used(num_vertex, -1);

            fin = std::ifstream(hyperedge_file);
            while (getline(fin, line)) {
                auto edge = parse(line, ",");
                std::vector<int> current_hyperedge;
                std::vector<int> current_signature;
                for (auto &elem : edge) {
                    current_hyperedge.push_back(stoi(elem)-1);
                }
                std::sort(current_hyperedge.begin(), current_hyperedge.end());
                current_hyperedge.erase(std::unique(current_hyperedge.begin(), current_hyperedge.end()), current_hyperedge.end());
                if (current_hyperedge.size() == 1) { continue; }
                for (auto &elem : current_hyperedge) {
                    current_signature.push_back(tmp_vertex_label[elem]);
                }
                std::sort(current_signature.begin(), current_signature.end());
                int l = P.GetMappedHyperedgeLabel(current_signature);
                if (l == -1) continue;
                for (auto &elem : current_hyperedge) {
                    vertex_used[elem] = 0;
                }
                total_arity += current_hyperedge.size();
                hyperedges.push_back(current_hyperedge);
            }
            num_vertex = 0;
            for (int i = 0; i < vertex_used.size(); i++) {
                if (vertex_used[i] >= 0) {
                    vertex_label.push_back(tmp_vertex_label[i]);
                    vertex_used[i] = num_vertex++;
                }
            }
            for (int i = 0; i < hyperedges.size(); i++) {
                for (int j = 0; j < hyperedges[i].size(); j++) {
                    hyperedges[i][j] = vertex_used[hyperedges[i][j]];
                }
            }
            std::sort(hyperedges.begin(), hyperedges.end());
            hyperedges.erase(std::unique(hyperedges.begin(), hyperedges.end()), hyperedges.end());
            num_edge = hyperedges.size();
            for (auto &E : hyperedges) {
                hyperedge_signatures.push_back(std::vector<int>());
                for (int x : E) {
                    hyperedge_signatures[hyperedge_signatures.size()-1].push_back(vertex_label[x]);
                }
                std::sort(hyperedge_signatures[hyperedge_signatures.size()-1].begin(), hyperedge_signatures[hyperedge_signatures.size()-1].end());
                hyperedge_label.push_back(P.GetMappedHyperedgeLabel(hyperedge_signatures.back()));
            }

            num_vertex_labels = P.GetNumVertexLabels();
            num_hyperedge_labels = P.GetNumHyperedgeLabels();

            hyperedges_by_label.resize(GetNumHyperedgeLabels());
            for (int i = 0; i < GetNumHyperedges(); i++) {
                hyperedges_by_label[GetHyperedgeLabel(i)].push_back(i);
            }
            vertex_by_labels.resize(GetNumVertexLabels());
            for (int i = 0; i < GetNumVertices(); i++) {
                vertex_by_labels[GetVertexLabel(i)].push_back(i);
            }
            BuildIncidenceList();
            BuildNeighborIndex();
        }
    }
}