#include <iostream>
#include <format>
#include "DataStructure/HyperGraph/HyperGraph.h"
#include "SubhypergraphMatching/DataHyperGraph.h"
#include "SubhypergraphMatching/PatternHyperGraph.h"
#include "SubgraphMatching/DataGraph.h"
#include "SubgraphMatching/PatternGraph.h"

#include "Base/Metrics.h"
#include "Base/Timer.h"
#include "DataStructure/Graph.h"
#include "SpecialSubgraphs/SmallCycle.h"
#include "SubgraphMatching/DataGraph.h"
#include "SubgraphMatching/PatternGraph.h"
#include "SubgraphMatching/CandidateSpace.h"
#include "SubgraphMatching/CandidateFilter.h"
#include "SubgraphCounting/Option.h"
#include "SubgraphMatching/Backtrack.h"

#include "SubhypergraphMatching/HyperCandidateSpace.h"
#include "SubhypergraphMatching/Preprocess.h"
using namespace std;
using namespace GraphLib;

const int MAX_NUM_VERTICES = 300;
Timer timer;
int32_t main(int argc, char *argv[]) {
    std::string dataset = "amazon-reviews";
    std::string query_name = "query_3_0", query;
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
                case 'd':
                    dataset = argv[i + 1];
                    break;
                case 'q':
                    query_name = argv[i + 1];
                    break;
            }
        }
    }
    query = "../dataset/hypergraphs/"+dataset+"/queries/"+query_name+".txt";

    GraphLib::SubHyperGraphMatching::SubHyperGraphMatchingOption opt_;
    GraphLib::SubHyperGraphMatching::PatternHyperGraph PG;
    PG.ReadPatternHyperGraph(query);
    PG.PrintStatistics("PatternGraph");
    if (PG.GetNumVertices() >= MAX_NUM_VERTICES) {
        return 5;
    }
    GraphLib::SubHyperGraphMatching::DataHyperGraph HG;
    HG.LoadDataGraph(dataset, "../dataset/hypergraphs/", PG);
    HG.PrintStatistics("Extracted DataGraph");

    timer.Start();
    GraphLib::SubHyperGraphMatching::HyperCandidateSpace HCS(&HG, &PG, opt_);
    timer.Stop();
    fprintf(stderr, "FilteringTime: %.02lf\n", timer.GetTime());

//    DataGraph D(HG.BipartiteRepresentation());
//    D.Preprocess();
//    D.EnumerateLocalTriangles();
//    D.EnumerateLocalFourCycles();
//    PatternGraph P(PG.BipartiteRepresentation());
//    GraphLib::SubgraphMatching::SubgraphMatchingOption opt;
//    GraphLib::SubgraphMatching::BacktrackEngine backtrack(&D, opt);
//    P.ProcessPattern(D);
//    timer.Start();
//    P.EnumerateLocalTriangles();
//    P.EnumerateLocalFourCycles();
//    backtrack.Match(&P);
//    cout << backtrack.num_embeddings << endl;
//    timer.Stop();
}