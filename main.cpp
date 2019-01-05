#include <iostream>
#include "Mutation.h"
#include "vector"
#include "DataFetcher.h"
#include "string"
#include <fstream>
#include <iostream>
#include "Functions.h"
#include "Index.h"
#include "Align.h"

using namespace std;

const int w = 11;
const int k = 11;

int main() {

    const string genomePath = "./data/lambda.fasta";
    const string mutatedPath = "./data/lambda_simulated_reads.fasta";

    DataFetcher fetcher (genomePath, mutatedPath);
    fetcher.loadData();

    unordered_map<string,set<int>> genomeMinimizers = Index::buildMinimizerIndex(fetcher.genome, w, k);
    unordered_map<string,set<int>> readMinimizers;
    vector<Mutation> mutations;

    for(string &mutatedGenomeRead : fetcher.mutatedGenomeReads) {
        readMinimizers = Index::buildMinimizerIndex(mutatedGenomeRead, w, k);
        tuple<int, int> mappedAreaBorders = Index::getBestMatch(genomeMinimizers, readMinimizers);
        string mappedArea = fetcher.genome.substr(get<0>(mappedAreaBorders), get<1>(mappedAreaBorders) - get<0>(mappedAreaBorders));
        vector<Mutation> areaMutations = BIOINFORMATIKA_ALIGN_H::align(mappedArea, mutatedGenomeRead, get<0>(mappedAreaBorders));
        mutations.insert(mutations.end(), areaMutations.begin(), areaMutations.end());
    }

    sort(mutations.begin(), mutations.end());

    for(auto &m : mutations) {
        cout << m;
    }

    return 0;
}
