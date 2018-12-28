#include <iostream>
#include "Mutation.h"
#include "vector"
#include "DataFetcher.h"
#include "string"
#include <fstream>
#include <iostream>
#include "Functions.h"

using namespace std;

int main() {

    /*
    vector<Mutation> mutations;
    Mutation first (Substitution, 45, 'G');
    Mutation second (Insertion, 2, 'T');
    Mutation third (Deletion, 223, 'C');
    Mutation fourth (Substitution, 200, 'A');
    mutations.push_back(first);
    mutations.push_back(second);
    mutations.push_back(third);
    mutations.push_back(fourth);

    for(vector<Mutation>::iterator it = mutations.begin(); it != mutations.end(); ++it) {
        cout << *it;
    }
     */
    /*
    const string genomePath = "./data/ecoli.fasta";
    const string mutatedPath = "./data/ecoli_simulated_reads.fasta";

    DataFetcher fetcher (genomePath, mutatedPath);
    fetcher.loadData();

    cout << "Original genome:" << "\n";
    cout << fetcher.genome << "\n";

    cout << "Mutated genome:" << "\n";
    for (const auto &mutatedGenomeRead : fetcher.mutatedGenomeReads) {
        cout << "Next read:" << "\n";
        cout << mutatedGenomeRead << "\n";
    }

    */

    string seq = "abcdefghijkl";
    vector<string> kmers = Functions::getKMers(seq, 3);

    for(const auto &kmer : kmers){
        cout << kmer << "\n";
    }

    return 0;
}
