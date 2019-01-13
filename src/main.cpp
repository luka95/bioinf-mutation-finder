#include <iostream>
#include "Mutation.h"
#include "vector"
#include "string"
#include "DataLoader.h"
#include <fstream>
#include <iostream>
#include <set>
#include <unordered_map>
#include <mapping/index/Index.h>
#include <algorithm>
#include <mapping/inverter/Inverter.h>
#include "Hirschberg.cpp"

using namespace std;

const int w = 5;
const int k = 15;

const string genome_path = "../data/lambda.fasta";
const string mutated_path = "../data/lambda_simulated_reads.fasta";


int main() {

    DataLoader data_loader(genome_path, mutated_path);
    data_loader.loadData();

    unordered_map<string, set<tuple<int, int>>> genome_index = Index::index(data_loader.genome, w, k);
    unordered_map<string, set<tuple<int, int>>> read_index;

    //TODO Open MP for parallelization here
    for (string &read : data_loader.mutated_genome_reads) {
        read_index = Index::index(read, w, k);
        tuple<tuple<int, int, int, int>, int> mapping = Index::getBestMatch(genome_index, read_index);

        int strand_xor = get<1>(mapping);
        if (strand_xor == -1) continue;

        tuple<int, int, int, int> positions = get<0>(mapping);
        int genome_start = get<0>(positions);
        int genome_end = get<1>(positions);
        int read_start = get<2>(positions);
        int read_end = get<3>(positions);

        string mapped_read = read.substr(read_start, read_end - read_start + 1);
        string mapped_genome = data_loader.genome.substr(genome_start, genome_end - genome_start + 1);

        if(strand_xor == 1){
            mapped_read = Inverter::inverse(mapped_read);
        }

        zw alignment = Hirschberg(mapped_genome, mapped_read);
    }

    return 0;
}
