#include "Mutation.h"
#include "DataLoader.h"
#include "string"
#include "vector"
#include <iostream>
#include <fstream>
#include <iostream>
#include <set>
#include <unordered_map>
#include <mapping/index/Index.h>
#include <algorithm>
#include <mapping/inverter/Inverter.h>
#include <alignment/Hirschberg.h>

using namespace std;

const int w = 5;
const int k = 15;

const string genome_path = "../data/ecoli.fasta";
const string mutated_path = "../data/ecoli_simulated_reads.fasta";
const string results_path = "../data/results_ecoli.csv";
const int COLLECTION_LIMIT = 4;
const int INSERTION_COLLECTION_LIMIT = 8;
const int MAX_READ_LENGTH = 10000;

tuple<short, char> getMax(map<char, short> &map);


int main() {
    DataLoader data_loader(genome_path, mutated_path);
    data_loader.loadData();

    unordered_map<string, set<tuple<int, int>>> genome_index = Index::index(data_loader.genome, w, k);

    // a counter for A,C,G,T,- for each position in the genome
    vector<map<char, short>> alignments;
    vector<map<char, short>> insertions;
    for (int i = 0, n = data_loader.genome.size(); i < n; i++) {
        alignments.push_back(map<char, short>());
        insertions.push_back(map<char, short>());
    }
    vector<Mutation> mutations;

    int total = data_loader.mutated_genome_reads.size();
    int processed = 0;

    int j;
    #pragma omp parallel for
    for (j = 0; j < total; j++) {

        string read = data_loader.mutated_genome_reads[j];

        #pragma omp critical(processed)
        {
            processed++;
            cout << "Proccesing " << processed << " of " << total << endl;
            cout <<"Genome index size : "<<genome_index.size()<<endl;
        };

        if (read.length() > MAX_READ_LENGTH) continue;

        unordered_map<string, set<tuple<int, int>>> read_index = Index::index(read, w, k);
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

        if (strand_xor == 1) {
            mapped_read = Inverter::inverse(mapped_read);
        }

        tuple<string, string> alignment = Hirschberg(mapped_genome, mapped_read);
        string reg_align = get<0>(alignment);
        string read_align = get<1>(alignment);

        #pragma omp critical(update)
        {
            int genome_pos = genome_start;
            int last_insertion = -1;
            for (int i = 0, n = reg_align.length(); i < n; i++) {
                char c = reg_align[i];

                if (c == '-') {
                    if (last_insertion == genome_pos) {
                        //taking only one symbol insertions
                        continue;
                    }
                    insertions[genome_pos][read_align[i]]++;
                    last_insertion = genome_pos;
                } else {
                    alignments[genome_pos][read_align[i]]++;
                    genome_pos++;
                }
            }
        };
    }

    //collect mutations - deletions and supstitutions
    for (int i = 0, n = data_loader.genome.length(); i < n; i++) {
        map<char, short> position_alignments = alignments[i];
        int total_hits = 0;
        for (auto &tup : position_alignments) {
            total_hits += tup.second;
        }

        if (total_hits < COLLECTION_LIMIT) {
            continue;
        }

        tuple<short, char> res = getMax(position_alignments);
        short occurences = get<0>(res);
        char c = get<1>(res);

        if (occurences <= total_hits / 2) {
            continue;
        }

        if (c == '-') {
            //deletion
            mutations.push_back(Mutation(MutationType::Deletion, i, c));
        } else if (c != data_loader.genome[i]) {
            //supstitution
            mutations.push_back(Mutation(MutationType::Substitution, i, c));
        }
    }

    //collect insertions
    for (int i = 0, n = data_loader.genome.length(); i < n; i++) {
        map<char, short> position_insertions = insertions[i];

        int total_hits = 0;
        for(auto &tup : position_insertions){
            total_hits+=tup.second;
        }

        if (total_hits < INSERTION_COLLECTION_LIMIT) {
            continue;
        }

        tuple<short, char> res = getMax(position_insertions);
        short occurences = get<0>(res);
        char c = get<1>(res);

        if (occurences <= total_hits / 2) {
            continue;

        }
        mutations.push_back(Mutation(MutationType::Insertion, i, c));
    }

    //output mutations
    ofstream file;
    file.open(results_path);
    sort(mutations.begin(), mutations.end());
    for (auto &mutation : mutations) {
        file << mutation;
    }

    return 0;
}

tuple<short, char> getMax(map<char, short> &mp) {
    short max = 0;
    char c;

    for (auto &tup : mp) {
        if (tup.second > max) {
            max = tup.second;
            c = tup.first;
        }
    }

    return {max, c};
}