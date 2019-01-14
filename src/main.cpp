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

const string genome_path = "../data/lambda.fasta";
const string mutated_path = "../data/lambda_simulated_reads.fasta";
const string results_path = "../data/results.csv";


tuple<int, char> findMostFrequentValue(vector<char> &values);

int main() {
    DataLoader data_loader(genome_path, mutated_path);
    data_loader.loadData();

    unordered_map<string, set<tuple<int, int>>> genome_index = Index::index(data_loader.genome, w, k);
    unordered_map<string, set<tuple<int, int>>> read_index;
    map<int, vector<char>> alignments;
    vector<Mutation> mutations;

    int total = data_loader.mutated_genome_reads.size();
    int processed = 0;
    //TODO Open MP for parallelization here
    int j;

    #pragma omp parallel
    {
        #pragma omp for
        for (j = 0; j < total; j++) {

            string read = data_loader.mutated_genome_reads[j];
            //cout << "Processing read " <<j+1<< " of " << total << endl;

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

            if (strand_xor == 1) {
                mapped_read = Inverter::inverse(mapped_read);
            }

            zw alignment = Hirschberg(mapped_genome, mapped_read);
            string reg_align = alignment.z;
            string read_align = alignment.w;

            #pragma omp critical
            {
                int genome_pos = genome_start;
                for (int i = 0, n = reg_align.length(); i < n; i++) {
                    char c = reg_align[i];

                    if (c == '-') {
                        alignments[genome_pos].push_back(static_cast<char>(tolower(read_align[i])));
                    } else {
                        alignments[genome_pos].push_back(read_align[i]);
                        genome_pos++;
                    }
                };

                processed++;
                cout << "Proccesed " << processed << " of " << total << endl;
            }
        }
    }

    //collect mutations
    int limit = 3;
    for (int i = 0, n = data_loader.genome.length(); i < n; i++) {
        vector<char> position_alignments = alignments[i];
        if (position_alignments.size() <= limit) {
            continue;
        }

        tuple<int, char> res = findMostFrequentValue(position_alignments);
        int occurences = get<0>(res);
        char c = get<1>(res);

        if (occurences <= position_alignments.size() / 2) {
            continue;
        }

        if (islower(c)) {
            //insertion
            mutations.push_back(Mutation(MutationType::Insertion, i, static_cast<char>(toupper(c))));
        } else if (c == '-') {
            //deletion
            mutations.push_back(Mutation(MutationType::Deletion, i, c));
        } else if (c != data_loader.genome[i]) {
            //supstitution
            mutations.push_back(Mutation(MutationType::Substitution, i, c));
        }
    }

    //output mutations
    ofstream file;
    file.open(results_path);
    for (auto &mutation : mutations) {
        file << mutation;
    }

    return 0;
}

tuple<int, char> findMostFrequentValue(vector<char> &values) {
    if (values.empty()) return {0, '\0'};

    int max = 1;
    char mf_value = values[0];
    int cnt = 1;

    sort(values.begin(), values.end());
    for (int i = 1, n = values.size(); i < n; i++) {
        if (values[i - 1] == values[i]) {
            cnt++;
        } else {
            if (cnt > max) {
                max = cnt;
                mf_value = values[i - 1];
                cnt = 1;
            }
        }
    }

    return {max, mf_value};
}

