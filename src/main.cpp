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


char findMostFrequentValue(vector<char> &values);

int main() {

    DataLoader data_loader(genome_path, mutated_path);
    data_loader.loadData();

    unordered_map<string, set<tuple<int, int>>> genome_index = Index::index(data_loader.genome, w, k);
    unordered_map<string, set<tuple<int, int>>> read_index;
    map<int, vector<char>> alignments;
    vector<Mutation> mutations;

    int processed = 1;
    int total = data_loader.mutated_genome_reads.size();

    //TODO Open MP for parallelization here
    for (string &read : data_loader.mutated_genome_reads) {

        cout<<"Processing read "<<processed<<" of "<<total<<endl;
        processed++;

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

        for (int i = 0, n = reg_align.length(); i < n; i++) {
            char c = reg_align[i];
            if (c == '-') {
                //insertion - mark = lower case base
                alignments[genome_start + i].push_back(static_cast<char &&>(tolower(read_align[i])));
            } else {
                // '-' in read_align is for deletion, a letter is for supstitution or match
                alignments[genome_start + i].push_back(read_align[i]);
            }
        }
    }


    //collect mutations
    for (int i = 0, n = data_loader.genome.length(); i < n; i++) {
        vector<char> position_alignments = alignments[i];
        if(position_alignments.size() == 0) continue;
        char c = findMostFrequentValue(position_alignments);
        if (islower(c)) {
            //insertion
            mutations.push_back(Mutation(MutationType::Insertion, i, c));
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
    for(auto &mutation : mutations){
        file<<mutation;
    }

    return 0;
}

char findMostFrequentValue(vector<char> &values) {

    int max = -1;
    char mf_value;
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

    return mf_value;
}
