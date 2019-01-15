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

const char bases[] = {'A', 'C', 'G', 'T', '-'};
const int BASE_COUNT = 4;

const int w = 5;
const int k = 15;

const string genome_path = "../data/ecoli.fasta";
const string mutated_path = "../data/ecoli_simulated_reads.fasta";
const string results_path = "../data/results_ecoli.csv";
const int COLLECTION_LIMIT = 4;
const int INSERTION_COLLECTION_LIMIT = 8;
const int MAX_READ_LENGTH = 10000;

class BaseCounter {
private:
    short cnt[5] = {0, 0, 0, 0, 0};
public:
    void increaseCount(char base) {
        switch (base) {
            case 'A':
                cnt[0] += 1;
                break;
            case 'C':
                cnt[1] += 1;
                break;
            case 'G':
                cnt[2] += 1;
                break;
            case 'T':
                cnt[3] += 1;
                break;
            case '-':
                cnt[4] += 1;
                break;
            default:
                break;
        }
    }

    int getTotalHits() {
        int total_hits = 0;
        for (short s : cnt) {
            total_hits += s;
        }
        return total_hits;
    }

    tuple<short, char> getMax() {
        short max = 0;
        int idx;

        for (int i = 0; i < (BASE_COUNT + 1); i++) {
            short s = cnt[i];
            if (s > max) {
                max = s;
                idx = i;
            }
        }
        return {max, bases[idx]};
    }
};


int main() {
    DataLoader data_loader(genome_path, mutated_path);
    data_loader.loadData();

    vector<BaseCounter> alignments;
    vector<BaseCounter> insertions;
    for (int i = 0, n = data_loader.genome.size(); i < n; i++) {
        alignments.push_back(BaseCounter());
        insertions.push_back(BaseCounter());
    }

    vector<Mutation> mutations;

    cout<<"Building genome index"<<endl;
    unordered_map<string, set<tuple<int, int>>> genome_index = Index::index(data_loader.genome, w, k);

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
            //cout<<"Genome index size : "<<data_loader.genome.size()<<endl;
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
                    insertions[genome_pos].increaseCount(read_align[i]);
                    last_insertion = genome_pos;
                } else {
                    alignments[genome_pos].increaseCount(read_align[i]);
                    genome_pos++;
                }
            }
        };
    }

    //collect mutations - deletions and supstitutions
    for (int i = 0, n = data_loader.genome.length(); i < n; i++) {
        BaseCounter position_alignments = alignments[i];
        int total_hits = position_alignments.getTotalHits();
        if (total_hits < COLLECTION_LIMIT) {
            continue;
        }

        tuple<short, char> res = position_alignments.getMax();
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
        BaseCounter position_insertions = insertions[i];

        int total_hits = position_insertions.getTotalHits();
        if (total_hits < INSERTION_COLLECTION_LIMIT) {
            continue;
        }

        tuple<short, char> res = position_insertions.getMax();
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