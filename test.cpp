#include<omp.h>
#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include <iostream>
#include <DataLoader.h>
#include <set>
#include <unordered_map>
#include <mapping/index/Index.h>
#include <map>

using namespace std;
const string genome_path = "../data/ecoli.fasta";
const string mutated_path = "../data/ecoli_simulated_reads.fasta";

class BaseCounter{
private:
    short cnt[5] = {0,0,0,0,0};
public:
    void increaseCount(char base){
        switch(base) {
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
};


int main(){
    DataLoader data_loader(genome_path, mutated_path);
    data_loader.loadData();
//    BaseCounter* alignments = new BaseCounter[data_loader.genome.size()];
//    BaseCounter* insertions = new BaseCounter[data_loader.genome.size()];

    unordered_map<string, set<tuple<int, int>>> genome_index = Index::index(data_loader.genome, 5, 15);

    int i;
    int n = data_loader.mutated_genome_reads.size();
#pragma omp parallel for
    for(i=0; i<n;i++){
        string read = data_loader.mutated_genome_reads[i];
        auto index = Index::index(read, 5, 15);
        auto mapping = Index::getBestMatch(genome_index, index);
        cout<<i<<endl;
    }
//
//    delete[] alignments;
//    delete[] insertions;
}