#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>
#include <cmath>
#include "Index.h"
#include "LIS.h"
#include "Inverter.h"

typedef function<bool(pair<string, set<int>>, pair<string, set<int>>)> Comparator;
static Comparator cmp = [](pair<string, set<int>> const &pair1, pair<string, set<int>> const &pair2) {
    return *pair1.second.begin() < *pair2.second.begin();
};


unordered_map<string, set<tuple<int, int>>> Index::index(string &inputString, int w, int k){
    unordered_map<string, set<tuple<int, int>>> index;
    //length of window
    int l = w+k-1;

    for(int i=0,n=inputString.length()-l+1;i<n;i++){
        string window = inputString.substr(i,l);

        int minimizer_hash = getDoubleStrandMinimizer(window, k);

        //collect minimizers
        for (int j = 0; j < w; j++) {
            string kmer = window.substr(j, k);
            string inverse_kmer = Inverter::inverse(kmer);

            int kmer_hash = hash(kmer);
            int inv_hash = hash(inverse_kmer);

            if (kmer_hash<inv_hash && kmer_hash == minimizer_hash) {
                index[kmer].insert({i+j, 0});
            }else if(inv_hash < kmer_hash && inv_hash == minimizer_hash){
                index[inverse_kmer].insert({i+j, 1});
            }
        }
    }

    return index;
}

unsigned long long Index::hash(string& sequence){
    unsigned long long x = 0;
    int k = sequence.length();
    for(int i=0;i<k;i++){
        char c = sequence[i];
        switch (c) {
            case 'A':
                break;
            case 'C':
                x+=pow(4, k-i-1);
                break;
            case 'G':
                x+=2*pow(4,k-i-1);
                break;
            case 'T':
                x+=3*pow(4,k-i-1);
                break;
            default:
                break;
        }
    }
    unsigned long long m = static_cast<long>(pow(2, 2*k) - 1);
    x = (~x + (x<<21)) & m;
    x = x^x>>24;
    x = (x+(x<<3)+(x<<8))&m;
    x = x^x>>14;
    x = (x + (x<<2) + (x<<4)) & m;
    x = x^x>>28;
    x = (x + (x<<31))&m;
    return x;
}

int Index::getDoubleStrandMinimizer(string &window, int k) {
    int minimizer_hash = -1;

    for (int i = 0, n = window.length() - k + 1; i < n; i++) {
        string kmer = window.substr(i, k);
        string inverse_kmer = Inverter::inverse(kmer);

        int kmer_hash = hash(kmer);
        int inv_hash = hash(inverse_kmer);

        if (minimizer_hash == -1 || kmer_hash!=inv_hash) {
            minimizer_hash = min(kmer_hash, inv_hash);
        }
    }

    return minimizer_hash;
}



unordered_map<string, set<int>> Index::buildMinimizerIndex(string &inputString, int w, int k) {
    // length of window
    int l = w + k - 1;
    unordered_map<string, set<int>> index;

    for (int i = 0, n = inputString.length() - l + 1; i < n; i++) {
        string window = inputString.substr(i, l);

        tuple<string, int> tup = getMinimizer(window, k);
        string minimizer = get<0>(tup);
        int offset = i + get<1>(tup);

        index[minimizer].insert(offset);
    }

    unordered_map<string, set<int>> endMinimizers = getEndMinimizers(inputString, w, k);

    for (const auto &it : endMinimizers) {
        set<int> endIndices = it.second;
        for (const auto &idx : endIndices) {
            index[it.first].insert(idx);
        }
    }

    return index;
}

tuple<string, int> Index::getMinimizer(string &window, int k) {
    string minimizer;
    int offset = -1;

    for (int i = 0, n = window.length() - k + 1; i < n; i++) {
        string kmer = window.substr(i, k);
        if (minimizer.empty() || kmer.compare(minimizer) <= 0) {
            minimizer = kmer;
            offset = i;
        }
    }

    return {minimizer, offset};
}

unordered_map<string, set<int>> Index::getEndMinimizers(string &inputString, int w, int k) {
    unordered_map<string, set<int>> endMinimizersIndex;
    int n = inputString.length();

    for (int u = 1; u <= w - 1; u++) {
        int l = u + k - 1;
        string start = inputString.substr(0, l);
        string end = inputString.substr(n - l, l);

        tuple<string, int> tup = getMinimizer(start, k);
        endMinimizersIndex[get<0>(tup)].insert(get<1>(tup));

        tup = getMinimizer(end, k);
        endMinimizersIndex[get<0>(tup)].insert(n - l + get<1>(tup));
    }


    return endMinimizersIndex;
}


tuple<tuple<int, int, int , int>, int>
Index::getBestMatch(unordered_map<string, set<tuple<int,int>>> &reference_index, unordered_map<string, set<tuple<int, int>>> &sequence_index) {
    //strand_xor, diff_pos, reference_positions
    vector<tuple<int, int, int>> index_hits;
    int k = reference_index.begin()->first.length();

    for (const auto &it : sequence_index) {
        string kmer = it.first;
        //offset, strand
        set<tuple<int,int>> positions = reference_index[kmer];
        set<tuple<int,int>> sequence_positions = it.second;

        for (const auto &seq_pos : sequence_positions) {
            for (const auto &pos : positions) {
                int seq_position = get<0>(seq_pos);
                int seq_strand = get<1>(seq_pos);
                int ref_position = get<0>(pos);
                int ref_strand = get<1>(pos);

                if(seq_strand == ref_strand){
                    index_hits.push_back({0, ref_position-seq_position, ref_position});
                }else{
                    index_hits.push_back({1, ref_position+seq_position, ref_position});
                }
            }
        }

    }

    //sort index_hits (by strand_xor and then by diff_pos ascending)
    sort(index_hits.begin(), index_hits.end(), [](const tuple<int, int, int> &tuple1, const tuple<int, int, int> &tuple2) {
        if(get<0>(tuple1) == get<0>(tuple2)){
            return get<1>(tuple1) < get<1>(tuple2);
        } else{
            return get<0>(tuple1) < get<0>(tuple2);
        }
    });

    //clustering minimizer hits
    int b = 0;
    vector<tuple<vector<tuple<int,int>>, int>> groups;
    //one group is a vector of (seq_pos, ref_pos) + a strand_xor mark

    for (int e = 0, n = index_hits.size(); e < n; e++) {
        if (e == n - 1 || get<0>(index_hits[e+1])!=get<0>(index_hits[e])
                || get<1>(index_hits[e + 1]) - get<1>(index_hits[e]) >= INDEX_HIT_MARGIN) {

            vector<tuple<int, int>> group;
            int strand_xor = get<0>(index_hits[b]);

            for (int i = b; i <= e; i++) {
                tuple<int, int, int> index_hit = index_hits[i];
                if(strand_xor==0){
                    group.push_back({get<2>(index_hit) - get<1>(index_hit), get<2>(index_hit)});
                }else{
                    group.push_back({get<1>(index_hit) - get<2>(index_hit), get<2>(index_hit)});
                }
            }
            b = e + 1;

            //sort group by seq_pos
            sort(group.begin(), group.end(), [](const tuple<int, int> &tuple1, const tuple<int, int> &tuple2) {
                return get<0>(tuple1) < get<0>(tuple2);
            });

            vector<tuple<int, int>> lis_group = LIS::findBySecond(group, strand_xor);
            if(lis_group.size() < MINIMUM_MINIMIZER_MATCHES ||
               get<1>(lis_group.back()) - get<1>(lis_group.front()) < MINIMUM_MATCH_LENGTH){
                continue;
            }

            groups.push_back({lis_group, strand_xor});
        }
    }

    if(groups.size() == 0){
        return {{0,0,0,0}, -1};
    }

    //find the largest group
    sort(groups.begin(), groups.end(), [](const tuple<vector<tuple<int,int>>, int> &group1, const tuple<vector<tuple<int,int>>, int> &group2) {
        vector<tuple<int,int>> vector1 = get<0>(group1);
        vector<tuple<int,int>> vector2 = get<0>(group2);

        if(vector1.size() == vector2.size()){
            return (get<1>(vector1.back()) - get<1>(vector1.front())) > (get<1>(vector2.back()) - get<1>(vector2.front()));
        }else{
            return vector1.size() > vector2.size();
        }
    });


    tuple<vector<tuple<int, int>>, int> best_group = groups[0];
    vector<tuple<int,int>> group = get<0>(best_group);
    int strand_xor = get<1>(best_group);

    //reference_begin ,reference_end, seq_begin, seq_end   +  strand_xor
    //extend k - 1 positions after last hit
    if(strand_xor == 0) {
        return {{get<1>(group.front()), get<1>(group.back()) + k - 1, get<0>(group.front()),
                 get<0>(group.back()) + k - 1}, strand_xor};
    }else{
        //reverse sequence positions
        return {{get<1>(group.front()), get<1>(group.back()) + k - 1, get<0>(group.back()),
                 get<0>(group.front()) + k - 1}, strand_xor};
    }
}


