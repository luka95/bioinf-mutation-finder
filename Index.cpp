#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>
#include "Index.h"
#include "LIS.h"

typedef function<bool(pair<string, set<int>>, pair<string, set<int>>)> Comparator;
static Comparator cmp = [](pair<string, set<int>> const &pair1, pair<string, set<int>> const &pair2) {
    return *pair1.second.begin() < *pair2.second.begin();
};

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

/**
 * Finds the best match between reference and sequence given their indexes
 * @param reference_index
 * @param sequence_index
 * @return (int, int, int, int) reference begin, reference end, sequence begin, sequence end
 */
 //TODO return 4 values instead of 1
tuple<int, int>
Index::getBestMatch(unordered_map<string, set<int>> &reference_index, unordered_map<string, set<int>> &sequence_index) {
    vector<tuple<int, int>> index_hits;
    int k = reference_index.begin()->first.length();

    //sorting the sequence_index by kmer index (by appeareance in sequence)
    set<pair<string, set<int>>, Comparator> ordered_seq_index(sequence_index.begin(), sequence_index.end(), cmp);

    for (const auto &it : ordered_seq_index) {
        string kmer = it.first;
        set<int> positions = reference_index[kmer];
        set<int> sequence_positions = it.second;

        for (const auto &seq_pos : sequence_positions) {
            for (const auto &pos : positions) {
                index_hits.push_back({pos-seq_pos, pos});
            }
        }

    }

    //sort index_hits (by pos-seq_pos ascending)
    sort(index_hits.begin(), index_hits.end(), [](const tuple<int, int> &tuple1, const tuple<int, int> &tuple2) {
        return get<0>(tuple1) < get<0>(tuple2);
    });

    //clustering minimizer hits
    int b = 0;
    vector<vector<int>> groups;

    for (int e = 0, n = index_hits.size(); e < n; e++) {
        if (e == n - 1 || get<0>(index_hits[e + 1]) - get<0>(index_hits[e]) >= INDEX_HIT_MARGIN) {
            vector<tuple<int,int>> group;
            for (int i = b; i <= e; i++) {
                group.push_back(index_hits[i]);
            }
            b = e + 1;

            //sort group by seq_pos
            sort(group.begin(), group.end(), [](const tuple<int, int> &tuple1, const tuple<int, int> &tuple2) {
                return (get<1>(tuple1) - get<0>(tuple1)) < (get<1>(tuple2) - get<0>(tuple2));
            });

            vector<int> group_t;
            for(int i=0,m=group.size();i<m;i++){
                group_t.push_back(get<1>(group[i]));
            }

            groups.push_back(LIS::find(group_t));
        }
    }

    //find the largest group
    sort(groups.begin(), groups.end(), [](const vector<int> &vector1, const vector<int> &vector2) {
        if(vector1.size() == vector2.size()){
            return (vector1.back() - vector1.front()) > (vector2.back() - vector2.front());
        }else{
            return vector1.size() > vector2.size();
        }
    });

    int size = groups[0].size();
    //extend k - 1 positions after last hit
    return {groups[0][0], groups[0][groups[0].size() - 1] + k - 1};
}