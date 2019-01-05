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

tuple<int, int>
Index::getBestMatch(unordered_map<string, set<int>> &reference_index, unordered_map<string, set<int>> &sequence_index) {
    vector<int> index_hits;
    int k = reference_index.begin()->first.length();

    //sorting the sequence_index by kmer index (by appeareance in sequence)
    set<pair<string, set<int>>, Comparator> ordered_seq_index(sequence_index.begin(), sequence_index.end(), cmp);

    for (const auto &it : ordered_seq_index) {
        string kmer = it.first;
        set<int> positions = reference_index[kmer];

        for (const auto &pos : positions) {
            index_hits.push_back(pos);
        }
    }

    vector<int> group = groupByMargin(index_hits);
    vector<int> lis = LIS::find(group);

    //extend k - 1 positions after last hit
    return {lis[0], lis[lis.size() - 1] + k - 1};
}

vector<int> Index::groupByMargin(vector<int> positions) {
    vector<vector<int>> groups;

    for(int pos : positions){
        bool added = false;
        for(int i = 0, n = groups.size();i<n;i++){
            vector<int>& group = groups[i];
            if(abs(group.back() - pos) < INDEX_HIT_MARGIN){
                group.push_back(pos);
                added = true;
                break;
            }
        }

        if(!added){
            vector<int> group;
            group.push_back(pos);
            groups.push_back(group);
        }
    }

    //return the largest group
    sort(groups.begin(), groups.end(), [](const vector<int> &vector1, const vector<int> &vector2){
        return vector1.size() > vector2.size();
    });

    return groups[0];
}
