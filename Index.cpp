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

    //sort index_hits (ascending)
    sort(index_hits.begin(), index_hits.end(), [](const tuple<int, int> &tuple1, const tuple<int, int> &tuple2) {
        if (get<0>(tuple1) == get<0>(tuple2)) {
            return get<1>(tuple1) < get<1>(tuple2);
        } else {
            return get<0>(tuple1) < get<0>(tuple2);
        }
    });

    //clustering minimizer hits
    int b = 0;
    vector<vector<int>> groups;

    for (int e = 0, n = index_hits.size(); e < n; e++) {
        if (e == n - 1 || get<0>(index_hits[e + 1]) - get<0>(index_hits[e]) >= INDEX_HIT_MARGIN) {
            vector<int> group;
            for (int i = b; i <= e; i++) {
                group.push_back(get<1>(index_hits[i]));
            }
            b = e + 1;
            groups.push_back(LIS::find(group));
        }
    }

    //find the largest group
    sort(groups.begin(), groups.end(), [](const vector<int> &vector1, const vector<int> &vector2) {
        return (vector1.back() - vector1.front()) > (vector2.back() - vector2.front());
    });

    //extend k - 1 positions after last hit
    return {groups[0][0], groups[0][groups[0].size() - 1] + k - 1};
}

vector<int> Index::groupByMargin(vector<int> positions) {
    vector<vector<int>> groups;

    for (int pos : positions) {
        bool added = false;
        for (int i = 0, n = groups.size(); i < n; i++) {
            vector<int> &group = groups[i];
            if (abs(group.back() - pos) < INDEX_HIT_MARGIN) {
                group.push_back(pos);
                added = true;
                break;
            }
        }

        if (!added) {
            vector<int> group;
            group.push_back(pos);
            groups.push_back(group);
        }
    }

    //merge groups
    while (true) {
        vector<vector<int>> merged_groups;
        set<int> merged;
        int n = groups.size();

        for (int i = 0; i < n; i++) {
            vector<int> group1 = groups[i];
            if (merged.find(i) != merged.end()) continue;

            for (int j = 0; j < n; j++) {
                vector<int> group2 = groups[j];
                if (i == j || merged.find(j) != merged.end()) continue;

                if (abs(group1.back() - group2.front()) < INDEX_HIT_MARGIN) {
                    vector<int> merged_group;

                    merged_group.insert(merged_group.begin(), group1.begin(), group1.end());
                    merged_group.insert(merged_group.end(), group2.begin(), group2.end());
                    merged_groups.push_back(merged_group);

                    merged.insert(i);
                    merged.insert(j);
                    break;
                } else if (abs(group2.back() - group1.front()) < INDEX_HIT_MARGIN) {
                    vector<int> merged_group;

                    merged_group.insert(merged_group.begin(), group2.begin(), group2.end());
                    merged_group.insert(merged_group.end(), group1.begin(), group1.end());
                    merged_groups.push_back(merged_group);

                    merged.insert(i);
                    merged.insert(j);
                    break;
                }
            }
        }

        if (merged.size() == 0) break;

        //add all not merged
        for (int i = 0; i < n; i++) {
            if (merged.find(i) != merged.end()) continue;
            merged_groups.push_back(groups[i]);
        }
        groups = merged_groups;
    }

    int a = 3;
    //return the largest group
    sort(groups.begin(), groups.end(), [](const vector<int> &vector1, const vector<int> &vector2) {
        return vector1.size() > vector2.size();
    });

    return groups[0];
}
