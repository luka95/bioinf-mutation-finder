#include "Index.h"

unordered_map<string, set<int>> Index::buildMinimizerIndex(string& inputString, int w, int k) {
    // length of window
    int l = w + k - 1;
    unordered_map<string,set<int>> index;

    for (int i = 0, n = inputString.length() - l + 1; i < n; i++) {
        string window = inputString.substr(i, l);

        tuple<string,int> tup = getMinimizer(window, k);
        string minimizer = get<0>(tup);
        int offset = i + get<1>(tup);

        index[minimizer].insert(offset);
    }

    unordered_map<string,set<int>> endMinimizers = getEndMinimizers(inputString, w, k);

    for(const auto &it : endMinimizers){
        set<int> endIndices = it.second;
        for(const auto &idx : endIndices){
            index[it.first].insert(idx);
        }
    }

    return index;
}

tuple<string,int> Index::getMinimizer(string& window, int k) {
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

unordered_map<string, set<int>> Index::getEndMinimizers(string& inputString, int w, int k) {
    unordered_map<string,set<int>> endMinimizersIndex;
    int n = inputString.length();

    for (int u = 1; u <= w - 1; u++) {
        int l = u + k - 1;
        string start = inputString.substr(0,l);
        string end = inputString.substr(n-l,l);

        tuple<string,int> tup = getMinimizer(start, k);
        endMinimizersIndex[get<0>(tup)].insert(get<1>(tup));

        tup = getMinimizer(end,k);
        endMinimizersIndex[get<0>(tup)].insert(n-l+get<1>(tup));
    }


    return endMinimizersIndex;
}


