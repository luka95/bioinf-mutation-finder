#include "Minimizer.h"

set<string> Minimizer::getMinimizers(string inputString, int w, int k) {
    // length of window
    int l = w + k - 1;
    set<string> minimizers;

    for (int i = 0, n = inputString.length() - l + 1; i < n; i++) {
        string window = inputString.substr(i, l);
        minimizers.insert(getMinimizer(window, k));
    }

    set<string> endMinimizers = getEndMinimizers(inputString, w, k);
    minimizers.insert(endMinimizers.begin(), endMinimizers.end());

    return minimizers;
}

string Minimizer::getMinimizer(string window, int k) {
    string minimizer;

    for (int i = 0, n = window.length() - k + 1; i < n; i++) {
        string kmer = window.substr(i, k);
        if (minimizer.empty() || kmer.compare(minimizer) <= 0) {
            minimizer = kmer;
        }
    }

    return minimizer;
}

set<string> Minimizer::getEndMinimizers(string inputString, int w, int k) {
    set<string> endMinimizers;

    for (int u = 1; u <= w - 1; u++) {
        int l = u + k - 1;
        string start = inputString.substr(0,l);
        string end = inputString.substr(inputString.length()-l,l);

        endMinimizers.insert(getMinimizer(start, k));
        endMinimizers.insert(getMinimizer(end,k));
    }


    return endMinimizers;
}
