//
// Created by dorian on 13/01/2019.
//

#ifndef BIOINF_MUTATION_FINDER_HIRSCHBERG_H
#define BIOINF_MUTATION_FINDER_HIRSCHBERG_H

#include <string>
#include <map>

using namespace std;

struct zw {
    string z;
    string w;
};

struct mutation {
    char type;
    char oldBase;
    char newBase;
};

struct result {
    int startIndex;
    int endIndex;
    map<int, mutation> mutations;
};

zw Hirschberg(string x, string y);
result getHirschbergAlignmentMutations(string x, string y, int inputOffset);

#endif //BIOINF_MUTATION_FINDER_HIRSCHBERG_H

