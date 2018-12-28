//
// Created by Sime on 10.12.2018..
//

#include "Functions.h"

vector<string> Functions::getKMers(const string &sequence, unsigned k) {
    vector<string> kmers;

    for(unsigned i=0 ; i<sequence.size() - k + 1; ++i) {
        kmers.push_back(sequence.substr(i, k));
    }
    return kmers;
}