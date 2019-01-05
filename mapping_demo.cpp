//
// Created by bivankovic on 5.1.2019..
//
#include <string>
#include <iostream>
#include "Index.h"

using namespace std;

int main(){
    string reference = "AAGATCCCATGCCATGCCTGAAGTCCCTGAAT";
    string sequence = "GTGCAATGCCTGAAT";
    auto reference_index = Index::buildMinimizerIndex(reference, 5,5);
    auto sequence_index = Index::buildMinimizerIndex(sequence, 5,5);
    tuple<int,int> positions = Index::getBestMatch(reference_index, sequence_index);
    cout<<reference.substr(get<0>(positions), get<1>(positions) - get<0>(positions)+1);
}
