//
// Created by Sime on 9.12.2018..
//

#ifndef BIOINFORMATIKA_DATAFETCHER_H
#define BIOINFORMATIKA_DATAFETCHER_H

#include <iostream>
#include "vector"


using namespace std;

class DataLoader {

public:
    string genome;
    vector<string> mutated_genome_reads;
    string genome_path;
    string mutated_path;


    DataLoader(const string &genome_path, const string &mutated_path);
    void loadData();

};


#endif //BIOINFORMATIKA_DATAFETCHER_H
