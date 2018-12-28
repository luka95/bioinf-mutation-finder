//
// Created by Sime on 9.12.2018..
//

#ifndef BIOINFORMATIKA_DATAFETCHER_H
#define BIOINFORMATIKA_DATAFETCHER_H

#include <iostream>
#include "vector"


using namespace std;

class DataFetcher {

public:
    string genome;
    vector<string> mutatedGenomeReads;
    string genomePath;
    string mutatedPath;


    DataFetcher(const string &genomePath, const string &mutatedPath);
    void loadData();

};


#endif //BIOINFORMATIKA_DATAFETCHER_H
