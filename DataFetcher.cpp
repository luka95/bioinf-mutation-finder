//
// Created by Sime on 9.12.2018..
//

#include <string>
#include "DataFetcher.h"
#include <fstream>

using namespace std;

DataFetcher::DataFetcher(const string &genomePath, const string &mutatedPath) : genomePath(genomePath), mutatedPath(mutatedPath) {

}

void DataFetcher::loadData() {
    ifstream genomeFile (genomePath.c_str());
    string line;

    if(genomeFile.is_open())
    {
        cout << "Genome file opened";
        getline(genomeFile, line);
        getline(genomeFile, genome, string::traits_type::to_char_type(string::traits_type::eof()));
    }
    else {
        cerr << "Genome file couldnt be opened.";
        exit(-1);
    }
    genomeFile.close();

    ifstream mutatedFile (mutatedPath.c_str());

    if(mutatedFile.is_open()){
        cout << "Mutated file opened";
        while(getline(mutatedFile, line)){
            getline(mutatedFile, line);
            mutatedGenomeReads.push_back(line);
        }
    }
    else {
        cerr << "Mutated file couldnt be opened.";
        exit(-1);
    }
    mutatedFile.close();
}