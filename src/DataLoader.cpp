//
// Created by Sime on 9.12.2018..
//

#include <string>
#include "DataLoader.h"
#include <fstream>

using namespace std;

DataLoader::DataLoader(const string &genome_path, const string &mutated_path) : genome_path(genome_path), mutated_path(mutated_path) {

}

void DataLoader::loadData() {
    ifstream genomeFile (genome_path.c_str());
    string line;

    if(genomeFile.is_open())
    {
        cout << "Genome file opened";
        getline(genomeFile, line);
        while(getline(genomeFile, line)) {
            genome.append(line);
        }
    }
    else {
        cerr << "Genome file couldnt be opened.";
        exit(-1);
    }
    genomeFile.close();

    ifstream mutatedFile (mutated_path.c_str());

    if(mutatedFile.is_open()){
        cout << "Mutated file opened";
        while(getline(mutatedFile, line)){
            getline(mutatedFile, line);
            mutated_genome_reads.push_back(line);
        }
    }
    else {
        cerr << "Mutated file couldnt be opened.";
        exit(-1);
    }
    mutatedFile.close();
}