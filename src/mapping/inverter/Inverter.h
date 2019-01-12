//
// Created by bivankovic on 9.1.2019..
//
#include <string>


#ifndef BIOINF_MUTATION_FINDER_INVERTER_H
#define BIOINF_MUTATION_FINDER_INVERTER_H

using namespace std;

class Inverter {
public:
    static string inverse(string original);
private:
    static char complement(char c);
};


#endif //BIOINF_MUTATION_FINDER_INVERTER_H
