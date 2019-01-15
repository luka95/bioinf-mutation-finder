//
// Created by bivankovic on 9.1.2019..
//
#include <string>


#ifndef BIOINF_MUTATION_FINDER_INVERTER_H
#define BIOINF_MUTATION_FINDER_INVERTER_H

using namespace std;

class Inverter {
public:
    /**
     * Returns the inverse complement of original string.
     * @param original - string to reverse
     * @return inverse complement
     */
    static string inverse(string original);
private:
    /**
     * A->T, G->C
     * @param c
     * @return
     */
    static char complement(char c);
};


#endif //BIOINF_MUTATION_FINDER_INVERTER_H
