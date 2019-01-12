//
// Created by dorian on 1.1.2019..
//

#ifndef BIOINF_MUTATION_FINDER_LIS_H
#define BIOINF_MUTATION_FINDER_LIS_H

using namespace std;
#include <vector>

class LIS {
public:

    static vector<int> find(vector<int> positions);
    static vector<tuple<int,int>> findBySecond(vector<tuple<int,int>> vec, int strand_xor);
};


#endif //BIOINF_MUTATION_FINDER_LIS_H
