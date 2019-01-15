//
// Created by dorian on 1.1.2019..
//

#ifndef BIOINF_MUTATION_FINDER_LIS_H
#define BIOINF_MUTATION_FINDER_LIS_H

using namespace std;
#include <vector>

class LIS {
public:

    /**
     * Finds the longest increasing subsequence in positions
     * @param positions
     * @return longest increasing subsequence
     */
    static vector<int> find(vector<int> positions);

    /**
     * Finds the longest increasing subsequnce of second elements in vector
     * @param vec
     * @param strand_xor
     * @return
     */
    static vector<tuple<int,int>> findBySecond(vector<tuple<int,int>> vec, int strand_xor);
};


#endif //BIOINF_MUTATION_FINDER_LIS_H
