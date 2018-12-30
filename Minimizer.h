
#include <set>
#include <string>
using namespace std;

#ifndef BIOINF_MUTATION_FINDER_MINIMIZER_H
#define BIOINF_MUTATION_FINDER_MINIMIZER_H

class Minimizer {

public:
    static set<string> getMinimizers(string inputString, int w, int k);

private:
    static string getMinimizer(string window, int k);
    static set<string> getEndMinimizers(string inputString, int w, int k);
};


#endif //BIOINF_MUTATION_FINDER_MINIMIZER_H
