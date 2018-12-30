
#include <string>
#include <set>
#include <unordered_map>

using namespace std;

#ifndef BIOINF_MUTATION_FINDER_MINIMIZER_H
#define BIOINF_MUTATION_FINDER_MINIMIZER_H

/**
 * Class that contains methods for extracting minimizers which are used to replace
 * all k-mers in a string.
 * Implemented as described in https://academic.oup.com/bioinformatics/article/20/18/3363/202143
 * @author Dorian Ivanković
 */
class Index {

public:
    /**
     * The method extracts all minimizers from the input string using a window of size w
     * and k-mers of size k. The algorithm is described in https://academic.oup.com/bioinformatics/article/20/18/3363/202143
     * @param inputString - string to extract minimizers from
     * @param w - window size
     * @param k - kmer's size
     * @return minimizer index
     */
    static unordered_map<string, set<int>> buildMinimizerIndex(string& inputString, int w, int k);

private:
    /**
     *
     * @param window - window to extract minimizer from
     * @param k - kmer's size
     * @return minimizer of window and its offset in a window
     */
    static tuple<string,int> getMinimizer(string& window, int k);

    /**
     * Extracts end minimizers from the input string.
     * @param inputString - string to extract end minimizers from
     * @param w - window size
     * @param k - kmer's size
     * @return end minimizers
     */
    static unordered_map<string, set<int>> getEndMinimizers(string& inputString, int w, int k);
};


#endif //BIOINF_MUTATION_FINDER_MINIMIZER_H
