
#include <string>
#include <set>
#include <unordered_map>
#include <vector>

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

    static const int INDEX_HIT_MARGIN = 500;
    static const int MINIMUM_MINIMIZER_MATCHES = 4;
    static const int MINIMUM_MATCH_LENGTH = 40;

    /**
     * double strand minimizers
     * @param inputString
     * @param w
     * @param k
     * @return
     */
    static unordered_map<string, set<tuple<int, int>>> index(string &inputString, int w, int k);

    /**
     * The method extracts all minimizers from the input string using a window of size w
     * and k-mers of size k. The algorithm is described in https://academic.oup.com/bioinformatics/article/20/18/3363/202143
     * @param inputString - string to extract minimizers from
     * @param w - window size
     * @param k - kmer's size
     * @return minimizer index
     */
    static unordered_map<string, set<int>> buildMinimizerIndex(string &inputString, int w, int k);

    /**
    * Finds the best match between reference and sequence given their indexes
    * @param reference_index
    * @param sequence_index
    * @return ((int, int, int, int) , int) reference begin, reference end, sequence begin, sequence end + strand_xor(same strands = 0, different strands = 1)
    */
    static tuple<tuple<int, int, int, int>, int>
    getBestMatch(unordered_map<string, set<tuple<int, int>>> &reference_index,
                 unordered_map<string, set<tuple<int, int>>> &sequence_index);

private:
    /**
     *
     * @param window - window to extract minimizer from
     * @param k - kmer's size
     * @return minimizer of window and its offset in a window
     */
    static tuple<string, int> getMinimizer(string &window, int k);

    /**
     * Extracts end minimizers from the input string.
     * @param inputString - string to extract end minimizers from
     * @param w - window size
     * @param k - kmer's size
     * @return end minimizers
     */
    static unordered_map<string, set<int>> getEndMinimizers(string &inputString, int w, int k);

    /**
     * Extract the double strand minimizer from the window of minimizers and returns its hash
     * @param window - window of minimizers
     * @param k - minimizer's length
     * @return double strand minimizer for window
     */
    static int getDoubleStrandMinimizer(string &window, int k);


    static unsigned long long hash(string &sequence);
};


#endif //BIOINF_MUTATION_FINDER_MINIMIZER_H
