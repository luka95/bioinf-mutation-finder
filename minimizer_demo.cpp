#include <iostream>
#include <string>
#include <set>
#include "Minimizer.h"

using namespace std;

/**
 * Demonstration of minimizers.
 * @return 0
 * @author Dorian Ivanković
 */
int main() {
    string input = "231032101233101";
    int w = 3;
    int k = 3;
    set<string> minimizers = Minimizer::getMinimizers(input, w, k);
    printf("Minimizers for string %s, w = %d, k = %d : \n", input.c_str(), w, k);

    for (const auto &minimizer : minimizers) {
        cout<< minimizer <<endl;
    }
}