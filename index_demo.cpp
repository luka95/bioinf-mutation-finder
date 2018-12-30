#include <iostream>
#include <string>
#include <set>
#include "Index.h"

using namespace std;

/**
 * Demonstration of minimizer index.
 * @return 0
 * @author Dorian Ivanković
 */
int main() {
    string input = "231032101233101";
    int w = 3;
    int k = 3;
    unordered_map<string,set<int>> index = Index::buildMinimizerIndex(input, w, k);
    printf("Minimizers for string %s, w = %d, k = %d : \n", input.c_str(), w, k);

    for (const auto &it : index) {
        cout<< it.first <<" : ";
        for(const auto &pos_it : it.second){
            cout<<pos_it<<" ";
        }
        cout<<endl;
    }
}