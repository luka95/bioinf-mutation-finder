#include <iostream>
#include <string>
#include <set>
#include <algorithm>
#include <functional>
#include <vector>
#include "Index.h"
#include "LIS.h"

using namespace std;

typedef function<bool(pair<string, set<int>>, pair<string, set<int>>)> Comparator;
Comparator cmp = [](pair<string,set<int>> const &pair1, pair<string,set<int>> const &pair2){
    return *pair1.second.begin() < *pair2.second.begin();
};

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

    set<pair<string, set<int>>, Comparator> sorted_index(index.begin(), index.end(), cmp);

    for (const auto &it : sorted_index) {
        cout<< it.first <<" : ";
        for(const auto &pos_it : it.second){
            cout<<pos_it<<" ";
        }
        cout<<endl;
    }
}