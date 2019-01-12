//
// Created by bivankovic on 2.1.2019..
//

#include <iostream>
#include <string>
#include <set>
#include <algorithm>
#include <functional>
#include <vector>
#include "Index.h"
#include "LIS.h"

using namespace std;

/**
 * Demonstration of LIS algorithm.
 * @return 0
 * @author Dorian Ivanković
 */
int main() {

    vector<int> positions;
    positions.push_back(10);
    positions.push_back(100);
    positions.push_back(22);
    positions.push_back(9);
    positions.push_back(33);
    positions.push_back(21);
    positions.push_back(50);
    positions.push_back(41);
    positions.push_back(60);
    positions.push_back(80);

    vector<int> lis = LIS::find(positions);
    for(const auto &it : lis){
        cout<<it<<endl;
    }

}