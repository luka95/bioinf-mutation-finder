//
// Created by bivankovic on 1.1.2019..
//

#include "LIS.h"

vector<int> LIS::find(vector<int> positions) {
    int len = positions.size();
    int T[len];
    int act_solution[len];
    vector<int> solution;

    for(int i=0;i<len;i++){
        T[i] = 1;
        act_solution[i] = i;
    }

    for(int i=1; i < len;i++){
        for(int j = 0;j<i;j++){
            if(positions[i] > positions[j]){
                if(T[j] + 1 > T[i]){
                    T[i] = T[j] + 1;
                    //point to position of the previous
                    act_solution[i] = j;
                }
            }
        }
    }

    int maxIndex = 0;
    for(int i=1;i<len;i++){
        if(T[i] > T[maxIndex]){
            maxIndex = i;
        }
    }

    //reconstruct actual sequence
    int t = maxIndex;
    int t1 = t;
    do {
        t = t1;
        solution.insert(solution.begin(), positions[t]);
        t1 = act_solution[t];
    }while(t != t1);

    return solution;
}
