//
// Created by divankovic on 9.1.2019..
//

#include <vector>
#include <algorithm>
#include "Mutation.h"
#include "Align.h"

using namespace std;

vector<Mutation> mutateSequence(string &sequence, double pm);
char symbols[] = {'A','C','G','T'};

int main(){
    string original = "ACGTA";
    string mutated = "AGGTA";

    vector<Mutation> mutations = align(original, mutated, 0);
    sort(mutations.begin(), mutations.end());
    for(auto &m : mutations) {
        cout << m;
    }

    return 0;
}

vector<Mutation> mutateSequence(string &sequence, double pm) {
    srand(unsigned(time(0)));
    vector<Mutation> mutations;
    int pos = 0;

    for(int i=0,n=sequence.length();i<n;i++) {
        double r = (double) rand() / RAND_MAX;

        if (r <= pm) {
            int mode = rand()%3;
            switch(mode){
                case 0: {
                    //supstitution
                    char replacement = symbols[rand() % 4];
                    sequence[pos] = replacement;
                    mutations.push_back(Mutation(MutationType::Substitution, i, replacement));
                    pos++;
                    break;
                }case 1: {
                    //deletion
                    sequence.erase(pos, 1);
                    mutations.push_back(Mutation(MutationType::Deletion, i, '-'));
                    break;
                }case 2: {
                    //insertion
                    string insertion(1, symbols[rand() % 4]);
                    sequence.insert(pos, insertion);
                    mutations.push_back(Mutation(MutationType::Insertion, i, insertion[0]));
                    pos+=2;
                    break;
                }default:
                    break;
            }
        }else{
            pos++;
        }
    }

    return mutations;
}