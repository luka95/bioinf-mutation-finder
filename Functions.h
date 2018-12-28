//
// Created by Sime on 10.12.2018..
//

#ifndef BIOINFORMATIKA_FUNCTIONS_H
#define BIOINFORMATIKA_FUNCTIONS_H

#include "string"
#include "vector"

using namespace std;

class Functions {

public:
    static vector<string> getKMers(const string &sequence, unsigned k);

};


#endif //BIOINFORMATIKA_FUNCTIONS_H
