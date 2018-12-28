//
// Created by Sime on 31.10.2018..
//

#ifndef BIOINFORMATIKA_MUTATION_H
#define BIOINFORMATIKA_MUTATION_H

#include "string"
#include <iostream>
#include <fstream>

using namespace std;

enum MutationType { Substitution = 'X', Insertion = 'I', Deletion = 'D' };

class Mutation {

public:
    MutationType type;
    int position;
    char base;

    Mutation(MutationType type, int position, char base);
private:
    friend ostream& operator<<(ostream&, const Mutation& mt);

};


#endif //BIOINFORMATIKA_MUTATION_H
