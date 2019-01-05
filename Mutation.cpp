//
// Created by Sime on 31.10.2018..
//

#include "Mutation.h"

Mutation::Mutation() {}

Mutation::Mutation(MutationType type, int position, char base) : type(type), position(position), base(base){}

bool Mutation::operator<(const Mutation &other) const {
    return position < other.position;
}

ostream& operator<<(ostream& strm, const Mutation& mt) {
    return strm << (char)mt.type << " " << mt.position << " " << mt.base << "\n";
}