//
// Created by divankovic on 9.1.2019..
//

#include <string>
#include <iostream>
#include <algorithm>
#include <cassert>
#include "Inverter.h"

using namespace std;


string Inverter::inverse(string original) {
    transform(begin(original), end(original), begin(original), complement);
    reverse(original.begin(), original.end());
    return original;
}

char Inverter::complement(char c) {
    switch (c) {
        case 'A':
            return 'T';
        case 'T':
            return 'A';
        case 'G':
            return 'C';
        case 'C':
            return 'G';
        default:
            break;
    }
}
