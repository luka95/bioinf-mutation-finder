//
// Created by bivankovic on 9.1.2019..
//
#include <string>
#include <iostream>
#include "Inverter.h"

using namespace std;

int main(){
    string sequence = "ACCCGT";
    string complement = Inverter::inverse(sequence);
    cout<<complement;
}
