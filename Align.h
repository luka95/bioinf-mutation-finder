//
// Created by Sime on 5.1.2019..
//

#include <array>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>
#include <cstring>
#include "Mutation.h"

#ifndef BIOINFORMATIKA_ALIGN_H
#define BIOINFORMATIKA_ALIGN_H

using namespace std;

vector<Mutation> align(string& mappedArea, string& mutatedGenomeRead, int regionOffset);

#endif //BIOINFORMATIKA_ALIGN_H

