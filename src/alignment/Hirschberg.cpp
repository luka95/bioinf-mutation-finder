#include <array>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <iterator>

using namespace std;

const int D = -2;
const int MATCH = 2;
const int SUBSTITUTION = -1;

struct zw {
    string z;
    string w;
};

struct mutation {
    char type;
    char oldBase;
    char newBase;
};

struct result {
    int startIndex;
    int endIndex;
    map<int, mutation> mutations;
};

int ins(char y) {
    return D;
}

int del(char x) {
    return D;
}

int sub(char x, char y) {
    if (x == y) {
        return MATCH;
    } else {
        return SUBSTITUTION;
    }
}

int maximum(int a, int b, int c) {
    int max = ( a < b ) ? b : a;
    max = ( max < c ) ? c : max;

    return max;
}

string revString(string original) {
    reverse(original.begin(), original.end());
    return original;
}

int maxBetweenTwoLines(const int* firstLine, const int* secondLine, int length) {
    int max = firstLine[0] + secondLine[0];
    int maxIndex = 0;

    for (int i = 1; i < length; i++) {
        int tmp = firstLine[i] + secondLine[i];
        if (tmp > max) {
            max = tmp;
            maxIndex = i;
        }
    }
    return maxIndex;
}

int* NWScore(string x, string y) {
    unsigned long xlen = x.length();
    unsigned long ylen = y.length();
    int* firstLine = new int[ylen + 1];
    int* secondLine = new int[ylen + 1];

    firstLine[0] = 0;
    for (int j = 1; j <= ylen; j++) {
        firstLine[j] = firstLine[j-1] + ins(y.at(j - 1));
    }

    for (int i = 1; i <= xlen; i++) {
        secondLine[0] = firstLine[0] + del(x.at(i - 1));

        for (int j = 1; j <= ylen; j++) {
            int scoreSub = firstLine[j-1] + sub(x.at(i - 1), y.at(j - 1));
            int scoreDel = firstLine[j] + del(x.at(i - 1));
            int scoreIns = secondLine[j-1] + ins(y.at(j - 1));
            secondLine[j] = maximum(scoreSub, scoreIns, scoreDel);
        }

        if (i < xlen) {
            copy(secondLine, secondLine + ylen + 1, firstLine);
        }
    }
    delete[] firstLine;

    return secondLine;
}

int getAlignmentStartIndex(string w) {
    int start = 0;
    for(char& c : w) {
        if (c == '-') {
            start++;
        } else {
            break;
        }
    }
    return start;
}

int getAlignmentEndIndex(string w) {
    int end = w.length() - 1;
    for(auto c = w.rbegin(); c != w.rend(); ++c) {
        if (*c == '-') {
            end--;
        } else {
            break;
        }
    }
    return end;
}

map<int, mutation> getMutations(zw alignment, int startIndex, int endIndex, int inputOffset) {
    map<int, mutation> mutations;

    // TODO paralelize
    for(int i = startIndex; i <= endIndex; i++) {
        char x = alignment.z.at(i);
        char y = alignment.w.at(i);
        if(x == y) {
            continue;
        } else {
            auto mutationItem = new struct mutation;
            mutationItem->oldBase = x;
            mutationItem->newBase = y;

            if (x == '-') {
                mutationItem->type = 'I';
            } else if (y == '-') {
                mutationItem->type = 'D';
            } else { // x != y
                mutationItem->type = 'S';
            }

            mutations[i + inputOffset] = *mutationItem;
        }
    }

    return mutations;
}

zw NeedlemanWunsch(string x, string y) {
    unsigned long xlen = x.length();
    unsigned long ylen = y.length();

    zw alignment;

    if (xlen == 1 && ylen == 1) {
        alignment.z = x;
        alignment.w = y;
    }
    else {
        string::size_type n;

        if (xlen == 1) {
            n = y.find(x);

            if (n == string::npos) {
                alignment.z = x + string((ylen-1), '-');
            }
            else {
                alignment.z = string(n, '-') + x + string(ylen - n - 1, '-');
            }
            alignment.w = y;
        }
        else if (ylen == 1) {
            n = x.find(y);

            alignment.z = x;
            if (n == string::npos) {
                alignment.w = y + string((xlen -1), '-');
            }
            else {
                alignment.w = string(n, '-') + y + string(xlen - n - 1, '-');
            }
        }
    }

    return alignment;
}

zw Hirschberg(string x, string y) {
    zw alignment;
    alignment.z = "";
    alignment.w = "";

    unsigned long xlen = x.length();
    unsigned long ylen = y.length();

    if (xlen == 0) {
        for (int i = 1; i <= ylen; i++) {
            alignment.z = alignment.z + '-';
            alignment.w = alignment.w + y.at(i - 1);
        }
    } else if (ylen == 0) {
        for (int i = 1; i <= xlen; i++) {
            alignment.z = alignment.z + x.at(i - 1);
            alignment.w = alignment.w + '-';
        }
    } else if (xlen == 1 || ylen == 1) {
        alignment = NeedlemanWunsch(x, y);
    } else {
        int xmid = xlen / 2;
        int* scoreLeft = NWScore(x.substr(0, xmid), y);
        int* scoreRight = NWScore(revString(x.substr(xmid)), revString(y));

        reverse(scoreRight, &scoreRight[ylen + 1]);

        int ymid = maxBetweenTwoLines(scoreLeft, scoreRight, ylen);

        delete[] scoreLeft;
        delete[] scoreRight;

        zw leftZW = Hirschberg(x.substr(0, xmid), y.substr(0, ymid));
        zw rightZW = Hirschberg(x.substr(xmid), y.substr(ymid));

        alignment.z = alignment.z + leftZW.z + rightZW.z;
        alignment.w = alignment.w + leftZW.w + rightZW.w;
    }

    return alignment;
}

result getHirschbergAlignmentMutations(string x, string y, int inputOffset) {
    zw alignment = Hirschberg(x, y);

    int alignmentStartIndex = inputOffset + getAlignmentStartIndex(alignment.w);
    int alignmentEndIndex = inputOffset + getAlignmentEndIndex(alignment.w);

    map<int, mutation> mutations = getMutations(alignment, alignmentStartIndex, alignmentEndIndex, inputOffset);

    auto res = new struct result;
    res->startIndex = alignmentStartIndex;
    res->endIndex = alignmentEndIndex;
    res->mutations = mutations;

    cout << "Alignment:" << endl;
    cout << ">" << alignment.z << endl;
    cout << ">" << alignment.w << endl;

    return *res;
}

//int main() {
//
//    string x = "AGTACGCA";
//    string y = "TATGC";
//
//    result res = getHirschbergAlignmentMutations(x, y, 1000);
//
//    return 0; // TODO return res
//}