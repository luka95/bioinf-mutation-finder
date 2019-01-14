#include <array>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <iterator>
#include <chrono>
#include <ctime>

using namespace std;

const int DELETION = -2;
const int INSERTION = -2;
const int MATCH = 2;
const int SUBSTITUTION = -1;

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
        firstLine[j] = firstLine[j-1] + INSERTION;
    }

    for (int i = 1; i <= xlen; i++) {
        secondLine[0] = firstLine[0] + DELETION;

        for (int j = 1; j <= ylen; j++) {
            int scoreSub = firstLine[j-1] + sub(x.at(i - 1), y.at(j - 1));
            int scoreDel = firstLine[j] + DELETION;
            int scoreIns = secondLine[j-1] + INSERTION;
            secondLine[j] = maximum(scoreSub, scoreIns, scoreDel);
        }

        if (i < xlen) {
            copy(secondLine, secondLine + ylen + 1, firstLine);
        }
    }
    delete[] firstLine;

    return secondLine;
}

tuple<string, string> NeedlemanWunsch(string x, string y) {
    unsigned long xlen = x.length();
    unsigned long ylen = y.length();

    string z;
    string w;

    if (xlen == 1 && ylen == 1) {
        z = x;
        w = y;
    }
    else {
        string::size_type n;

        if (xlen == 1) {
            n = y.find(x);

            if (n == string::npos) {
                z = x + string((ylen-1), '-');
            }
            else {
                z = string(n, '-') + x + string(ylen - n - 1, '-');
            }
            w = y;
        }
        else if (ylen == 1) {
            n = x.find(y);

            z = x;
            if (n == string::npos) {
                w = y + string((xlen -1), '-');
            }
            else {
                w = string(n, '-') + y + string(xlen - n - 1, '-');
            }
        }
    }

    return {z, w};
}

tuple<string, string> Hirschberg(string x, string y) {
    unsigned long xlen = x.length();
    unsigned long ylen = y.length();

    string z;
    string w;

    if (xlen == 0) {
        for (int i = 1; i <= ylen; i++) {
            z = z + '-';
            w = w + y.at(i - 1);
        }
    } else if (ylen == 0) {
        for (int i = 1; i <= xlen; i++) {
            z = z + x.at(i - 1);
            w = w + '-';
        }
    } else if (xlen == 1 || ylen == 1) {
        auto zw =  NeedlemanWunsch(x, y);
        z = get<0>(zw);
        w = get<1>(zw);

    } else {
        int xmid = xlen / 2;
        int* scoreLeft = NWScore(x.substr(0, xmid), y);
        int* scoreRight = NWScore(revString(x.substr(xmid)), revString(y));

        reverse(scoreRight, &scoreRight[ylen + 1]);

        int ymid = maxBetweenTwoLines(scoreLeft, scoreRight, ylen);

        delete[] scoreLeft;
        delete[] scoreRight;

        auto leftZW = Hirschberg(x.substr(0, xmid), y.substr(0, ymid));
        auto rightZW = Hirschberg(x.substr(xmid), y.substr(ymid));

        z = z + get<0>(leftZW) + get<0>(rightZW);
        w = w + get<1>(leftZW) + get<1>(rightZW);
    }

    return {z, w};
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