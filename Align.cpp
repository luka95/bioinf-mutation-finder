#include "Align.h"

const int D = -2;
const int MATCH = 4;
const int SUBSTITUTION = -1;

const int NUM_NEIGHBOUR_NODES = 3;

struct matrixNode {
    short value;
    matrixNode *src[NUM_NEIGHBOUR_NODES];
    int x;
    int y;
};

struct depthData {
    int depth;
    vector<short int> path;
    matrixNode* firstNode;
};

struct mutation {
    char type;
    char base;
};

struct result {
    int startIndex;
    int endIndex;
    vector<Mutation> mutations;
};

int checkIfMatches(char a, char b) {
    if (a==b) {
        return MATCH;
    }
    return SUBSTITUTION;
}

matrixNode* getMatrixNode(int x, int y, int value = 0) {
    auto node = new matrixNode;
    node->value = value;
    node->x = x;
    node->y = y;
    return node;
}

depthData getLongestPath(matrixNode* node, depthData data) {
    if (node->value == 0) {
        return data;
    } else {
        data.depth++;
        data.path.push_back(-1); // increase vector size

        for (short int i = 0; i < NUM_NEIGHBOUR_NODES; i++) {
            if (node->src[i] == nullptr) {
                continue;
            }
            data.path[data.depth] = i;
            data.firstNode = node;
            depthData newData = getLongestPath(node->src[i], data);
            if (newData.depth > data.depth) {
                data = newData;
            }
        }
        return data;
    }
}

depthData getLongestPath(matrixNode* node) {
    depthData data = {-1, vector<short int>()};
    data = getLongestPath(node, data);

    data.path.shrink_to_fit();

    if (data.path[data.depth] == -1) {
        data.path.pop_back();
    }

    return data;
}

/*
void printNode(matrixNode node) {
    cout << setw(10) << node.name << "=" << node.value << "[";
    for (int i = 0; i < NUM_NEIGHBOUR_NODES; i++) {
        if(node.src[i] == nullptr) {
            continue;
        }
        cout << node.src[i]->name;
    }
    cout << "]";
}
 */
/*
void printPath(matrixNode* startNode, depthData data){
    cout << "Longest path inversed:";
    cout << startNode->name;
    matrixNode* tmpNode = startNode;
    for (short int p : data.path) {
        tmpNode = tmpNode->src[p];
        cout << "-(" << p << ")->" << tmpNode->name;
    }
    cout << endl;

    cout << "first node: ";
    printNode(*tmpNode);
    cout << endl;
}

 */

vector<Mutation> align(string& mappedArea, string& mutatedGenomeRead, int regionOffset) {
    long s = mutatedGenomeRead.length();
    long t = mappedArea.length();

    // matrix initialization
    int M = 0;
    int Mi = 0;
    int Mj = 0;

    int counter = 0;
    //matrixNode * ary = new matrixNode[(s+1)*(t+1)];
    auto ** matrix = new matrixNode*[s+1];
    for(int i = 0; i < s+1; ++i){
        matrix[i] = new matrixNode[t+1];
        counter++;
    }
    matrix[0][0] = *getMatrixNode(0, 0, 0);
    for (int i = 1; i <= s; i++) {
        matrix[i][0] = *getMatrixNode(i, 0, i * D);
    }
    for (int j = 1; j <= t; j++) {
        matrix[0][j] = *getMatrixNode(0, j);
    }

    // matrix calculation
    for (int i = 1; i <= s; ++i) {
        for (int j = 1; j <= t ; ++j) {
            matrixNode* matchNode = &matrix[i-1][j-1];
            matrixNode* insertionNode = &matrix[i][j-1];
            matrixNode* deletionNode = &matrix[i-1][j];

            int match = matchNode->value + checkIfMatches(mutatedGenomeRead.at(i-1), mappedArea.at(j-1));
            int insertion = insertionNode->value + D;
            int deletion = deletionNode->value + D;

            int tmp[] = {match, insertion, deletion};
            int v = *std::max_element(tmp, tmp + NUM_NEIGHBOUR_NODES);

            matrixNode* newNode = getMatrixNode(i, j);

            newNode->value = v;

            if (match == v && matchNode->value > 0) {
                newNode->src[0] = matchNode;
            }
            if (insertion == v && insertionNode->value > 0) {
                newNode->src[1] = insertionNode;
            }
            if (deletion == v && deletionNode->value > 0) {
                newNode->src[2] = deletionNode;
            }

            matrix[i][j] = *newNode;
            //printNode(*newNode);
            if (v > M) {
                M = v;
                Mi = i;
                Mj =j;
            }
        }
        //cout << endl;
    }

    matrixNode* maxNode = &matrix[Mi][Mj];

    //cout << "MAX = " << M << endl;

    depthData depthData = getLongestPath(maxNode);

    vector<matrixNode*> nodesVector;
    nodesVector.push_back(maxNode);

    matrixNode* tmpNode = maxNode;
    for (short int p : depthData.path) {
        tmpNode = tmpNode->src[p];
        nodesVector.push_back(tmpNode);
    }

    //printPath(maxNode, depthData);

    vector<Mutation> mutations;
    tmpNode = maxNode;
    for (short int p : depthData.path) {
        tmpNode = tmpNode->src[p];
        if (p != 0) {
            auto mutationItem = new Mutation();
            if (p == 1) {
                mutationItem->type = Deletion;
                mutationItem->base = '-';
                mutationItem->position = tmpNode->x + regionOffset;
                mutations.push_back(*mutationItem);
            } else if (p == 2) {
                mutationItem->type = Insertion;
                mutationItem->base = mappedArea.at(tmpNode->y);
                mutationItem->position = tmpNode->x + regionOffset;
                mutations.push_back(*mutationItem);
            }
        }
    }
    for(int i = 0; i < s+1; ++i)
        delete [] matrix[i];
    delete [] matrix;
    return mutations;
}

