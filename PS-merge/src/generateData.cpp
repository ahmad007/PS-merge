#include "generateData.h"
#include <iostream>
#include <cmath>
#include <algorithm>


using namespace std;
void generateArray(vector<int >& arr, int n, const vector<int > segmentsEnds) {
    for (int i = 0; i < n; i++) {
        arr.push_back(std::rand() % 100);
    }

    int start = 0;
    //sort each segment separately
    for (int segmentIndex = 0; segmentsEnds.begin() + segmentIndex < segmentsEnds.end(); segmentIndex++) {
        if (segmentIndex > 0)
            start = *(segmentsEnds.begin() + segmentIndex - 1) + 1;

        sort(arr.begin() + start, arr.begin() + *(segmentsEnds.begin() + segmentIndex) + 1);
    }
}

void generateArray(vector<int >& arr, int n, int numOfPortions, const vector<int > segmentsEnds) {
    for (int i = 0; i < n; i++) {
        arr.push_back(std::rand() % 100);
    }

    int start = 0;
    //sort each segment separately
    for (int segmentIndex = 0; segmentsEnds.begin() + segmentIndex < segmentsEnds.end(); segmentIndex++) {
        if (segmentIndex > 0)
            start = *(segmentsEnds.begin() + segmentIndex - 1) + 1;
       
        sort(arr.begin() + start, arr.begin() + *(segmentsEnds.begin() + segmentIndex) + 1);
    }
}

void setSegmentsEndEqually(int size, int numOfSegments, vector<int >& segmentsEnds) {
    unsigned segmentSize = ceil((double) size / numOfSegments);

    for (unsigned segmentIndex = 0, currentEnd = segmentSize; segmentIndex < numOfSegments; segmentIndex++) {
        segmentsEnds.push_back((currentEnd - 1 < size) ? currentEnd - 1 : size - 1);
        currentEnd += segmentSize;
    }
}
