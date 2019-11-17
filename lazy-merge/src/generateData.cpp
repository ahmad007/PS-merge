#include "generateData.h"
#include <iostream>
#include <cmath>
#include <algorithm>


using namespace std;
/*This method generate a list of random elements for a given number of segments
 * @param arr the list to be created
 * @param n total list size
 * @param segmentsEnds list of each segment ending index
 * @param
 *
 */
void generateArray(vector<int >& arr, int n, const vector<int > segmentsEnds) {

  //  cout << "start generating data....\n";
    //generate all data
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

   // cout << "ending generating data.\n";
}
/* This method calculate the ending index of each segment so that all the segments are of equal sizes
 * @param size the total size
 * @param numOfSegments number of segments
 * @param segmentsEnds list of each segment ending index
 * @param
 *
 */
void setSegmentsEndEqually(int size, int numOfSegments, vector<int >& segmentsEnds) {
    unsigned segmentSize = ceil((double) size / numOfSegments);

    for (unsigned segmentIndex = 0, currentEnd = segmentSize; segmentIndex < numOfSegments; segmentIndex++) {
        segmentsEnds.push_back((currentEnd - 1 < size) ? currentEnd - 1 : size - 1);
        currentEnd += segmentSize;
    }
}

/*This method generate a list of random elements for a given number of segments for Guang's algorithm
 * @param arr the list to be created
 * @param n total list size
 * @param  numOfSegments number of segments
 * @param segmentsEnds list of each segment ending index
 * @param
 *
 */
void generateDataForGuang(vector<vector<int> >& arr, int n, int numOfSegments, vector<vector<int > >& segmentsEnds) {
    int segmentSize = n / numOfSegments;

    //each row has two segments
    int rowIndex = -1;
    for (rowIndex = 0; rowIndex < (numOfSegments / 2); rowIndex++) {
        vector<int > v;
        arr.push_back(v);
        vector<int > s;
        segmentsEnds.push_back(s);
    }
    //add the last segment if the number of segments is odd.
    if (numOfSegments % 2 == 1) {
        vector<int > v;
        arr.push_back(v);
        vector<int > s;
        segmentsEnds.push_back(s);
    }

    rowIndex = -1;
    for (int segmentIndex = 0; segmentIndex < numOfSegments; segmentIndex++) {
        if (segmentIndex == (numOfSegments - 1))
            segmentSize = n - (numOfSegments - 1) * segmentSize;

        if (segmentIndex % 2 == 0)
            rowIndex++;

        for (int i = 0; i < segmentSize; i++) {
            arr[rowIndex].push_back(std::rand() % 100);
        }
    }

    for (int i = 0; i < arr.size(); i++) {
        if(!(numOfSegments%2== 1&& i == arr.size()-1)){
            sort(arr[i].begin(), arr[i].begin() + arr[i].size() / 2);
            sort( arr[i].begin() + ( (arr[i].size() / 2)), arr[i].begin() + arr[i].size());           
            segmentsEnds[i].push_back(arr[i].size() / 2 -1);
            segmentsEnds[i].push_back(arr[i].size() - 1);
        }
    }
    if (numOfSegments % 2 == 1) {
        sort( arr[arr.size() - 1].begin(), arr[arr.size() - 1].begin() +  arr[arr.size() - 1].size());
        segmentsEnds[segmentsEnds.size() - 1].push_back(arr[arr.size() - 1].size() - 1);
    }

}
