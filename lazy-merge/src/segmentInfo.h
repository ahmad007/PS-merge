/* 
 * File:   segmentInfo.h
 * Author: ahmad
 *
 * Created on December 26, 2013, 3:12 PM
 */

#ifndef SEGMENTINFO_H
#define	SEGMENTINFO_H

#include <vector>
#include <iostream>
#include <cstdlib>
#include <map>
#include <omp.h>
#include "IndexInfoCache.h"

using namespace std;

struct segmentInfo {
public:

    vector<std::pair<int, int> > boundries;
    int* first; //actual start
    int* second; //actual end
    int* virtualStart;
    int* virtualEnd;
    int noOfSegments;
    int* virtualIndexRanges;

    void init(int numOfSegments) {
        first = second = 0;
        second = (int*) malloc(numOfSegments * sizeof (int)); //new int [numOfSegments];
        first = (int*) malloc(numOfSegments * sizeof (int)); //new int [numOfSegments];
        virtualIndexRanges = (int*) malloc((numOfSegments + 1) * sizeof (int)); //new int [numOfSegments];
        virtualStart = (int*) malloc(numOfSegments * sizeof (int));
        virtualEnd = (int*) malloc(numOfSegments * sizeof (int));
        noOfSegments = numOfSegments;
    }

    void finalize() {
        free(first); //delete [] first;
        free(second); //delete [] second;
        free(virtualIndexRanges);
        free(virtualStart);
        free(virtualEnd);
    }

    void getSize() {

        size = 0;
        virtualIndexRanges[0] = 0;
        for (int i = 0; i < noOfSegments; i++) {
            size += second[i] - first[i] + 1;
            virtualIndexRanges[i + 1] = size;
        }

        //        if (omp_get_thread_num() == 0) {
        //            cout << "virtual indexes sizes: ";
        //            for (int i = 0; i <= noOfSegments; i++) {
        //                cout << virtualIndexRanges[i] << ",";
        //            }
        //            cout << endl;
        //        }
    }

    int convertVirtualIndextoActual(int virtualIndex) {

        //         cout << virtualIndex << ",";               
        int actualIndex = -1;
        int start = 0, end = noOfSegments, current;

        virtualIndex = virtualIndex - first[0];


        while ((end - start) > 1) {
            current = (start + end) / 2;
            if (virtualIndex > virtualIndexRanges[current])
                start = current;
            else
                end = current;
        }

        if (virtualIndex == virtualIndexRanges[end])
            current = end;
        else
            current = start;
        actualIndex = first[current] + (virtualIndex - virtualIndexRanges[current]);

        //        if (first[0] == 0)
        //            cout << virtualIndex << "*" << current << ">>" << first[current] + (virtualIndex - virtualIndexRanges[current]) << "true: ";


        //        int i = 0;
        //        virtualIndex = virtualIndex - first[i];
        //        for (; i < noOfSegments; i++) {
        //            if (virtualIndex < second[i] - first[i] + 1) {
        //                actualIndex = first[i] + virtualIndex;
        //                break;
        //            }
        //            virtualIndex -= second[i] - first[i] + 1;
        //        }
        //
        //        if (first[0] == 0)
        //           cout << actualIndex << endl;

        return actualIndex;
    }

    int convertVirtualIndextoActual(int virtualIndex, IndexCache& cache) {

        //         cout << virtualIndex << ",";               
        int actualIndex = -1;
        int start = 0, end = noOfSegments, current;

        virtualIndex = virtualIndex - first[0];


        while ((end - start) > 1) {
            current = (start + end) / 2;
            if (virtualIndex > virtualIndexRanges[current])
                start = current;
            else
                end = current;
        }

        if (virtualIndex == virtualIndexRanges[end])
            current = end;
        else
            current = start;
        actualIndex = first[current] + (virtualIndex - virtualIndexRanges[current]);

        cache.actualStart[0] = first[current];
        cache.actualEnd[0] = second[current];
        cache.virtualStart[0] = first[0] + virtualIndexRanges[current];
        cache.virtualEnd[0] = first[0] + virtualIndexRanges[current + 1] - 1;

        return actualIndex;
    }
    unsigned size;
};

//        typedef vector<std::pair<int, int> > my_vector;
//        my_vector::iterator it = boundries.begin();
//        virtualIndex = virtualIndex - it->first;
//
//        for (; it < boundries.end(); it++) {
//            if (virtualIndex < it->second - it->first + 1) {
//                actualIndex = it->first + virtualIndex;
//                break;
//            }
//            virtualIndex -= it->second - it->first + 1;
//        }


#endif	/* SEGMENTINFO_H */

