#include <iostream>
#include "Francis.h"
#include <omp.h>
#include <vector>
#include <algorithm>
#include <sys/time.h>

using namespace std;

/*
 * This method finds the partitions ends of each segment according to the number of threads
 * number of threads is used to determine the number of partitions
 * and each partition of equals number of elements.
 * @param a input vector
 * @param noOfSegments number of segments (k-segment) array
 * @param noOfThreads number of partitions
 * @param segmentsEnds list of each segment end
 * @param boundaries list of each parttion range in each segment, this list size equals number of partitions,
 * each element of this list is a vector of ranges represent this partiion.
 *
 */
void FrancisPartitioning(vector<int > a, const int noOfSegments, const int noOfThreads, const vector<int > segmentsEnds, vector<std::pair<int, vector<int> > >& boundaries) {
    int totalArrayLength = a.size();
    int partitionLen[noOfThreads]; // D == rank
    vector<int > b; //last index  + 1
    b.push_back(0);
    vector<double> times;

    for (int segmentIndex = 0; segmentsEnds.begin() + segmentIndex < segmentsEnds.end(); segmentIndex++)
        b .push_back(*(segmentsEnds.begin() + segmentIndex) + 1);

    vector< vector<int > > temps;
    for( int i = 0 ; i < noOfThreads; i++)
    {
    	vector<int> temp;
    	temps.push_back(temp);
    }

    #pragma omp parallel num_threads(noOfThreads)
    {
        int thread_id = omp_get_thread_num();

        int partitionLen = (totalArrayLength / noOfThreads) * (thread_id + 1);

        int deltaMax, deltaMin, delta;
            int rmin = 100, lmax = -10;
            int segContainLmax = -1, segContainRmin = 100;
            int temp;

            vector<int> s(noOfSegments, 0);
            vector<int> l(noOfSegments, 0);
            vector<int> r(noOfSegments, 0);

            for (int i = 0; i < noOfSegments; i++) {
                l[i] = b[i];
                r[i] = b[i + 1];
                temp = ((r[i] - l[i]) < (partitionLen / noOfSegments)) ? (r[i] - l[i]) : (partitionLen / noOfSegments);
                s[i] = l[i] + temp + ((i + 1 <= partitionLen % noOfSegments) ? 1 : 0);
            }



            while (true) {

                segContainLmax = -1;
                segContainRmin = 100;
                rmin = 100;
                lmax = -100;

                for (int i = 0; i < noOfSegments; i++) {
                    if (s[i] > l[i] && s[i] - 1 >= l[i] && s[i] - 1 < r[i] && a[s[i] - 1] >= lmax) {
                        lmax = a[s[i] - 1];
                        if (a[s[i] - 1] == lmax && i > segContainLmax)
                            segContainLmax = i;
                        else
                            segContainLmax = i;
                    }
                    if (s[i] < r[i] && s[i] >= l[i] && a[s[i]] <= rmin) {
                        rmin = a[s[i]];
                        if (a[s[i]] == rmin && i < segContainRmin)
                            segContainRmin = i;
                        else
                            segContainRmin = i;
                    }
                }

                if (s[segContainLmax] == l[segContainLmax] ||
                        s[segContainRmin] == r[segContainRmin] ||
                        lmax < rmin ||
                        (lmax == rmin && segContainLmax < segContainRmin)) {
                    break;
                }

                r[segContainLmax] = s[segContainLmax] - 1;
                l[segContainRmin] = s[segContainRmin] + 1;

                deltaMax = ((r[segContainLmax] - l[segContainLmax]) / 2);
                deltaMin = ((r[segContainRmin] - l[segContainRmin]) / 2);


                delta = (deltaMax < deltaMin) ? deltaMax : deltaMin;

                s[segContainLmax] = r[segContainLmax] - delta;
                s[segContainRmin] = l[segContainRmin] + delta;
            }


        for (int i = 0; i < noOfSegments; i++) {
                temps[thread_id].push_back(s[i] - 1);
            }

    }

    for( int i = 0 ; i < noOfThreads; i++)
    boundaries.push_back(std::make_pair(i, temps[i]));

    for (int i = 0; segmentsEnds.begin() + i < segmentsEnds.end(); i++)
        if ((boundaries.end() - 1)->second[i] != *(segmentsEnds.begin() + i))
            (boundaries.end() - 1)->second[i] = *(segmentsEnds.begin() + i);
}

