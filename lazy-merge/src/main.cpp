/* 
 * File:   main.cpp
 * Author: ahmad
 *
 * Created on December 2, 2013, 2:13 PM
 */

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include<parallel/algorithm>
#include <vector>
#include <limits>
#include <ctime>
#include <sys/time.h>
#include <limits.h>
#include "Francis.h"
#include "splitMerge.h"
#include "splitMergeOnDistance.h"
#include "generateData.h"

//linker-> libraries-> libraries (-l) ->  boost_serialization

//#include <boost/archive/text_oarchive.hpp>
//#include <boost/archive/text_iarchive.hpp>
//#include <boost/serialization/vector.hpp>
vector<int> realIndex;

using namespace std;

typedef uint16_t _ThreadIndex;

ofstream myfile;

void printTime() {
	time_t now;
	struct tm *current;
	now = time(0);
	current = localtime(&now);
	cout << current->tm_hour << ":" << current->tm_min << ":" << current->tm_sec
			<< " millisecond: " << clock() << endl;
}

bool lazy_merge(int size, unsigned numOfSegments, int noOfThreads) {
	vector<int> a;

	vector<int> segmentsEnds; //actual index

	//actual index + 1, if it contains index less than segment start then
	//this segment in not included
	vector<std::pair<int, vector<int> > > boundaries;

	struct timeval startwtime1, endwtime1;
	double seq_time1;

	//generating list
	cout << "\nGenerating data is started...." << flush;
	setSegmentsEndEqually(size, numOfSegments, segmentsEnds);
	generateArray(a, size, segmentsEnds);
	cout << "Done!\n";

	vector<segmentInfo> segmentsInfos;
	gettimeofday(&startwtime1, NULL);

	//create partitions
	cout << "Merging (Partitioning step) is started..." << flush;
	FrancisPartitioning(a, numOfSegments, noOfThreads, segmentsEnds,
			boundaries);
	cout << "Done!" << endl;

	cout << "Merging (tournament step) is started..." << flush;
	vector<int> currentSegmentStarts;

	//set segments starts
	currentSegmentStarts.push_back(0);
	for (int i = 0; i < segmentsEnds.size() - 1; i++) {
		currentSegmentStarts.push_back(*(segmentsEnds.begin() + i) + 1);
	}

	for (int i = 0; i < boundaries.size(); i++) {
		//initializing the segment_info
		int tempSize = 0;
		for (int b = 0; b < (boundaries.begin() + i)->second.size(); b++) {
			if (*((boundaries.begin() + i)->second.begin() + b)
					>= *(currentSegmentStarts.begin() + b)) {
				tempSize++;
			}
		}
		segmentInfo s1;
		s1.init(tempSize);

		for (int b = 0, k = 0; b < (boundaries.begin() + i)->second.size();
				b++) {
			if (*((boundaries.begin() + i)->second.begin() + b)
					>= *(currentSegmentStarts.begin() + b)) {
				s1.first[k] = *(currentSegmentStarts.begin() + b);
				s1.second[k++] =
						*((boundaries.begin() + i)->second.begin() + b);
				currentSegmentStarts[b] =
						*((boundaries.begin() + i)->second.begin() + b) + 1;
			}
		}
		s1.getSize();
		segmentsInfos.push_back(s1);
	}

	long comparisons = -1;
	multiwaySplitMergeOnDistanceCache(a, segmentsInfos, comparisons);
	gettimeofday(&endwtime1, NULL);
	cout << "Done!" << endl;

	seq_time1 = (double) ((endwtime1.tv_usec - startwtime1.tv_usec) / 1.0e6
			+ endwtime1.tv_sec - startwtime1.tv_sec);

	cout << "Lazy-Merge size:" << size << ", # of seg: " << numOfSegments
			<< ", # of threads: " << noOfThreads << flush;

#if IS_DEBUG == 0
	cout << ", wall clock time = " << seq_time1 << " seconds." << endl;
#endif

#if IS_DEBUG == 1
	cout << ", # of swaps: " << getNumOfSwapsSplitMergeOnDistance() << endl;

	if (!checkArray(realIndex, a, segmentsInfos)) {
		cout << "errror!!!\n";
		return false;
	}

#endif
	for (int i = 0; i < segmentsInfos.size(); i++)
		(segmentsInfos.begin() + i)->finalize();

	return true;
}

int main(int argc, char** argv) {

	int array_size = 50000000; //std::pow(2.0, 24);

	int numberOfthreads = 2;
	int noOfSegments = 64;
	lazy_merge(array_size, noOfSegments, numberOfthreads);

	cout << "\n\nAll tests are Passed!\n";
	return 0;
}

