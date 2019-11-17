//============================================================================
// Name        : Bitonic-Merge.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <sys/time.h>
#include <omp.h>
#include <limits.h>

using namespace std;

//IS_DEBUG = 0 to print the actual running time of the algorithm
//IS_DEBUG = 1 to print the number of swaps and comparisons

#define IS_DEBUG 0

unsigned long numofBitonicComparisons = 0;
unsigned long bitonicCounter = 0;
unsigned long numOfBitonicSwaps = 0;

void setSegmentsEndEqually(int size, int numOfSegments,
		vector<int>& segmentsEnds) {

	unsigned segmentSize = ceil((double) size / numOfSegments);

	for (unsigned segmentIndex = 0, currentEnd = segmentSize;
			segmentIndex < numOfSegments; segmentIndex++) {
		segmentsEnds.push_back(
				(currentEnd - 1 < size) ? currentEnd - 1 : size - 1);
		currentEnd += segmentSize;
	}
}

void generateArray(vector<int>& arr, int n, const vector<int> segmentsEnds) {

	cout << "Generating data is started...." << flush;
	//generate all data
	for (int i = 0; i < n; i++) {
		arr.push_back(std::rand() % 100);
	}

	int start = 0;
	//sort each segment separately
	for (int segmentIndex = 0;
			segmentsEnds.begin() + segmentIndex < segmentsEnds.end();
			segmentIndex++) {
		if (segmentIndex > 0)
			start = *(segmentsEnds.begin() + segmentIndex - 1) + 1;

		sort(arr.begin() + start,
				arr.begin() + *(segmentsEnds.begin() + segmentIndex) + 1);
	}

	cout << " generating data is done!\n";
}

bool checkIt(vector<int> a) {
	cout << "Checking if array is merged correctly started..." << flush;
	for (int i = 1; i < a.size(); i++)
		if (a[i] < a[i - 1]) {
			cout << "error>>>>>>>>" << a[i] << "," << a[i - 1] << endl;
			return false;
		}
	cout << "the check is done with no errors!" << endl;
	return true;
}

int getNumOfPhases(int x) {
	int counter = 0;
	while (x > 0) {
		counter++;
		x /= 2;
	}
	return counter;
}

void reverseIt(vector<int>& a, int startIndex, int endIndex, int numOfThreads) {

	int s = ((endIndex - startIndex - 1) / 2) / numOfThreads;

#pragma omp parallel num_threads(numOfThreads)
	{
		int id = omp_get_thread_num();
		int size = ((endIndex - startIndex - 1) / 2) / numOfThreads;

		int i = startIndex + id * size;
		int j = endIndex - id * size;

		int endingIndex = (i + size);
		if (id == numOfThreads - 1)
			endingIndex = ((endIndex - startIndex + 1)) / 2 + startIndex;

		for (; i < j && i < endingIndex; i++, j--) {
			int temp = a[i];
			a[i] = a[j];
			a[j] = temp;
#if IS_DEBUG == 1
#pragma omp critical
			{
				if (numOfBitonicSwaps == ULONG_MAX) {
					bitonicCounter++;
					numOfBitonicSwaps = 1;
				} else
				numOfBitonicSwaps++;
			}
#endif
		}
	}
}

/* This method performs the 2-way bitonic merging
 * @param a input list of k-segment
 * @param firstSegmentBoundaries first segment range
 * @param secondSegmentBoundaries second segment range
 * @param numOfAvailableThreads number of used threads
 *
 */
void bitonicMerge(vector<int>& a, std::pair<int, int> firstSegmentBoundaries,
		std::pair<int, int> secondSegmentBoundaries, int numOfThreads) {

	int startIndex = firstSegmentBoundaries.first;
	int size = firstSegmentBoundaries.second - firstSegmentBoundaries.first + 1;
	size += secondSegmentBoundaries.second - secondSegmentBoundaries.first + 1;

//reversing the second segment
	reverseIt(a, startIndex + (size / 2), startIndex + (size - 1),
			numOfThreads);

	int numOfPhases = getNumOfPhases(size);
	int temp;

	for (int phase = 1; phase < numOfPhases; phase++) {
		int shift = size / (std::pow(2.0, phase));
#pragma omp parallel num_threads(numOfThreads) shared(shift)
		{
			int threadID = omp_get_thread_num();

			int blockSize = (int) (std::ceil((size / 2.0) / numOfThreads));
			int comparatorIDStart = threadID * (blockSize);
			int comparatorIDEnd =
					(((threadID + 1) * (blockSize) - 1) > (size / 2 - 1)) ?
							(size / 2 - 1) : ((threadID + 1) * (blockSize) - 1);

			for (int comparatorID = comparatorIDStart;
					comparatorID <= comparatorIDEnd; comparatorID++) {
				int blockID = (comparatorID / shift);
				int entryIndex = (blockID * 2 * shift)
						+ (comparatorID - (blockID * shift));

				if (a[startIndex + entryIndex]
						> a[startIndex + entryIndex + shift]) {
					int temp = a[startIndex + entryIndex];
					a[startIndex + entryIndex] = a[startIndex + entryIndex
							+ shift];
					a[startIndex + entryIndex + shift] = temp;
#if IS_DEBUG == 1
#pragma omp critical
					{
						if (numOfBitonicSwaps == ULONG_MAX) {
							bitonicCounter++;
							numOfBitonicSwaps = 1;
						} else
						numOfBitonicSwaps++;
					}
#endif

				}
#if IS_DEBUG == 1
#pragma omp critical
				{
					numofBitonicComparisons++;
				}
#endif

			}
		}
	}
}

/* This method performs the k-way bitonic merging
 * @param a input list of k-segment
 * @param segmentsEnds list of each segment end
 * @param numOfAvailableThreads number of used threads
 *
 */
void mutliwayBitonicMerge(vector<int>& a, vector<int> segmentsEnds,
		int numOfAvailableThreads) {

	vector<std::pair<int, int> > limits; //actual index
	int size = a.size();

	for (int i = 0; i < segmentsEnds.size(); i++) {
		if (i == 0)
			limits.push_back(std::make_pair(0, segmentsEnds[i]));
		else
			limits.push_back(
					std::make_pair(segmentsEnds[i - 1] + 1, segmentsEnds[i]));
	}

	int numOfSegments = limits.size();

	int noOfMergingTasks, tempStart, tempEnd;

	int phases = std::ceil(std::log(limits.size()) / std::log(2));

	struct timeval startwtime1, endwtime1;
	double elapsed_time1;
	gettimeofday(&startwtime1, NULL);

	for (int i = 0; i < phases; i++) {

		int noOfMergingTasks = limits.size() / 2;

		for (int mergingTaskIndex = 0; mergingTaskIndex < noOfMergingTasks;
				mergingTaskIndex++) {
			bitonicMerge(a, limits[2 * mergingTaskIndex],
					limits[2 * mergingTaskIndex + 1], numOfAvailableThreads);
		}

		for (int segmentIndex = 0; segmentIndex < limits.size();
				segmentIndex++) {
			if (segmentIndex == limits.size() - 1)
				break;

			tempStart = limits[segmentIndex].first;
			tempEnd = limits[segmentIndex + 1].second;

			for (int i = tempStart + 1; i < tempEnd; i++)
				if (a[i] < a[i - 1])
					cout << "error\n";

			limits.erase(limits.begin() + segmentIndex,
					limits.begin() + (segmentIndex + 2));
			limits.insert(limits.begin() + segmentIndex,
					std::make_pair(tempStart, tempEnd));
		}
	}

	gettimeofday(&endwtime1, NULL);

	elapsed_time1 = (double) ((endwtime1.tv_usec - startwtime1.tv_usec) / 1.0e6
			+ endwtime1.tv_sec - startwtime1.tv_sec);
#if IS_DEBUG == 0
	cout << " the wall clock time = " << elapsed_time1 << " seconds." << endl;
#endif

#if IS_DEBUG == 1
	cout << " # of swaps: " << numOfBitonicSwaps;
	cout << ", and # of comparisons: " << numofBitonicComparisons << endl;
	cout << flush;
#endif
}

int main() {
	vector<int> array;
	vector<int> segmentsEnds;

	//parameters tuning
	int noOfSegments = 128; //k value in the k-way merge
	int array_size = std::pow(2.0, 24);
	int numberOfthreads = 2;

	if (!((array_size != 0) && ((array_size & (array_size - 1)) == 0))) {
		cerr << "Error: The array array_size is not a power of two." << endl;
		return 0;
	}

	setSegmentsEndEqually(array_size, noOfSegments, segmentsEnds);
	generateArray(array, array_size, segmentsEnds);
	cout << "The k-way merging is started for Bitonic Merge.. array size: "
			<< array.size() << ", # of segments: " << segmentsEnds.size()
			<< ", # of threads: " << numberOfthreads << " .... " << endl;
	mutliwayBitonicMerge(array, segmentsEnds, numberOfthreads);
	cout << " merging is done!\n" << endl;

#if IS_DEBUG == 1
	numofBitonicComparisons = 0;
	bitonicCounter = 0;
	numOfBitonicSwaps = 0;
	checkIt(array);
#endif

	cout << "\n\nAll tests are Passed!\n";
	return 0;
}
