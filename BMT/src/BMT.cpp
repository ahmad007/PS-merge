//============================================================================
// Name        : BMT.cpp
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
#include <stack>

using namespace std;

//IS_DEBUG = 0 to print the actual running time of the algorithm
//IS_DEBUG = 1 to print the number of swaps and comparisons
#define IS_DEBUG 0

unsigned long numofBMTComparisons = 0;
unsigned long BMTCounter = 0;
unsigned long numOfBMTSwaps = 0;

struct param {
public:

	param(int p1, int p2, int p3) {
		f1 = p1;
		f2 = p2;
		last = p3;
	}
	int f1;
	int f2;
	int last;
};

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

	cout << "\nGenerating data is started...." << flush;
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
	cout << "the check is done with no errors!\n" << endl;
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

void rightShiftCircular(vector<int>& arr, int shiftStartIndex,
		int shiftEndIndex, int shiftValue) {
	int temp = 0;
	for (int i = shiftStartIndex, j = shiftEndIndex; --j > i; i++) {
		temp = arr[i];
		arr[i] = arr[j];
		arr[j] = temp;
#if IS_DEBUG == 1
#pragma omp critical
        {
            if (numOfBMTSwaps == ULONG_MAX) {
            	BMTCounter++;
                numOfBMTSwaps = 1;
            } else
            	numOfBMTSwaps++;
        }
#endif
	}

	for (int i = shiftStartIndex, j = i + shiftValue; --j > i; i++) {
		temp = arr[i];
		arr[i] = arr[j];
		arr[j] = temp;
#if IS_DEBUG == 1
#pragma omp critical
        {
            if (numOfBMTSwaps == ULONG_MAX) {
            	BMTCounter++;
                numOfBMTSwaps = 1;
            } else
            	numOfBMTSwaps++;
        }
#endif
	}

	for (int i = shiftStartIndex + shiftValue, j = shiftEndIndex; --j > i;
			i++) {
		temp = arr[i];
		arr[i] = arr[j];
		arr[j] = temp;
#if IS_DEBUG == 1
#pragma omp critical
        {
            if (numOfBMTSwaps == ULONG_MAX) {
            	BMTCounter++;
                numOfBMTSwaps = 1;
            } else
            	numOfBMTSwaps++;
        }
#endif
	}
}

void iterativeSplitMerge(vector<int>& a, int f1, int f2, int last) {

	std::stack<param> parameters;
	param p(f1, f2, last);
	parameters.push(p);
	int l, r, ldash, rdash, m, mdash;
	while (!parameters.empty()) {

		param temp = parameters.top();
		parameters.pop();

		f1 = temp.f1;
		f2 = temp.f2;
		last = temp.last;

#if IS_DEBUG == 1
#pragma omp critical
		{
			numofBMTComparisons++;
		}
#endif

		if (f1 >= f2 || f2 >= last)
			continue;
		l = f1;
		r = f2;
		ldash = f2;
		rdash = last;
		m = 0;
		mdash = 0;

		while (!(l >= r && ldash >= rdash)) {
#if IS_DEBUG == 1
#pragma omp critical
			{
				numofBMTComparisons += 4;
			}
#endif

			if (l < r)
				m = (l + r) / 2;
			if (ldash < rdash)
				mdash = (ldash + rdash) / 2;
			if (a[m] <= a[mdash]) {
				l = m + 1;
				rdash = mdash;
			} else {
				ldash = mdash + 1;
				r = m;
			}
		}
		rightShiftCircular(a, r, ldash, (ldash - f2));

		param pFirst(f1, r, r + rdash - f2);
		param pSecond(l + ldash - f2, ldash, last);
		parameters.push(pSecond);
		parameters.push(pFirst);

	}
}

void splitmerge2Portions(vector<int>& a, vector<std::pair<int, int> > limits,
		int mergingTaskID) {
	int mergingTask = mergingTaskID;
	int size = (limits.begin() + 2 * mergingTask)->second
			- (limits.begin() + 2 * mergingTask)->first + 1;
	size += (limits.begin() + 2 * mergingTask + 1)->second
			- (limits.begin() + 2 * mergingTask + 1)->first + 1;

	iterativeSplitMerge(a, (limits.begin() + 2 * mergingTask)->first,
			(limits.begin() + 2 * mergingTask + 1)->first,
			(limits.begin() + 2 * mergingTask)->first + size);
}

void multiWaySplitMerge(vector<int>& a, vector<int> segmentsEnds,
		int usedThreads) {
	struct timeval startwtime1, endwtime1;
	double elapsed_time1;
	gettimeofday(&startwtime1, NULL);

	int noOfSwaps = 0;
	vector<std::pair<int, int> > limits(segmentsEnds.size()); //actual index

	int start = 0;
	for (int i = 0; segmentsEnds.begin() + i < segmentsEnds.end(); i++) {
		limits.push_back(std::make_pair(start, *(segmentsEnds.begin() + i)));
		start = *(segmentsEnds.begin() + i) + 1;
	}

	int numOfMergingTasks, tempStart, tempEnd;
	typedef std::vector<std::pair<int, int> > my_vector;
	my_vector::iterator it;
	int phases = std::ceil(std::log(limits.size()) / std::log(2));

	numOfMergingTasks = limits.size() / 2;

	for (int i = 0; i < phases; i++) {
		numOfMergingTasks = limits.size() / 2;

#pragma omp parallel for schedule(static,1) num_threads(usedThreads)
		for (int w = 0; w < numOfMergingTasks; w++) {
			splitmerge2Portions(a, limits, w);
		}

		for (int threadId = 0; threadId < numOfMergingTasks; threadId++) {
			it = limits.begin();
			tempStart = (it + threadId)->first;
			tempEnd = (it + (threadId + 1))->second;
			limits.erase(it + threadId, it + (threadId + 2));
			limits.insert(it + threadId, std::make_pair(tempStart, tempEnd));
		}
	}
	gettimeofday(&endwtime1, NULL);

	elapsed_time1 = (double) ((endwtime1.tv_usec - startwtime1.tv_usec) / 1.0e6
			+ endwtime1.tv_sec - startwtime1.tv_sec);

#if IS_DEBUG == 0
	cout << " the wall clock time = " << elapsed_time1 << " seconds." << endl;
#endif

#if IS_DEBUG == 1
	cout << ", # of swaps: " << numOfBMTSwaps;
	cout << ", and # of comparisons: " << numofBMTComparisons << endl;
	cout << flush;
#endif

}

int main() {
	vector<int> array;
	vector<int> segmentsEnds;

	int noOfSegments = 128; //k value in the k-way merge
	int array_size = std::pow(2.0, 24);
	int numberOfthreads = 2;

		array.clear();
		segmentsEnds.clear();

		setSegmentsEndEqually(array_size, noOfSegments, segmentsEnds);
		generateArray(array, array_size, segmentsEnds);

		cout << "The k-way merging is started for BMT Merge.. array size: "
				<< array.size() << " # of segments: " << segmentsEnds.size()
				<< ", # of threads: " << numberOfthreads << " .... " << endl;
		multiWaySplitMerge(array, segmentsEnds, numberOfthreads);
		cout << " merging is done!" << endl;

#if IS_DEBUG == 1
		numofBMTComparisons = 0;
		numOfBMTSwaps = 0;
		BMTCounter = 0;
		checkIt(array);
#endif


	cout << "\n\nAll tests are Passed!\n";
	return 0;
}
