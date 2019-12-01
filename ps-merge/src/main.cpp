/* 
 * File:   main.cpp
 * Author: ahmad
 *
 * Created on July 27, 2014, 3:48 PM
 */

#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <ctime>
#include <sys/time.h>
#include <sstream>
#include <string>
#include "generateData.h"
#include "Francis.h"
#include "partitionInfo.h"
#include "splitMerge.h"

using namespace std;

bool debug1 = false;
bool fineGrained = true;
bool parallel = true;
//bool numOfThreads = -1;
int numOfThreads = -1;
ofstream threadLog;
ofstream overallLog;
double accuSavedTime = 0;
double numOfSwaps = 0;
double mergeSwaps = 0;
unsigned long mergeComparsons = 0;
double shuffleSwaps = 0;

double shuffleTime = 0;
double mergeTime = 0;
double partTime = 0;
double partitioningTime = 0;
double elapsed_time;

bool checkIt(vector<int> a) {
	for (int i = 1; i < a.size(); i++)
		if (a[i] < a[i - 1]) {
			cout << "error>>>>>>>>" << a[i] << "," << a[i - 1] << endl;
			return false;
		}
	return true;
}

/**
 * 
 * @param arr
 * @param firstPartStartIndex inclusive
 * @param secondPartStartIndex inclusive
 */

int isPowerOfTwo(unsigned int x) {
	return ((x != 0) && !(x & (x - 1)));
}

void cleanShuffle(vector<int>& arr) {
	int temp;
	for (int i = 2; i < arr.size(); i += 2)
		if (arr[i] < arr[i - 1]) {
			temp = arr[i];
			arr[i] = arr[i - 1];
			arr[i - 1] = temp;
		}
}

void rightShiftCircular2(vector<int>& arr, int shiftStartIndex,
		int shiftEndIndex, int shiftValue) {
	const clock_t begin_time = clock();
	int temp = 0;
	int batchSize = (shiftEndIndex - shiftStartIndex) / numOfThreads;
#pragma omp parallel num_threads(numOfThreads)
	// for (int i = 0 ; i  < numOfThreads; i++)
	{
		int x = omp_get_thread_num();

		int start, end;
		start = shiftStartIndex + x * batchSize;
		if (x == (numOfThreads - 1)) {
			end = ((shiftStartIndex + (x + 1) * batchSize) > shiftEndIndex) ?
					shiftEndIndex : (shiftStartIndex + (x + 1) * batchSize);

		} else
			end = shiftStartIndex + (x + 1) * batchSize;

		for (int i = start, j = end - 1; i < (start + (batchSize / 2));
				i++, j--) {
			temp = arr[i];
			arr[i] = arr[j];
			arr[j] = temp;
		}
#if IS_DEBUG == 1
#pragma omp critical
		{
			numOfSwaps++;
		}
#endif
	}

	batchSize = ((shiftStartIndex + shiftValue) - shiftStartIndex)
			/ numOfThreads;
#pragma omp parallel num_threads(numOfThreads)
	// for (int i = 0 ; i  < numOfThreads; i++)
	{
		int x = omp_get_thread_num();

		int start, end;
		start = shiftStartIndex + x * batchSize;
		if (x == (numOfThreads - 1)) {
			end = ((shiftStartIndex + (x + 1) * batchSize) > shiftEndIndex) ?
					shiftEndIndex : (shiftStartIndex + (x + 1) * batchSize);

		} else
			end = shiftStartIndex + (x + 1) * batchSize;

		for (int i = start, j = end - 1; i < (start + (batchSize / 2));
				i++, j--) {
			temp = arr[i];
			arr[i] = arr[j];
			arr[j] = temp;
		}
#if IS_DEBUG == 1
#pragma omp critical
		{
			numOfSwaps++;
		}
#endif
	}

	batchSize = (shiftEndIndex - (shiftStartIndex + shiftValue)) / numOfThreads;
#pragma omp parallel num_threads(numOfThreads)
	// for (int i = 0 ; i  < numOfThreads; i++)
	{
		int x = omp_get_thread_num();

		int start, end;
		start = (shiftStartIndex + shiftValue) + x * batchSize;
		if (x == (numOfThreads - 1)) {
			end = (((shiftStartIndex + shiftValue) + (x + 1) * batchSize)
					> shiftEndIndex) ?
					shiftEndIndex :
					((shiftStartIndex + shiftValue) + (x + 1) * batchSize);

		} else
			end = (shiftStartIndex + shiftValue) + (x + 1) * batchSize;

		for (int i = start, j = end - 1; i < (start + (batchSize / 2));
				i++, j--) {
			temp = arr[i];
			arr[i] = arr[j];
			arr[j] = temp;

		}
#if IS_DEBUG == 1
#pragma omp critical
		{
			numOfSwaps++;
		}
#endif
	}
}

void rightShiftCircular1(vector<int>& arr, int shiftStartIndex,
		int shiftEndIndex, int shiftValue) {

	int temp = 0;
	for (int i = shiftStartIndex, j = shiftEndIndex; --j > i; i++) {
		temp = arr[i];
		arr[i] = arr[j];
		arr[j] = temp;
#if IS_DEBUG == 1
#pragma omp critical
		{
			numOfSwaps++;
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
			numOfSwaps++;
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
			numOfSwaps++;
		}
#endif
	}
}
int gcd(int a, int b) {
	if (b == 0)
		return a;
	else
		return gcd(b, a % b);
}

void par_rightShiftCircular2(vector<int>& arr, int shiftStartIndex,
		int shiftEndIndex, int shiftValue) {
	struct timeval start1, end1;
	gettimeofday(&start1, NULL);

	int n = shiftEndIndex - shiftStartIndex;
	int k = n - shiftValue;

	int d = -1, temp, j;
	int max = gcd(n, k);
	int threadSize = ((max / numOfThreads) == 0) ? 1 : max / numOfThreads;

	if (max >= numOfThreads) {
#pragma omp parallel for  schedule(static ,threadSize) num_threads(numOfThreads)
		for (int i = 0; i < max; i++) {
			int j = i;
			int temp = arr[shiftStartIndex + i];
			int d = 0;

			while (1) {
				d = (j + k) % n;
				if (d == i)
					break;
#if IS_DEBUG == 1
#pragma omp critical
				{
					numOfSwaps++;
					shuffleSwaps++;
				}
#endif
				arr[shiftStartIndex + j] = arr[shiftStartIndex + d];
				j = d;
			}

			arr[shiftStartIndex + j] = temp;
#if IS_DEBUG == 1
#pragma omp critical
			{
				numOfSwaps++;
				shuffleSwaps++;
			}
#endif
		}
	} else {

		vector<int> thread_temp_value(numOfThreads);
		vector<int> thread_temp_index(numOfThreads);

		vector<int> cycle_temp(max);

		for (int threadID = 0; threadID < numOfThreads; threadID++) {
			thread_temp_index[threadID] = -1;
		}

		for (int y = 0; y < max; y++)
			cycle_temp[y] = -1;

		vector<int> thread_iterations(numOfThreads);

		int threadSize = (n / numOfThreads == 0) ? 1 : n / numOfThreads;

		for (int i = 0; i < numOfThreads; i++) {
			thread_iterations[i] = threadSize;
		}

		if (numOfThreads * threadSize != n)
			thread_iterations[numOfThreads - 1] = n
					- (numOfThreads - 1) * threadSize;

#pragma omp parallel num_threads(numOfThreads)
		//	for (int threadID = 0; threadID < numOfThreads; threadID++)
		{
			int threadID = omp_get_thread_num();

			// cycle i
			int i = (threadID * threadSize) / (n / max);
			//iteration j of cycle i
			int j = ((threadID * threadSize) % (n / max));
			j = i + ((long) j * k) % n;

			if (j == i) {
				cycle_temp[i] = arr[shiftStartIndex + j];
			} else {
				thread_temp_value[threadID] = arr[shiftStartIndex + j];
				thread_temp_index[threadID] = j;
			}
			int d = -1;
			for (int x = 1; x < thread_iterations[threadID]; x++) {
				d = (j + k) % n;
				//move to a new cycle
				if (d == i) {
					i = i + 1;
					cycle_temp[i] = arr[shiftStartIndex + i];
					j = i;
				} else {
					arr[shiftStartIndex + j] = arr[shiftStartIndex + d];
					j = d;
				}
			}
		}

		//connecting the element within cycles and
		for (int threadID = 0; threadID < numOfThreads; threadID++) {
			if (thread_temp_index[threadID] != -1) {
				int j = thread_temp_index[threadID] - k;
				if (j < 0)
					j += n;
				arr[shiftStartIndex + j] = thread_temp_value[threadID];
			}
		}

		//connecting cyclyes
		for (int cycle = 0; cycle < max; cycle++) {
			int j = cycle - k;
			if (j < 0)
				j += n;
			arr[shiftStartIndex + j] = cycle_temp[cycle];
		}

	}
}

void par_rightShiftCircular1(vector<int>& arr, int shiftStartIndex,
		int shiftEndIndex, int shiftValue) {

	int temp = 0;
	int j = shiftEndIndex;
	int endIndex =
			((shiftEndIndex - shiftStartIndex - 1) % 2 == 0) ?
					((shiftEndIndex - shiftStartIndex - 1) / 2) - 1 :
					((shiftEndIndex - shiftStartIndex - 1) / 2);

	int threadSize =
			((endIndex / numOfThreads) == 0) ? 1 : endIndex / numOfThreads;

#pragma omp parallel for private( j,temp) schedule(static ,threadSize) num_threads(numOfThreads)
	for (int i = shiftStartIndex; i <= shiftStartIndex + endIndex; i++) {

		j = shiftEndIndex - (i - shiftStartIndex + 1);
		temp = arr[i];
		arr[i] = arr[j];
		arr[j] = temp;
#if IS_DEBUG == 1
#pragma omp critical
		{
			numOfSwaps++;
			shuffleSwaps++;
		}
#endif
	}

	j = shiftStartIndex + shiftValue;
	endIndex =
			((shiftStartIndex + shiftValue - shiftStartIndex - 1) % 2 == 0) ?
					((shiftStartIndex + shiftValue - shiftStartIndex - 1) / 2)
							- 1 :
					((shiftStartIndex + shiftValue - shiftStartIndex - 1) / 2);
	threadSize = ((endIndex / numOfThreads) == 0) ? 1 : endIndex / numOfThreads;

#pragma omp parallel for private( j,temp) schedule(static ,threadSize) num_threads(numOfThreads)
	for (int i = shiftStartIndex; i <= shiftStartIndex + endIndex; i++) {

		j = shiftStartIndex + shiftValue - (i - shiftStartIndex + 1);

		temp = arr[i];
		arr[i] = arr[j];
		arr[j] = temp;
#if IS_DEBUG == 1
#pragma omp critical
		{
			numOfSwaps++;
			shuffleSwaps++;
		}
#endif
	}

	j = shiftEndIndex;
	endIndex =
			((shiftEndIndex - (shiftStartIndex + shiftValue) - 1) % 2 == 0) ?
					((shiftEndIndex - (shiftStartIndex + shiftValue) - 1) / 2)
							- 1 :
					((shiftEndIndex - (shiftStartIndex + shiftValue) - 1) / 2);
	threadSize = ((endIndex / numOfThreads) == 0) ? 1 : endIndex / numOfThreads;

#pragma omp parallel for private( j,temp) schedule(static ,threadSize) num_threads(numOfThreads)
	for (int i = shiftStartIndex + shiftValue;
			i <= shiftStartIndex + shiftValue + endIndex; i++) {
		j = shiftEndIndex - (i - (shiftStartIndex + shiftValue) + 1);
		temp = arr[i];
		arr[i] = arr[j];
		arr[j] = temp;
#if IS_DEBUG == 1
#pragma omp critical
		{
			numOfSwaps++;
			shuffleSwaps++;
		}
#endif
	}
}

bool shuffleRanges(vector<int>& arr, vector<std::pair<int, int> > & twoLists) {

	int listSize = twoLists.size();

	vector<int> swapPositions;
	vector<int> tempSwapPositions;
	swapPositions.push_back(listSize / 2 - 1);
	int integerPart, level = 1, swapSize;
	float swapSizeFl;
	int levelsNo = 0;
	int v = listSize;

	struct timeval start1, end1;
	gettimeofday(&start1, NULL);

	while (v >>= 1) // unroll for more speed...
	{
		levelsNo++;
	}
	if (isPowerOfTwo(listSize))
		levelsNo--;

	int count = 0;

	do {

		if (numOfThreads == 1) {
			break;
		}

		tempSwapPositions.clear();
		swapSizeFl = listSize / pow(2.0, level + 1);
		integerPart = (swapSizeFl - (int) swapSizeFl) * 100;

		if (integerPart > 50)
			swapSize = (int) swapSizeFl + 1;
		else
			swapSize = (int) swapSizeFl;

		if (integerPart == 50 && swapSize > 0) {
			count++;
			if (count > 1) {
				swapSize++;
			}
		}

		if (swapSize == 0)
			swapSize = 1;

		int shiftStartPosition, shiftEndPosition, shiftValue;

		//do the swapping for each position
		for (int i = 0; i < swapPositions.size(); i++) {

			int secondPartStartRangeIndex = -1, secondPartEndRangeIndex = -1; //exclude range of -1,-1 from the two ends
			int firstPartStartRangeIndex = -1, firstPartEndRangeIndex = -1; //exclude range of -1,-1 from the two ends

			for (int rID = (swapPositions[i] + swapSize);
					rID >= (swapPositions[i] + 1); rID--)
				if (twoLists[rID].first != -1) {
					secondPartEndRangeIndex = rID;
					break;
				}

			for (int rID = (swapPositions[i] + 1);
					rID <= (swapPositions[i] + swapSize); rID++)
				if (twoLists[rID].first != -1) {
					secondPartStartRangeIndex = rID;
					break;
				}

			for (int rID = swapPositions[i];
					rID >= (swapPositions[i] - swapSize + 1); rID--)
				if (twoLists[rID].first != -1) {
					firstPartEndRangeIndex = rID;
					break;
				}

			for (int rID = (swapPositions[i] - swapSize + 1);
					rID <= swapPositions[i]; rID++)
				if (twoLists[rID].first != -1) {
					firstPartStartRangeIndex = rID;
					break;
				}

			shiftStartPosition = twoLists[firstPartStartRangeIndex].first;
			shiftEndPosition = twoLists[secondPartEndRangeIndex].second + 1;

			if (firstPartEndRangeIndex == -1 && secondPartEndRangeIndex == -1) {
				continue;
			}

			if (firstPartEndRangeIndex != -1 && secondPartEndRangeIndex != -1) {
				int reducingValue = twoLists[firstPartEndRangeIndex].second
						- twoLists[firstPartStartRangeIndex].first + 1;
				shiftValue = twoLists[secondPartEndRangeIndex].second
						- twoLists[secondPartStartRangeIndex].first + 1;

				//Right circular Shifty by the right side value
				if (!parallel || !fineGrained) {
					rightShiftCircular1(arr, shiftStartPosition,
							shiftEndPosition, shiftValue);
				} else {
					int n = shiftEndPosition - shiftStartPosition;
					int k = n - shiftValue;
					int max = gcd(n, k);
					//cout<< n << ", k= " << k << " max = " << max << " ";
					int threadSize =
							((max / numOfThreads) == 0) ?
									1 : max / numOfThreads;

					if (max >= numOfThreads) {
						par_rightShiftCircular2(arr, shiftStartPosition,
								shiftEndPosition, shiftValue);
					} else
						par_rightShiftCircular1(arr, shiftStartPosition,
								shiftEndPosition, shiftValue);

				}

				//updating ranges;
				std::pair<int, int> temp;
				for (int firstPartIndex = (swapPositions[i] - swapSize + 1),
						secondPartIndex = (swapPositions[i] + 1);
						firstPartIndex <= swapPositions[i];
						firstPartIndex++, secondPartIndex++) {

					if (twoLists[firstPartIndex].first != -1) {
						twoLists[firstPartIndex].first += shiftValue;
						twoLists[firstPartIndex].second += shiftValue;
					}

					if (twoLists[secondPartIndex].first != -1) {
						twoLists[secondPartIndex].first -= reducingValue;
						twoLists[secondPartIndex].second -= reducingValue;
					}
					//swap pairs
					temp = twoLists[firstPartIndex];
					twoLists[firstPartIndex] = twoLists[secondPartIndex];
					twoLists[secondPartIndex] = temp;
				}
			} //end of both parts has non dummy ranges
			else {
				//if one side has all dummy values
				//just swap parts
				std::pair<int, int> temp;
				for (int firstPartIndex = (swapPositions[i] - swapSize + 1),
						secondPartIndex = (swapPositions[i] + 1);
						firstPartIndex <= swapPositions[i];
						firstPartIndex++, secondPartIndex++) {

					//swap pairs
					temp = twoLists[firstPartIndex];
					twoLists[firstPartIndex] = twoLists[secondPartIndex];
					twoLists[secondPartIndex] = temp;
				}
			}
		}

		//update positions, replacing each position by new two positions
		for (int i = 0; i < swapPositions.size(); i++) {
			tempSwapPositions.push_back(swapPositions[i] - swapSize);
			tempSwapPositions.push_back(swapPositions[i] + swapSize);
		}
		swapPositions.clear();
		for (int i = 0; i < tempSwapPositions.size(); i++) {
			swapPositions.push_back(tempSwapPositions[i]);
		}
		level++;
	} while (level <= levelsNo);

	gettimeofday(&end1, NULL);
	shuffleTime += (double) ((end1.tv_usec - start1.tv_usec) / 1.0e6
			+ end1.tv_sec - start1.tv_sec);

	struct timeval startwtime2, endwtime2;
	gettimeofday(&startwtime2, NULL);
	if (!parallel || fineGrained == false) {
		threadLog << "new tasks:\n";
		for (int i = 0; i < twoLists.size(); i += 2) {
			struct timeval startwtime1, endwtime1;
			double seq_time1;
			double swaps = 0;

			if (twoLists[i].first != -1 && twoLists[i + 1].first != -1)
				iterativeSplitMerge(arr, twoLists[i].first,
						twoLists[i + 1].first, twoLists[i + 1].second + 1,
						swaps);
#if IS_DEBUG == 1
#pragma omp critical
			{
				mergeSwaps += swaps;
				mergeComparsons += getNumOfComparisonsSplitMerge();
			}
#endif
			gettimeofday(&endwtime1, NULL);
			seq_time1 = (double) ((endwtime1.tv_usec - startwtime1.tv_usec)
					/ 1.0e6 + endwtime1.tv_sec - startwtime1.tv_sec);
			threadLog << " time in seconds: " << seq_time1 << endl;
		}
	} else {

		threadLog << "new tasks:\n";
		//	cout<<"used threads here>>>>" << twoLists.size()/2 <<endl;
		resetNumOfComparisonsSplitMerge();
#pragma omp parallel num_threads((twoLists.size()/2))
		//for (int i = 0; i < twoLists.size(); i += 2)
		{
			int i = omp_get_thread_num();
			double swaps = 0;
			struct timeval startwtime1, endwtime1;
			double seq_time1;
			//gettimeofday(&startwtime1, NULL);

			i *= 2;
			if (twoLists[i].first != -1 && twoLists[i + 1].first != -1)
				iterativeSplitMerge(arr, twoLists[i].first,
						twoLists[i + 1].first, twoLists[i + 1].second + 1,
						swaps);
#if IS_DEBUG == 1
#pragma omp critical
			{
				mergeSwaps += swaps;
				mergeComparsons += getNumOfComparisonsSplitMerge();
			}
#endif
		}
	}

	gettimeofday(&endwtime2, NULL);

	mergeTime += (double) ((endwtime2.tv_usec - startwtime2.tv_usec) / 1.0e6
			+ endwtime2.tv_sec - startwtime2.tv_sec);

	//update ranges
	for (int i = 0; i < twoLists.size(); i++) {
		if (twoLists[i].first != -1 && twoLists[i + 1].first != -1) {
			twoLists[i].second = twoLists[i + 1].second;
			twoLists.erase(twoLists.begin() + (i + 1));
		} else {
			twoLists[i].second =
					(twoLists[i].second > twoLists[i + 1].second) ?
							twoLists[i].second : twoLists[i + 1].second;
			twoLists[i].first =
					(twoLists[i].first > twoLists[i + 1].first) ?
							twoLists[i].first : twoLists[i + 1].first;
			twoLists.erase(twoLists.begin() + (i + 1));
		}
	}
	return true;
}

bool kwayMerge(vector<int>& a,
		vector<vector<std::pair<int, int> > > segmentsWithPartitions) {
	int phases = std::ceil(
			std::log(segmentsWithPartitions.size()) / std::log(2));

	accuSavedTime = 0;
	for (int phaseIndex = 0; phaseIndex < phases; phaseIndex++) {

		if (!parallel || segmentsWithPartitions.size() == 2
				|| fineGrained == true) {
			for (int internalMergeTaskIndex = 0;
					internalMergeTaskIndex < segmentsWithPartitions.size();
					internalMergeTaskIndex++) {

				if (internalMergeTaskIndex == segmentsWithPartitions.size() - 1)
					break;

				vector<std::pair<int, int> > twoLists;

				for (int i = 0;
						i
								< segmentsWithPartitions[internalMergeTaskIndex].size();
						i++)
					twoLists.push_back(
							segmentsWithPartitions[internalMergeTaskIndex][i]);
				for (int i = 0;
						i
								< segmentsWithPartitions[internalMergeTaskIndex
										+ 1].size(); i++)
					twoLists.push_back(
							segmentsWithPartitions[internalMergeTaskIndex + 1][i]);
				segmentsWithPartitions[internalMergeTaskIndex] = twoLists;
				segmentsWithPartitions.erase(
						segmentsWithPartitions.begin() + internalMergeTaskIndex
								+ 1);
				shuffleRanges(a,
						segmentsWithPartitions[internalMergeTaskIndex]);
			}
		} else {
			vector<vector<std::pair<int, int> > > temp;

			for (int internalMergeTaskIndex = 0;
					internalMergeTaskIndex < (segmentsWithPartitions.size() / 2);
					internalMergeTaskIndex++) {
				vector<std::pair<int, int> > tempVector;
				temp.push_back(tempVector);
			}

			struct timeval startwtime1, endwtime1;
			double seq_time1;
			gettimeofday(&startwtime1, NULL);
			double maxThread = 0;

			//	cout << "phase " << phaseIndex << " used thread here is:" << (segmentsWithPartitions.size() / 2) <<endl;
#pragma omp parallel num_threads((segmentsWithPartitions.size() / 2))
//            for (int internalMergeTaskIndex = 0; internalMergeTaskIndex < (segmentsWithPartitions.size() / 2); internalMergeTaskIndex++)
			{
				struct timeval startwtime2, endwtime2;
				double seq_time2;
				gettimeofday(&startwtime2, NULL);
				int internalMergeTaskIndex = omp_get_thread_num();

				vector<std::pair<int, int> > twoLists;

				for (int i = 0;
						i
								< segmentsWithPartitions[2
										* internalMergeTaskIndex].size(); i++)
					twoLists.push_back(
							segmentsWithPartitions[2 * internalMergeTaskIndex][i]);
				for (int i = 0;
						i
								< segmentsWithPartitions[2
										* internalMergeTaskIndex + 1].size();
						i++)
					twoLists.push_back(
							segmentsWithPartitions[2 * internalMergeTaskIndex
									+ 1][i]);
				segmentsWithPartitions[2 * internalMergeTaskIndex] = twoLists;
				shuffleRanges(a,
						segmentsWithPartitions[2 * internalMergeTaskIndex]);
				temp[internalMergeTaskIndex] = segmentsWithPartitions[2
						* internalMergeTaskIndex];
				gettimeofday(&endwtime2, NULL);
				seq_time2 = (double) ((endwtime1.tv_usec - startwtime1.tv_usec)
						/ 1.0e6 + endwtime1.tv_sec - startwtime1.tv_sec);
				if (seq_time2 > maxThread)
					maxThread = seq_time2;
			}

			gettimeofday(&endwtime1, NULL);
			seq_time1 = (double) ((endwtime1.tv_usec - startwtime1.tv_usec)
					/ 1.0e6 + endwtime1.tv_sec - startwtime1.tv_sec);
			accuSavedTime += seq_time1 - maxThread;

			if ((segmentsWithPartitions.size() % 2) == 1)
				temp.push_back(
						segmentsWithPartitions[segmentsWithPartitions.size() - 1]);

			segmentsWithPartitions.clear();
			for (int addIndex = 0; addIndex < temp.size(); addIndex++) {
				segmentsWithPartitions.push_back(temp[addIndex]);
			}
		}
	}
	return true;
}

bool createPartitionsFromSegments(vector<int> a, int size, int numOfSegments,
		int numOfPartitions) {

	cout << "\nPS-merge array size: " << size << ", # of segments: "
			<< numOfSegments << ", # of threads: " << numOfPartitions << endl;

	vector<int> segmentsEnds;
	cout << "Generating data is started...." << flush;
	setSegmentsEndEqually(size, numOfSegments, segmentsEnds);
	generateArray(a, size, numOfSegments, segmentsEnds);
	cout << "Done!" << endl;

	struct timeval startwtime1, endwtime1;
	gettimeofday(&startwtime1, NULL);

	//actual index + 1, if it contains index less than segment start then
	//this segment in not included
	vector<std::pair<int, vector<int> > > boundaries;

	struct timeval start1, end1;
	gettimeofday(&start1, NULL);

	mergeComparsons = 0;
	cout << "Partitioning data is started...." << flush;
	FrancisPartitioning(a, numOfSegments, numOfPartitions, segmentsEnds,
			boundaries);
	gettimeofday(&end1, NULL);
	cout << "Done!" << endl;

	partitioningTime = (double) ((end1.tv_usec - start1.tv_usec) / 1.0e6
			+ end1.tv_sec - start1.tv_sec);

	cout << "Merging (tournament step) is started..." << flush;
	vector<vector<std::pair<int, int> > > segmentsWithPartitions;
	int lowerBound, upperBound;

	for (int segmentIndex = 0; segmentIndex < numOfSegments; segmentIndex++) {
		vector<std::pair<int, int> > segmentRanges;
		segmentsWithPartitions.push_back(segmentRanges);
	}

	vector<int> currentSegmentStarts;

	//set segments starts
	currentSegmentStarts.push_back(0);
	for (int i = 0; i < segmentsEnds.size() - 1; i++) {
		currentSegmentStarts.push_back(*(segmentsEnds.begin() + i) + 1);
	}
	vector<partitionInfo> partitionsInfo;

	for (int i = 0; i < boundaries.size(); i++) {
		//initialize the semgmentinfo
		int tempSize = 0;
		for (int b = 0; b < (boundaries.begin() + i)->second.size(); b++) {
			if (*((boundaries.begin() + i)->second.begin() + b)
					>= *(currentSegmentStarts.begin() + b)) {
				tempSize++;
			}
		}
		partitionInfo s1;
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
		partitionsInfo.push_back(s1);
	}

	//convert each segment into sub-segments according to
	//the partitions, and filling dummy ranges (-1,-1) for partition with no share
	//at this segment.
	for (int partitionIndex = 0; partitionIndex < partitionsInfo.size();
			partitionIndex++) {
		for (int segmentIndex = 0, subSegmentIndex = 0;
				segmentIndex < numOfSegments; segmentIndex++) {

			if (segmentIndex == 0)
				lowerBound = 0;
			else
				lowerBound = segmentsEnds[segmentIndex - 1] + 1;
			upperBound = segmentsEnds[segmentIndex];
			if (subSegmentIndex < partitionsInfo[partitionIndex].noOfSegments
					&& partitionsInfo[partitionIndex].second[subSegmentIndex]
							>= lowerBound
					&& partitionsInfo[partitionIndex].second[subSegmentIndex]
							<= upperBound) {
				segmentsWithPartitions[segmentIndex].push_back(
						std::make_pair(
								partitionsInfo[partitionIndex].first[subSegmentIndex],
								partitionsInfo[partitionIndex].second[subSegmentIndex]));
				subSegmentIndex++;
			} else {
				segmentsWithPartitions[segmentIndex].push_back(
						std::make_pair(-1, -1));
			}
		}
	}

	if (!kwayMerge(a, segmentsWithPartitions)) {
		gettimeofday(&endwtime1, NULL);
		elapsed_time = (double) ((endwtime1.tv_usec - startwtime1.tv_usec)
				/ 1.0e6 + endwtime1.tv_sec - startwtime1.tv_sec);

		cout << "SplitMerge size:" << size << " # of seg: " << numOfSegments
				<< " # of threads: " << numOfPartitions << " wall clock time = "
				<< elapsed_time << endl;
		cout.flush();
		return false;
	}
	gettimeofday(&endwtime1, NULL);
	elapsed_time = (double) ((endwtime1.tv_usec - startwtime1.tv_usec) / 1.0e6
			+ endwtime1.tv_sec - startwtime1.tv_sec);

	cout << "Done!" << endl;

	return true;
}

bool ps_merge(int size, int numOfSegments, int numberOfthreads) {
	vector<int> a;

	resetNumOfComparisonsSplitMerge();
	int numOfPartitions = numberOfthreads;

	mergeTime = partTime = shuffleTime = 0;
	a.clear();
	numOfThreads = numOfPartitions;
	numOfSwaps = 0;
	mergeSwaps = 0;
	shuffleSwaps = 0;

	if (!createPartitionsFromSegments(a, size, numOfSegments,
			numOfPartitions)) {
		cout << flush;
		cout << "shuffle time = " << shuffleTime << endl;
		return false;
	}

	if (checkIt(a))
		cout << "The merged array is checked successfully!" << endl;
	return true;
}

/**
 * 
 * @return 
 */
int main() {

	//k value in the k-way merge
	int array_size = 50000000; //std::pow(2.0, 24);
	int numberOfthreads = 4;
	int noOfSegments = 128;
	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(numberOfthreads);

	ps_merge(array_size, noOfSegments, numberOfthreads);

#if IS_DEBUG == 0
	cout << "Partitioning time = " << partitioningTime
	<< ", shuffle time = " << shuffleTime << ", ";
	cout << " merge time = " << mergeTime << " seconds."<< endl;
	cout << "The wall clock total merging time = " << elapsed_time
	<< " seconds." << endl;
#endif

#if IS_DEBUG == 1
	numOfSwaps = mergeSwaps + shuffleSwaps;
	std::cout.precision(0);
	cout << std::fixed;
	cout << "All swaps = " << numOfSwaps << " , # of shuffle swaps = "
			<< shuffleSwaps << " , # of merge swaps = " << mergeSwaps << endl;
#endif

	cout << "\nAll tests are done!" << endl;
	return 0;
}
