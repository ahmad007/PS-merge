#include "splitMergeOnDistance.h"
#include <iostream>
#include <omp.h>
#include <cmath>
#include <map>
#include <stack>
#include <climits>
using namespace std;

int shiftTime1 = 0;
unsigned long noOfSwaps1 = 0;
int counter1 = 0;
vector<IndexCache> cache1;
vector<IndexCache> cache2;
long *threads1;

struct param {
public:

	param(int p1, int p2, int p3, int p4) {
		f1 = p1;
		f2 = p2;
		last = p3;
		portionOneLastIndex = p4;
	}
	int f1;
	int f2;
	int last;
	int portionOneLastIndex;
};

long getNumOfSwapsSplitMergeOnDistance() {
	if (counter1 > 0)
		cout << "********swaps multiple*******" << counter1 << endl;
	return noOfSwaps1;
}

/* This method convert virtual indexes of loop to real indexes it convert ranges rather than element indexes
 * @param start loop virtual start index
 * @param end loop virtual end index
 * @param segmentInfoData partition informaiton, virtual and corrsponding real ranges
 * @param iVector list of real range of loop i variable
 * @param jVector list of real range of loop j variable
 *
 */

void getLoops(int start, int end, segmentInfo segmentInfoData,
		vector<std::pair<int, int> >& iVector,
		vector<std::pair<int, int> >& jVector) {

	vector<std::pair<int, int> > realLoop;

	IndexCache startCache;
	IndexCache endCache;
	startCache.setSize(1);
	endCache.setSize(1);
	int i = 0, realStart, realEnd;

	while (start <= end) {
		segmentInfoData.convertVirtualIndextoActual(start, startCache);
		segmentInfoData.convertVirtualIndextoActual(end, endCache);
		realStart = startCache.actualStart[0]
				+ (start - startCache.virtualStart[0]);
		realEnd = endCache.actualStart[0] + (end - endCache.virtualStart[0]);

		if (startCache.actualStart[0] < endCache.actualStart[0]) {
			realLoop.insert(realLoop.begin() + i,
					std::make_pair(realStart, startCache.actualEnd[0]));
			realLoop.insert(realLoop.begin() + i + 1,
					std::make_pair(endCache.actualStart[0], realEnd));
			i++;

		} else if (startCache.actualStart[0] == endCache.actualStart[0]) {
			realLoop.insert(realLoop.begin() + i,
					std::make_pair(realStart, realEnd));
		}
		start = startCache.virtualEnd[0] + 1;
		end = endCache.virtualStart[0] - 1;
	}

	startCache.finalize();
	endCache.finalize();

	typedef std::vector<std::pair<int, int> > my_vector;
	my_vector::iterator it = realLoop.begin();
	int iSize = 0, jSize = 0;

	for (int i = 0, j = realLoop.size() - 1; i <= j;) {
		if (i == j) {
			iVector.push_back(*(it + i));
			jVector.push_back(*(it + i));
			break;
		}

		iSize = (it + i)->second - (it + i)->first + 1;
		jSize = (it + j)->second - (it + j)->first + 1;

		if (iSize < jSize) {
			iVector.push_back(*(it + i));
			jVector.push_back(
					std::make_pair((it + j)->second - iSize + 1,
							(it + j)->second));
			(it + j)->second -= iSize;
			i++;
		} else if (iSize > jSize) {
			iVector.push_back(
					std::make_pair((it + i)->first,
							(it + i)->first + jSize - 1));
			(it + i)->first += jSize;
			jVector.push_back(*(it + j));
			j--;
		} else if (iSize == jSize) {
			iVector.push_back(*(it + i));
			jVector.push_back(*(it + j));
			i++;
			j--;
		}

	}
}

/* This method performs the loop optimized rotation of items
 * @param segmentInfoData partition information; sub-segments ranges
 * @param arr input list of k-segment
 * @param  shiftStartIndex virtual index of shifting start index
 * @param  shiftEndIndex virtual index of shifting end index
 * @param  shiftValue shift value
 * @param  portionOneLastIndex virtual end index of sub-segment 1
 */

void rightShiftCircularOnDistanceOptimized(segmentInfo segmentInfoData,
		vector<int>& arr, int shiftStartIndex, int shiftEndIndex,
		int shiftValue, int portionOneLastIndex) {

	int temp = 0;
	int threadid = omp_get_thread_num();
	vector<std::pair<int, int> > vectorI;
	vector<std::pair<int, int> > vectorJ;

	vectorI.clear();
	vectorJ.clear();
	getLoops(shiftStartIndex, shiftEndIndex - 1, segmentInfoData, vectorI,
			vectorJ);

	for (int iRange = 0; iRange < vectorI.size(); iRange++) {

		for (int i = (vectorI.begin() + iRange)->first, j = (vectorJ.begin()
				+ iRange)->second + 1;
				--j > i && i <= (vectorI.begin() + iRange)->second; i++) {
			temp = arr[i];
			arr[i] = arr[j];
			arr[j] = temp;
#if IS_DEBUG == 1
#pragma omp critical
			{
				if (noOfSwaps1 == ULONG_MAX) {
					counter1++;
					noOfSwaps1 = 1;
				} else
				noOfSwaps1++;
			}
#endif
		}
	}

	vectorI.clear();
	vectorJ.clear();
	getLoops(shiftStartIndex, shiftStartIndex + shiftValue - 1, segmentInfoData,
			vectorI, vectorJ);

	for (int iRange = 0; iRange < vectorI.size(); iRange++) {

		for (int i = (vectorI.begin() + iRange)->first, j = (vectorJ.begin()
				+ iRange)->second + 1;
				--j > i && i <= (vectorI.begin() + iRange)->second; i++) {

			temp = arr[i];
			arr[i] = arr[j];
			arr[j] = temp;
#if IS_DEBUG == 1
#pragma omp critical
			{
				if (noOfSwaps1 == ULONG_MAX) {
					counter1++;
					noOfSwaps1 = 1;
				} else
				noOfSwaps1++;
			}
#endif
		}
	}

	vectorI.clear();
	vectorJ.clear();
	getLoops(shiftStartIndex + shiftValue, shiftEndIndex - 1, segmentInfoData,
			vectorI, vectorJ);

	for (int iRange = 0; iRange < vectorI.size(); iRange++) {
		for (int i = (vectorI.begin() + iRange)->first, j = (vectorJ.begin()
				+ iRange)->second + 1;
				--j > i && i <= (vectorI.begin() + iRange)->second; i++) {
			temp = arr[i];
			arr[i] = arr[j];
			arr[j] = temp;
#if IS_DEBUG == 1
#pragma omp critical
			{
				if (noOfSwaps1 == ULONG_MAX) {
					counter1++;
					noOfSwaps1 = 1;
				} else
				noOfSwaps1++;
			}
#endif
		}
	}
}

/* This method performs the Lazy-Merge  local 2-way merging routine for two sub-segments of this partition
 * @param segmentInfoData partition information; sub-segments ranges
 * @param aa input list of k-segment
 * @param  f1 virtual start index of sub-segment 1
 * @param  f2 virtual start index of sub-segment 2
 * @param  last virtual end index of sub-segment 2
 * @param  portionOneLastIndex virtual end index of sub-segment 1
 */

void iterativeSplitMergeOnDistanceCache(segmentInfo segmentInfoData,
		vector<int>& a, int f1, int f2, int last, int portionOneLastIndex) {
	std::stack<param> parameters;

	param p(f1, f2, last, portionOneLastIndex);
	parameters.push(p);

	int l, r, ldash, rdash, m, mdash;

	while (!parameters.empty()) {

		param temp = parameters.top();
		parameters.pop();

		f1 = temp.f1;
		f2 = temp.f2;
		last = temp.last;
		portionOneLastIndex = temp.portionOneLastIndex;

		if (f1 >= f2 || f2 >= last || f1 == -1 || f2 == -1) {
			continue;
		}
		int l = f1, r = f2, ldash = f2, rdash = last, m = 0, mdash = 0;

		while (!(l >= r && ldash >= rdash)) {
			if (l < r) {
				m = (l + r) / 2;
			}
			if (ldash < rdash) {
				mdash = (ldash + rdash) / 2;
			}
			if (a[segmentInfoData.convertVirtualIndextoActual(m)]
					<= a[segmentInfoData.convertVirtualIndextoActual(mdash)]) {
				l = m + 1;
				rdash = mdash;
			} else {
				ldash = mdash + 1;
				r = m;
			}
		}

		rightShiftCircularOnDistanceOptimized(segmentInfoData, a, r, ldash,
				(ldash - f2), portionOneLastIndex);
		param pFirst(f1, r, r + rdash - f2, portionOneLastIndex);
		param pSecond(l + ldash - f2, ldash, last, portionOneLastIndex);
		parameters.push(pSecond);
		parameters.push(pFirst);

	}
}

void SplitMergeOnDistanceCache(segmentInfo segmentInfoData, vector<int>& a,
		int f1, int f2, int last, int portionOneLastIndex) {

	if (f1 >= f2 || f2 >= last) {
		return;
	}
	int l = f1, r = f2, ldash = f2, rdash = last, m = 0, mdash = 0;

	while (!(l >= r && ldash >= rdash)) {
		if (l < r) {
			m = (l + r) / 2;
		}
		if (ldash < rdash) {
			mdash = (ldash + rdash) / 2;
		}
		if (a[segmentInfoData.convertVirtualIndextoActual(m)]
				<= a[segmentInfoData.convertVirtualIndextoActual(mdash)]) {
			l = m + 1;
			rdash = mdash;
		} else {
			ldash = mdash + 1;
			r = m;
		}
	}

	//rightShiftCircularOnDistanceCache(segmentInfoData, a, r, ldash, (ldash - f2), portionOneLastIndex);
	rightShiftCircularOnDistanceOptimized(segmentInfoData, a, r, ldash,
			(ldash - f2), portionOneLastIndex);
	SplitMergeOnDistanceCache(segmentInfoData, a, f1, r, r + rdash - f2,
			portionOneLastIndex);
	SplitMergeOnDistanceCache(segmentInfoData, a, l + ldash - f2, ldash, last,
			portionOneLastIndex);
}

void initializeAndStartCache(vector<segmentInfo> segmentsInfos, vector<int>& a,
		bool dir, int i) {

	int threadID = i; //omp_get_thread_num();
	segmentInfo segmentInfoData = *(segmentsInfos.begin() + threadID);

	vector<std::pair<int, int> > tempBoundaries;
	for (int i = 0; i < segmentInfoData.noOfSegments; i++)
		tempBoundaries.push_back(
				std::make_pair(segmentInfoData.first[i],
						segmentInfoData.second[i]));

	//#pragma omp critical
	//    cout << threadID << " ok"
	//            << "first: " << segmentInfoData.boundries.begin()->first
	//            << " size: " << segmentInfoData.size << endl;

	int tempStart, tempEnd;
	typedef std::vector<std::pair<int, int> > my_vector;
	my_vector::iterator it;

	unsigned virtualStart1, virtualStart2, virtualEnd1, virtualEnd2, distance;
	unsigned currentSize = 0;
	it = tempBoundaries.begin();
	currentSize = it->second - it->first + 1;

	while (tempBoundaries.size() > 1) {
		it = tempBoundaries.begin();

		virtualStart1 = it->first;
		virtualStart2 = currentSize + virtualStart1;
		virtualEnd1 = it->first + currentSize - 1;
		virtualEnd2 = ((it + 1)->second - (it + 1)->first + 1) + virtualStart2;

		//adding the new merged segment size
		currentSize += (it + 1)->second - (it + 1)->first + 1;

		tempStart = it->first;
		tempEnd = tempStart + currentSize - 1;

		//array, virtual f1, virtual f2, size +1, virtual end1, distance
		//SplitMergeOnDistanceCache(segmentInfoData, a, virtualStart1, virtualStart2, virtualEnd2, virtualEnd1);

		iterativeSplitMergeOnDistanceCache(segmentInfoData, a, virtualStart1,
				virtualStart2, virtualEnd2, virtualEnd1);
		tempBoundaries.erase(it, it + 2);
		tempBoundaries.insert(it, std::make_pair(tempStart, tempEnd));
	}
}

int getNumOfPhases1(int x) {

	int k = 1;
	while (k < x)
		k = k << 1;

	int counter = 0;
	while (k > 0) {
		counter++;
		k /= 2;
	}
	return counter - 1;
}

/* This method performs the Lazy-Merge  local 2-way merging routine for all sub-segment of this partition
 * @param segmentsInfos list of each partition information; sub-segments ranges
 * @param arr input list of k-segment
 * @param  bool direction of sorting
 */

void initializeAndStartCacheTree(vector<segmentInfo> segmentsInfos,
		vector<int>& a, bool dir) {
	int threadID = omp_get_thread_num();
	segmentInfo segmentInfoData = *(segmentsInfos.begin() + threadID);

	vector<std::pair<int, int> > tempBoundaries;
	unsigned virtualStart1, virtualStart2, virtualEnd1, virtualEnd2, distance;
	int numOfPhases = getNumOfPhases1(segmentInfoData.noOfSegments);
	unsigned tempStart, tempEnd;
	int size1, size2, accumalativeSize = 0;

	for (int i = 0; i < segmentInfoData.noOfSegments; i++) {
		tempBoundaries.push_back(
				std::make_pair(segmentInfoData.first[i],
						segmentInfoData.second[i]));
	}

	for (int phaseIndex = 0; phaseIndex < numOfPhases; phaseIndex++) {
		accumalativeSize = tempBoundaries[0].first;

		// accumalativeSize = 0;
		//merging at tree level equals phaseIndex
		for (int mergingIndex = 0; mergingIndex < tempBoundaries.size();
				mergingIndex += 2) {
			if (mergingIndex == tempBoundaries.size() - 1)
				break;
			size1 = tempBoundaries[mergingIndex].second
					- tempBoundaries[mergingIndex].first + 1;
			size2 = tempBoundaries[mergingIndex + 1].second
					- tempBoundaries[mergingIndex + 1].first + 1;
			virtualStart1 = accumalativeSize;
			virtualEnd1 = virtualStart1 + size1 - 1;
			virtualStart2 = virtualEnd1 + 1;
			virtualEnd2 = virtualStart2 + size2;
			iterativeSplitMergeOnDistanceCache(segmentInfoData, a,
					virtualStart1, virtualStart2, virtualEnd2, virtualEnd1);
			accumalativeSize += size1 + size2;
		}

		//reduce the segments by half
		//accumalativeSize = 0;
		accumalativeSize = tempBoundaries[0].first;

		for (int i = 0; i < tempBoundaries.size(); i++) {
			if (i == tempBoundaries.size() - 1)
				break;

			tempStart = accumalativeSize;
			size1 = tempBoundaries[i].second - tempBoundaries[i].first + 1;
			size2 = tempBoundaries[i + 1].second - tempBoundaries[i + 1].first
					+ 1;

			tempEnd = tempStart + size1 + size2 - 1;
			accumalativeSize += size1 + size2;
			tempBoundaries.erase(tempBoundaries.begin() + i,
					tempBoundaries.begin() + i + 2);
			tempBoundaries.insert(tempBoundaries.begin() + i,
					std::make_pair(tempStart, tempEnd));

		}

	}
}

/* This method performs the Lazy-Merge tournamennt tree using SplitMerge as local 2-way merging routine
 * @param arr input list of k-segment
 * @param segmentsInfos list of each partition information
 */
void multiwaySplitMergeOnDistanceCache(vector<int>& arr,
		vector<segmentInfo> segmentsInfos, long& threadTime) {
	noOfSwaps1 = 0;
	int noOfThreads = segmentsInfos.size();
	cache1.clear();
	cache2.clear();

	threads1 = new long[noOfThreads];

	for (int i = 0; i < noOfThreads; i++) {
		IndexCache c1;
		c1.setSize(1);
		cache1.push_back(c1);
		IndexCache c2;
		c2.setSize(1);
		cache2.push_back(c2);
		threads1[i] = 0;
	}

#pragma omp parallel num_threads(noOfThreads)
	{
		initializeAndStartCacheTree(segmentsInfos, arr, true);
	}

	for (int i = 0; i < noOfThreads; i++) {
		cache1[i].finalize();
		cache2[i].finalize();
	}

	threadTime = -1;
	for (int i = 0; i < noOfThreads; i++) {
		if (threads1[i] > threadTime)
			threadTime = threads1[i];
	}
	delete[] threads1;
}

map<int, int> getMap(segmentInfo & segmentInfoData) {
	map<int, int> map;
	int virtualIndex = segmentInfoData.first[0] - 1;
	//map[segmentInfoData.first[0]] = segmentInfoData.first[0];
	int tempSize = 0;
	for (int i = 0; i < segmentInfoData.noOfSegments; i++) {
		segmentInfoData.virtualStart[i] = virtualIndex + 1;
		tempSize += segmentInfoData.second[i] - segmentInfoData.first[i] + 1;
		virtualIndex = segmentInfoData.first[0] + tempSize - 1;
		segmentInfoData.virtualEnd[i] = virtualIndex;
		map[virtualIndex] = i;
	}
	int threadID = omp_get_thread_num();
	return map;
}

void rightShiftCircularOnDistanceMap(segmentInfo segmentInfoData,
		vector<int>& arr, int shiftStartIndex, int shiftEndIndex,
		int shiftValue, int portionOneLastIndex, map<int, int> map) {

	int threadid = omp_get_thread_num();
	int segmentIndex = -1, actualI = -1, actualJ = -1;

	int temp = 0;
	int tempi, tempj;
	for (int i = shiftStartIndex, j = shiftEndIndex; --j > i; i++) {
		segmentIndex = map.lower_bound(i)->second;
		actualI = segmentInfoData.first[segmentIndex] + i
				- segmentInfoData.virtualStart[segmentIndex];

		segmentIndex = map.lower_bound(j)->second;
		actualJ = segmentInfoData.first[segmentIndex] + j
				- segmentInfoData.virtualStart[segmentIndex];

		//#pragma omp critical
		//        cout << "thread: " << threadid << " obtained: " << tempi << " virIndex: " << i <<
		//                " virStart: " << cache1[threadid].virtualStart[0] <<
		//                " virEnd: " << cache1[threadid].virtualEnd[0] <<
		//                " actStart: " << cache1[threadid].actualStart [0] <<
		//                " actEnd: " << cache1[threadid].actualEnd[0] <<
		//                " first[0]: " << segmentInfoData.first[0] << endl;
		//
		//#pragma omp critical
		//        cout << "thread: " << threadid << " obtained: " << tempj << " virIndex: " << j <<
		//                " virStart: " << cache2[threadid].virtualStart[0] <<
		//                " virEnd: " << cache2[threadid].virtualEnd[0] <<
		//                " actStart: " << cache2[threadid].actualStart [0] <<
		//                " actEnd: " << cache2[threadid].actualEnd[0] <<
		//                " first[0]: " << segmentInfoData.first[0] << endl;

		//        if (tempi >= arr.size()) {
		//            tempi = cache1[threadid].actualStart[0] + (i - cache1[threadid].virtualStart[0]);
		//            cout << "&&&&&&&&&&&&&&&&&&&& tempi=" << tempi << ", i=" << i << endl;
		//            for (int k = 0; k < segmentInfoData.noOfSegments; k++) {
		//                cout << "(" << segmentInfoData.first[k] <<
		//                        "," << segmentInfoData.second[k] << ") ";
		//            }
		//            cout << endl;
		//        }
		//
		//        if (tempj >= arr.size()) {
		//            tempj = cache2[threadid].actualStart[0] + (j - cache2[threadid].virtualStart[0]);
		//            cout << "&&&&&&&&&&&&&&&&&&&& tempj=" << tempj << ", j=" << j << endl;
		//            for (int k = 0; k < segmentInfoData.noOfSegments; k++) {
		//                cout << "(" << segmentInfoData.first[k] <<
		//                        "," << segmentInfoData.second[k] << ") ";
		//            }
		//            cout << endl;
		//        }

		temp = arr[actualI];
		arr[actualI] = arr[actualJ];
		arr[actualJ] = temp;
//#pragma omp critical
//        noOfSwaps1++;
	}

	for (int i = shiftStartIndex, j = i + shiftValue; --j > i; i++) {
		segmentIndex = map.lower_bound(i)->second;
		actualI = segmentInfoData.first[segmentIndex] + i
				- segmentInfoData.virtualStart[segmentIndex];

		segmentIndex = map.lower_bound(j)->second;
		actualJ = segmentInfoData.first[segmentIndex] + j
				- segmentInfoData.virtualStart[segmentIndex];

		temp = arr[actualI];
		arr[actualI] = arr[actualJ];
		arr[actualJ] = temp;
//#pragma omp critical
//        noOfSwaps1++;
	}

	for (int i = shiftStartIndex + shiftValue, j = shiftEndIndex; --j > i;
			i++) {
		segmentIndex = map.lower_bound(i)->second;
		actualI = segmentInfoData.first[segmentIndex] + i
				- segmentInfoData.virtualStart[segmentIndex];

		segmentIndex = map.lower_bound(j)->second;
		actualJ = segmentInfoData.first[segmentIndex] + j
				- segmentInfoData.virtualStart[segmentIndex];

		temp = arr[actualI];
		arr[actualI] = arr[actualJ];
		arr[actualJ] = temp;
//#pragma omp critical
//        noOfSwaps1++;
	}
}

void SplitMergeOnDistanceMap(segmentInfo segmentInfoData, vector<int>& a,
		int f1, int f2, int last, int portionOneLastIndex, map<int, int> map) {

	int startIndex = segmentInfoData.first[0];
	int actualM = -1, actualMDash = -1, segmentIndex;

	if (f1 >= f2 || f2 >= last) {
		return;
	}
	int l = f1, r = f2, ldash = f2, rdash = last, m = 0, mdash = 0;

	while (!(l >= r && ldash >= rdash)) {
		if (l < r) {
			m = (l + r) / 2;
		}
		if (ldash < rdash) {
			mdash = (ldash + rdash) / 2;
		}

		int threadID = omp_get_thread_num();
		segmentIndex = map.lower_bound(m)->second;
		actualM = segmentInfoData.first[segmentIndex] + m
				- segmentInfoData.virtualStart[segmentIndex];

		segmentIndex = map.lower_bound(mdash)->second;
		actualMDash = segmentInfoData.first[segmentIndex] + mdash
				- segmentInfoData.virtualStart[segmentIndex];

		if (a[actualM] <= a[actualMDash]) {
			l = m + 1;
			rdash = mdash;
		} else {
			ldash = mdash + 1;
			r = m;
		}
	}

	rightShiftCircularOnDistanceMap(segmentInfoData, a, r, ldash, (ldash - f2),
			portionOneLastIndex, map);
	SplitMergeOnDistanceMap(segmentInfoData, a, f1, r, r + rdash - f2,
			portionOneLastIndex, map);
	SplitMergeOnDistanceMap(segmentInfoData, a, l + ldash - f2, ldash, last,
			portionOneLastIndex, map);
}

void initializeAndStartMap(vector<segmentInfo> segmentsInfos, vector<int>& a,
		bool dir) {

	int threadID = omp_get_thread_num();
	segmentInfo segmentInfoData = *(segmentsInfos.begin() + threadID);

	vector<std::pair<int, int> > tempBoundaries;
	for (int i = 0; i < segmentInfoData.noOfSegments; i++)
		tempBoundaries.push_back(
				std::make_pair(segmentInfoData.first[i],
						segmentInfoData.second[i]));

	//#pragma omp critical
	//    cout << threadID << " ok"
	//            << "first: " << segmentInfoData.boundries.begin()->first
	//            << " size: " << segmentInfoData.size << endl;

	int tempStart, tempEnd;
	typedef std::vector<std::pair<int, int> > my_vector;
	my_vector::iterator it;

	unsigned virtualStart1, virtualStart2, virtualEnd1, virtualEnd2, distance;
	unsigned currentSize = 0;
	it = tempBoundaries.begin();
	currentSize = it->second - it->first + 1;

	map<int, int> map = getMap(segmentInfoData);

	//        if (threadID == 0) {
	//            cout << "map info:\n";
	//            std::map<int, int>::iterator it1;
	//            it1 = map.begin();
	//
	//            for (; it1 != map.end(); it1++)
	//                cout << it1->first << "," << it1->second << endl;
	//        }
	//    return;

	while (tempBoundaries.size() > 1) {
		it = tempBoundaries.begin();

		virtualStart1 = it->first;
		virtualStart2 = currentSize + virtualStart1;
		virtualEnd1 = it->first + currentSize - 1;
		virtualEnd2 = ((it + 1)->second - (it + 1)->first + 1) + virtualStart2;

		//adding the new merged segment size
		currentSize += (it + 1)->second - (it + 1)->first + 1;

		tempStart = it->first;
		tempEnd = tempStart + currentSize - 1;

		//array, virtual f1, virtual f2, size +1, virtual end1, distance
		SplitMergeOnDistanceMap(segmentInfoData, a, virtualStart1,
				virtualStart2, virtualEnd2, virtualEnd1, map);
		tempBoundaries.erase(it, it + 2);
		tempBoundaries.insert(it, std::make_pair(tempStart, tempEnd));
	}
}

void multiwaySplitMergeOnDistanceMap(vector<int>& arr,
		vector<segmentInfo> segmentsInfos) {

	cout << "\n\nLazy merge based on split merging started...." << endl;
	const clock_t begin_time = clock();
	noOfSwaps1 = 0;
	int noOfThreads = segmentsInfos.size();

#pragma omp parallel num_threads(noOfThreads)
	{
		initializeAndStartMap(segmentsInfos, arr, true);
	}

	cout << "Lazy merge based on split merging ended." << endl;
//    std::cout << "Elapsed time:  " << float( clock() - begin_time) / CLOCKS_PER_SEC << " seconds\n";
//    std::cout << "Elapsed time:  " << float( clock() - begin_time) / 1000l << " milli-seconds\n\n";
}

bool checkArray(vector<int>& indexes, vector<int>& a,
		vector<segmentInfo> segmentsInfos) {

	int lastElement, firstElement;
	for (int i = 0; i < segmentsInfos.size(); i++) {
		if (i > 0) {
			lastElement =
					segmentsInfos[i - 1].second[segmentsInfos[i - 1].noOfSegments
							- 1];
			firstElement = segmentsInfos[i].first[0];
			if (a[lastElement] > a[firstElement]) {
				cout << "error at: a[" << lastElement << "]=" << a[lastElement]
						<< "[" << firstElement << "]=" << a[firstElement]
						<< endl;
				return false;
			}
		}

		for (int k = 0, index = 0; k < segmentsInfos[i].noOfSegments; k++) {
			if (k > 0) {

				if (a[segmentsInfos[i].second[k - 1]]
						> a[segmentsInfos[i].first[k]]) {
					cout << "error at: a[" << segmentsInfos[i].second[k - 1]
							<< "]=" << a[segmentsInfos[i].second[k - 1]] << "["
							<< segmentsInfos[i].first[k] << "]="
							<< a[segmentsInfos[i].first[k]] << endl;
					return false;
				}
			}
			for (int r = segmentsInfos[i].first[k];
					r <= segmentsInfos[i].second[k] - 1; r++) {
				//	indexes[index++] =r;
				if (a[r] > a[r + 1]) {
					cout << "error at: a[" << r << "]=" << a[r] << "[" << r + 1
							<< "]=" << a[r + 1] << endl;
					return false;
				}
			}
		}
	}
	cout << "check is done successfully!\n";
	return true;
}

void rightShiftCircularOnDistance(segmentInfo segmentInfoData, vector<int>& arr,
		int shiftStartIndex, int shiftEndIndex, int shiftValue,
		int portionOneLastIndex) {

	int temp = 0;
	int tempi, tempj;
	for (int i = shiftStartIndex, j = shiftEndIndex; --j > i; i++) {
		tempi = segmentInfoData.convertVirtualIndextoActual(i);
		tempj = segmentInfoData.convertVirtualIndextoActual(j);

		temp = arr[tempi];
		arr[tempi] = arr[tempj];
		arr[tempj] = temp;
//#pragma omp critical
//        noOfSwaps1++;
	}

	for (int i = shiftStartIndex, j = i + shiftValue; --j > i; i++) {
		tempi = segmentInfoData.convertVirtualIndextoActual(i);
		tempj = segmentInfoData.convertVirtualIndextoActual(j);
		temp = arr[tempi];
		arr[tempi] = arr[tempj];
		arr[tempj] = temp;
//#pragma omp critical
//        noOfSwaps1++;
	}

	for (int i = shiftStartIndex + shiftValue, j = shiftEndIndex; --j > i;
			i++) {
		tempi = segmentInfoData.convertVirtualIndextoActual(i);
		tempj = segmentInfoData.convertVirtualIndextoActual(j);
		temp = arr[tempi];
		arr[tempi] = arr[tempj];
		arr[tempj] = temp;
//#pragma omp critical
//        noOfSwaps1++;
	}
}

void SplitMergeOnDistance(segmentInfo segmentInfoData, vector<int>& a, int f1,
		int f2, int last, int portionOneLastIndex) {

	if (f1 >= f2 || f2 >= last) {
		return;
	}
	int l = f1, r = f2, ldash = f2, rdash = last, m = 0, mdash = 0;

	while (!(l >= r && ldash >= rdash)) {
		if (l < r) {
			m = (l + r) / 2;
		}
		if (ldash < rdash) {
			mdash = (ldash + rdash) / 2;
		}
		if (a[segmentInfoData.convertVirtualIndextoActual(m)]
				<= a[segmentInfoData.convertVirtualIndextoActual(mdash)]) {
			l = m + 1;
			rdash = mdash;
		} else {
			ldash = mdash + 1;
			r = m;
		}
	}

	rightShiftCircularOnDistance(segmentInfoData, a, r, ldash, (ldash - f2),
			portionOneLastIndex);
	SplitMergeOnDistance(segmentInfoData, a, f1, r, r + rdash - f2,
			portionOneLastIndex);
	SplitMergeOnDistance(segmentInfoData, a, l + ldash - f2, ldash, last,
			portionOneLastIndex);
}

void initializeAndStart(vector<segmentInfo> segmentsInfos, vector<int>& a,
		bool dir) {

	int threadID = omp_get_thread_num();
	segmentInfo segmentInfoData = *(segmentsInfos.begin() + threadID);

	vector<std::pair<int, int> > tempBoundaries;
	for (int i = 0; i < segmentInfoData.noOfSegments; i++)
		tempBoundaries.push_back(
				std::make_pair(segmentInfoData.first[i],
						segmentInfoData.second[i]));

	//#pragma omp critical
	//    cout << threadID << " ok"
	//            << "first: " << segmentInfoData.boundries.begin()->first
	//            << " size: " << segmentInfoData.size << endl;

	int tempStart, tempEnd;
	typedef std::vector<std::pair<int, int> > my_vector;
	my_vector::iterator it;

	unsigned virtualStart1, virtualStart2, virtualEnd1, virtualEnd2, distance;
	unsigned currentSize = 0;
	it = tempBoundaries.begin();
	currentSize = it->second - it->first + 1;

	while (tempBoundaries.size() > 1) {
		it = tempBoundaries.begin();

		virtualStart1 = it->first;
		virtualStart2 = currentSize + virtualStart1;
		virtualEnd1 = it->first + currentSize - 1;
		virtualEnd2 = ((it + 1)->second - (it + 1)->first + 1) + virtualStart2;

		//adding the new merged segment size
		currentSize += (it + 1)->second - (it + 1)->first + 1;

		tempStart = it->first;
		tempEnd = tempStart + currentSize - 1;

		//array, virtual f1, virtual f2, size +1, virtual end1, distance
		SplitMergeOnDistance(segmentInfoData, a, virtualStart1, virtualStart2,
				virtualEnd2, virtualEnd1);
		tempBoundaries.erase(it, it + 2);
		tempBoundaries.insert(it, std::make_pair(tempStart, tempEnd));
	}
}

void multiwaySplitMergeOnDistance(vector<int>& arr,
		vector<segmentInfo> segmentsInfos) {

	//cout << "\n\nLazy merge based on split merging started...." << endl;
	const clock_t begin_time = clock();
	noOfSwaps1 = 0;
	int noOfThreads = segmentsInfos.size();

#pragma omp parallel num_threads(noOfThreads)
	{
		initializeAndStart(segmentsInfos, arr, true);
	}

	//    cout << "Lazy merge based on split merging ended." << endl;
	//    std::cout << "Elapsed time:  " << float( clock() - begin_time) / CLOCKS_PER_SEC << " seconds\n";
	//    std::cout << "Elapsed time:  " << float( clock() - begin_time) / 1000l << " milli-seconds\n\n";

	//printOnDistance(arr);
}
