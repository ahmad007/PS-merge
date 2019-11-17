#include "Guan_langston.h"
#include "splitMerge.h"
#include <iostream>
#include <vector>
#include <omp.h>
#include <cmath>
#include <omp.h>
#include <algorithm>
#include <fstream>
#include <climits>
#include <sys/time.h>

unsigned long noOfGuangSwaps = 0;
double accumTime = 0;

bool debug = false;
int blockSize = -1;
int counter2 = 0;



struct Breaker {
public:
	int blockID;
	//this index presents the entry with a value greater than the tail
	//of the pervious block.
	//So elements less than the breaker at the breaker record should only be considered
	//conclusion, use ((((breakIndex -1))) for merging with the previous blocks.
	int breakerIndex;

	bool operator<(Breaker other) const {
		return blockID < other.blockID;
	}
};

struct DisplacmentEntry {
public:
	int blockID;
	int displacement;
	int startIndex;
	int breakerID;
	int breakerBlockID;
};

void printArray(vector<int> a, int start, int end) {
	for (int i = start, j = 1; i <= end; i++, j++) {
		cout << a[i] << " ";
		if (j % blockSize == 0) {
			cout << "|";
		}
	}
	cout << endl;
}
/* This method sorts block by thier tails
 * @param a input list of 2-segment
 * @param limits list of block tails
 * @param numOfThreads number of used threads
 *
 */
vector<pair<int, int> > sortByBlocksTails(vector<int>& a,
		vector<std::pair<int, int> > limits, int numOfThreads) {

	vector<std::pair<int, int> > tails(numOfThreads);

	int totSize = 0;
	int firstListSize, secondListsize;

	totSize = firstListSize = limits[0].second - limits[0].first + 1;
	totSize += secondListsize = limits[1].second - limits[1].first + 1;
	blockSize = (totSize / numOfThreads);

	//add tails information
#pragma omp parallel num_threads(numOfThreads)
	//for( int i = 0 ; i < numOfThreads; i++)
	{
		int thread_id = omp_get_thread_num();
		int tailIndex = blockSize * (thread_id + 1) - 1;
		tails[thread_id] = std::make_pair(a[tailIndex], tailIndex);
	}

	if (debug) {
		cout << "unsorted tails (value, index):\n";
		for (int i = 0; i < tails.size(); i++)
			cout << "(" << tails[i].first << "," << tails[i].second << ") ";
		cout << endl;
	}

	sort(tails.begin(), tails.end());
	//sorting tails

	if (debug) {
		cout << "sorted tails(value, index):\n";
		for (int i = 0; i < tails.size(); i++)
			cout << "(" << tails[i].first << "," << tails[i].second << ") ";
		cout << endl;
	}

	vector<pair<int, int> > blocksSrcDest(numOfThreads);

	//set blocks source and dest
#pragma omp parallel num_threads(numOfThreads)
//    for( int i = 0 ; i < numOfThreads; i++)
	{
		int thread_id = omp_get_thread_num();
		int srcID = tails[thread_id].second / blockSize;
		blocksSrcDest[srcID] = std::pair<int, int>(srcID, thread_id);
	}

	if (debug) {
		cout << "sorted tails(srcblockid, destblockid):\n";
		for (int i = 0; i < tails.size(); i++)
			cout << "(" << blocksSrcDest[i].first << ","
					<< blocksSrcDest[i].second << ") ";
		cout << endl;
	}

	return blocksSrcDest;
}
/* This method swaps two block data
 * @param a input list of 2-segment
 * @param blocksSwappingSrcDest contains the src and dest block id
 * @param numOfThreads number of used threads
 *
 */
void swapBlocks(vector<int>& a, vector<pair<int, int> > blocksSwappingSrcDest,
		int numOfThreads) {
	int * currentValues = new int[numOfThreads];
#pragma omp parallel num_threads(numOfThreads)
	{
		int thread_id = omp_get_thread_num();
		for (int i = blockSize * thread_id; i < blockSize * (thread_id + 1);
				i++) {
			currentValues[thread_id] = a[i];
#pragma omp barrier 

			int destIndex = blockSize * blocksSwappingSrcDest[thread_id].second
					+ (i - blockSize * thread_id);
			a[destIndex] = currentValues[thread_id];
#if IS_DEBUG  == 1
#pragma omp critical
			{
				if (noOfGuangSwaps == ULONG_MAX) {
					counter2++;
					noOfGuangSwaps = 1;
				} else
				noOfGuangSwaps++;
			}
#endif
#pragma omp barrier 
		}
	}
	delete[] currentValues;
}

int binarySearchFindNextLargerValue(vector<int> a, int key, int start,
		int end) {

	//    int mid;
	//    while ((end - start) > 1) {
	//        mid = (start + end) / 2;
	//
	//        if (a[mid] <= key )
	//            start = mid + 1;
	//        else
	//            end = (mid - 1 <= start) ? end - 1 : mid - 1;
	//    }
	//    return end;
	for (int i = start; i <= end; i++)
		if (a[i] > key)
			return i;
	return end;

}

int binarySearchFindNextSmallerValue(vector<int> a, int key, int start,
		int end) {

	int mid;

	while ((end - start) > 1) {
		mid = (start + end) / 2;

		if (a[mid] >= key)
			end = mid - 1;
		else
			start = (mid + 1 >= end) ? start - 1 : mid + 1;
	}

	if (key > a[end])
		return end;
	return start;
}

int countLessThan(vector<int> a, int key, int start, int end) {
	int count = 0;
	for (int i = start; i <= end; i++)
		if (a[i] < key)
			count++;
	return count;
}

/* This method is used to sort breaker by blockid
 * @param a input list of 2-segment
 * @param numOfThreads number of used threads
 *
 */
int compareBreakers(const void *s1, const void *s2) {
	struct Breaker *b1 = (struct Breaker *) s1;
	struct Breaker *b2 = (struct Breaker *) s2;

	if (b1->blockID >= b2->blockID)
		return 1;
	else
		return -1;
}

/* This method finds breaker at each block
 * @param a input list of 2-segment
 * @param numOfThreads number of used threads
 *
 */
vector<Breaker> findBreakers(vector<int> a, int numOfThreads) {

	vector<Breaker> breakers;

	//    for (int thread_id = 0; thread_id < numOfThreads - 1; thread_id++) {
	//        //int thread_id = omp_get_thread_num();
	//        int tailIndex = blockSize * (thread_id + 1) - 1;
	//        if (a[tailIndex] > a[tailIndex + 1]) {
	//            Breaker b;
	//            b.blockID = thread_id + 1;
	//            b.breakerIndex = binarySearchFindNextLargerValue(a, a[tailIndex], tailIndex + 1, tailIndex + blockSize);
	//            breakers.push_back(b);
	//        }
	//    }

#pragma omp parallel num_threads(numOfThreads-1)
	{
		int thread_id = omp_get_thread_num();
		int tailIndex = blockSize * (thread_id + 1) - 1;
		if (a[tailIndex] > a[tailIndex + 1]) {
			Breaker b;
			b.blockID = thread_id + 1;
			b.breakerIndex = binarySearchFindNextLargerValue(a, a[tailIndex],
					tailIndex + 1, tailIndex + blockSize);

#pragma omp critical
			breakers.push_back(b);
		}
	}

	if (breakers.size() == 0)
		return breakers;

	sort(breakers.begin(), breakers.end());

	if (debug) {
		cout << "Breakers(blockid, breaker index):\n";
		for (int i = 0; i < breakers.size(); i++)
			cout << "(" << breakers[i].blockID << ","
					<< breakers[i].breakerIndex << ") ";
		cout << endl;
	}

	//modify breakers, that when some of the breaker block less than the breaker will not
	//be moved to series one ex: 2 11 22 23 29 29 |15 21 26 26 27 35 | only 15 and 21 will move
	//moved to series one while the breaker is 35

#pragma omp parallel num_threads(breakers.size())
	{
		int i = omp_get_thread_num();
		int start = blockSize * breakers[i].blockID;
		int end = breakers[i].breakerIndex - 1;
		//this to search only the pervious block of the breaker block
		int blockID = breakers[i].blockID - 1;

		for (int index = end; index >= start; index--) {

			int greaterIndex = binarySearchFindNextLargerValue(a, a[index],
					blockID * blockSize, (blockID + 1) * blockSize - 1);

			if ((((blockID + 1) * blockSize - 1) - greaterIndex)
					< (breakers[i].breakerIndex
							- breakers[i].blockID * blockSize - 1))
				breakers[i].breakerIndex = index;
			else
				break;
		}
	}

	//seq
	//    for (int i = 0; i < breakers.size(); i++) {
	//        int start = blockSize * breakers[i].blockID;
	//        int end = breakers[i].breakerIndex - 1;
	//        //this to search only the pervious block of the breaker block
	//        int blockID = breakers[i].blockID - 1;
	//
	//        for (int index = end; index >= start; index--) {
	//
	//            int greaterIndex = binarySearchFindNextLargerValue(a, a[index], blockID*blockSize, (blockID + 1) * blockSize - 1);
	//
	//            if ((((blockID + 1) * blockSize - 1) - greaterIndex) < (breakers[i].breakerIndex - breakers[i].blockID * blockSize - 1))
	//                breakers[i].breakerIndex = index;
	//            else
	//                break;
	//        }
	//    }

	if (debug) {
		cout << "Breakers(blockid, breaker index):\n";
		for (int i = 0; i < breakers.size(); i++)
			cout << "(" << breakers[i].blockID << ","
					<< breakers[i].breakerIndex << ") ";
		cout << endl;
	}

	return breakers;
}
/* This method shifts the data right as a rotation step
 * @param arr input list of 2-segment
 * @param shiftStartIndex shift starting index
 * @param shiftEndIndex shift ending index
 * @param shiftValue shift value
 */

void displacmentTable(std::vector<int> a, std::vector<Breaker> breakers,
		std::vector<std::pair<int, std::vector<DisplacmentEntry> > >& displcTble) {

	for (int i = 0; i < breakers.size(); i++) {
		vector<DisplacmentEntry> v;
		displcTble.push_back(std::pair<int, vector<DisplacmentEntry> >(i, v));
	}

	//#pragma omp parallel num_threads(breakers.size())
	//    {
	//        int i = omp_get_thread_num();
	//        int startBlock = (i == 0) ? 0 : breakers[i - 1].blockID;
	//        int start = blockSize * breakers[i].blockID;
	//        int end = breakers[i].breakerIndex - 1;
	//        int accumlativeDisplacment = 0;
	//
	//        for (int blockIndex = startBlock; blockIndex <= breakers[i].blockID - 1; blockIndex++) {
	//            int tailIndex = blockSize * (blockIndex + 1) - 1 - accumlativeDisplacment;
	//
	//            int displacmentIndex = 0;
	//            end = breakers[i].breakerIndex - 1;
	//            do {
	//                if (a[start] >= a[tailIndex])
	//                    displacmentIndex = 0;
	//                else
	//                    displacmentIndex = binarySearchFindNextSmallerValue(a, a[tailIndex], start, end) - blockSize * breakers[i].blockID + 1;
	//                end = start + displacmentIndex - 2;
	//            } while (displacmentIndex != 0 && a [start + displacmentIndex - 1] >= a[tailIndex]);
	//
	//
	//            accumlativeDisplacment = displacmentIndex;
	//            DisplacmentEntry dEntry;
	//            dEntry.blockID = blockIndex;
	//            dEntry.displacement = accumlativeDisplacment;
	//            dEntry.breakerID = i;
	//            dEntry.breakerBlockID = breakers[i].blockID;
	//            displcTble[i].second.push_back(dEntry);
	//        }
	//    }

#pragma omp parallel num_threads(breakers.size())
	{
		int i = omp_get_thread_num();
		int startBlock = (i == 0) ? 0 : breakers[i - 1].blockID;
		int start = blockSize * breakers[i].blockID;
		int end = breakers[i].breakerIndex - 1;
		int accumlativeDisplacment = 0;

		for (int blockIndex = startBlock; blockIndex <= breakers[i].blockID - 1;
				blockIndex++) {
			int tailIndex = blockSize * (blockIndex + 1) - 1
					- accumlativeDisplacment;

			int displacmentIndex = 0;
			end = breakers[i].breakerIndex - 1;
			if (a[start] >= a[tailIndex])
				displacmentIndex = 0;
			else {
				for (int series2Index = start + accumlativeDisplacment,
						series1Index = tailIndex; series2Index <= end;
						series2Index++, series1Index--)
					if (a[series1Index] > a[series2Index])
						accumlativeDisplacment++;
					else
						break;
			}

			DisplacmentEntry dEntry;
			dEntry.blockID = blockIndex;
			dEntry.displacement = accumlativeDisplacment;
			dEntry.breakerID = i;
			dEntry.breakerBlockID = breakers[i].blockID;
			displcTble[i].second.push_back(dEntry);
		}
	}

	//    for (int i = 0; i < breakers.size(); i++) {
	//        int startBlock = (i == 0) ? 0 : breakers[i - 1].blockID;
	//        int start = blockSize * breakers[i].blockID;
	//        int end = breakers[i].breakerIndex - 1;
	//        int accumlativeDisplacment = 0;
	//
	//        for (int blockIndex = startBlock; blockIndex <= breakers[i].blockID - 1; blockIndex++) {
	//            int tailIndex = blockSize * (blockIndex + 1) - 1 - accumlativeDisplacment;
	//
	//            int displacmentIndex = 0;
	//            end = breakers[i].breakerIndex - 1;
	//            if (a[start] >= a[tailIndex])
	//                displacmentIndex = 0;
	//            else {
	//                for (int series2Index = start + accumlativeDisplacment, series1Index = tailIndex; series2Index <= end; series2Index++, series1Index--)
	//                    if (a[series1Index] > a[series2Index])
	//                        accumlativeDisplacment++;
	//                    else
	//                        break;
	//            }
	//            //            do {
	//            //                if (a[start] >= a[tailIndex])
	//            //                    displacmentIndex = 0;
	//            //                else {
	//            //                   displacmentIndex = binarySearchFindNextSmallerValue(a, a[tailIndex], start, end) - blockSize * breakers[i].blockID + 1;
	//            //
	//            //                }
	//            //                end = start + displacmentIndex - 2;
	//            //            } while (displacmentIndex != 0 && a [start + displacmentIndex - 1] >= a[tailIndex]);
	//
	//
	//            // accumlativeDisplacment = displacmentIndex;
	//
	//            DisplacmentEntry dEntry;
	//            dEntry.blockID = blockIndex;
	//            dEntry.displacement = accumlativeDisplacment;
	//            dEntry.breakerID = i;
	//            dEntry.breakerBlockID = breakers[i].blockID;
	//            displcTble[i].second.push_back(dEntry);
	//        }
	//    }

	if (debug) {
		cout
				<< "displacement tables one per line (blockid, displacment value):\n";
		for (int i = 0; i < displcTble.size(); i++) {
			for (int j = 0; j < displcTble[i].second.size(); j++) {
				cout << "(" << displcTble[i].second[j].blockID << ","
						<< displcTble[i].second[j].displacement << ") ";
			}
			cout << endl;
		}
	}
	if (debug) {
		cout
				<< "displacement tables one per line (blockid, displacment value,startIndex, breakerID):\n";
		for (int i = 0; i < displcTble.size(); i++) {
			for (int j = 0; j < displcTble[i].second.size(); j++) {
				cout << "(" << displcTble[i].second[j].blockID << ","
						<< displcTble[i].second[j].displacement << ","
						<< displcTble[i].second[j].startIndex << ","
						<< displcTble[i].second[j].breakerID << ") ";
			}
			cout << endl;
		}
	}

}
/* This method shifts the data right as a rotation step
 * @param arr input list of 2-segment
 * @param shiftStartIndex shift starting index
 * @param shiftEndIndex shift ending index
 * @param shiftValue shift value
 */
void guangRightShiftCircular(vector<int>& arr, int shiftStartIndex,
		int shiftEndIndex, int shiftValue) {
	const clock_t begin_time = clock();
	int temp = 0;
	for (int i = shiftStartIndex, j = shiftEndIndex; --j > i; i++) {
		temp = arr[i];
		arr[i] = arr[j];
		arr[j] = temp;
#if IS_DEBUG == 1
#pragma omp critical
		{
			if (noOfGuangSwaps == ULONG_MAX) {
				counter2++;
				noOfGuangSwaps = 1;
			} else
			noOfGuangSwaps++;
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
			if (noOfGuangSwaps == ULONG_MAX) {
				counter2++;
				noOfGuangSwaps = 1;
			} else
			noOfGuangSwaps++;
		}
#endif
	}

	for (int i = shiftStartIndex + shiftValue, j = shiftEndIndex; --j > i;
			i++) {
		temp = arr[i];
		arr[i] = arr[j];
		arr[j] = temp;
#if IS_DEBUG  == 1
#pragma omp critical
		{
			if (noOfGuangSwaps == ULONG_MAX) {
				counter2++;
				noOfGuangSwaps = 1;
			} else
			noOfGuangSwaps++;
		}
#endif
	}

}
/* This method rotates block as step before the data movement
 * @param a input list of 2-segment
 * @param breakers list of breakers
 * @param displcTble displacment table which shows which element should be moved from which to block to the destension block
 * @param numOfBlocks number of blocks
 */

void rotateBlocks(vector<int>& a, vector<Breaker> breakers,
		std::vector<std::pair<int, std::vector<DisplacmentEntry> > >& displcTble,
		int numOfBlocks) {

#pragma omp parallel num_threads(numOfBlocks)
	//for ( int i = 0 ; i < numOfBlocks; i++)
	{
		int i = omp_get_thread_num();

		int blockEnd = blockSize * (i + 1) - 1;
		int blockStart = blockSize * (i);

		//set block start to the breaker if the block contains a breaker
		for (int breakerIndex = 0; breakerIndex < breakers.size();
				breakerIndex++)
			if (i == breakers[breakerIndex].blockID)
				blockStart = breakers[breakerIndex].breakerIndex;
		// cout<<"breaker: "<<i<<" : ";
		//find displacment
		for (int displcmentTableIndex = 0;
				displcmentTableIndex < displcTble.size();
				displcmentTableIndex++)
			for (int entryIndex = 0;
					entryIndex < displcTble[displcmentTableIndex].second.size();
					entryIndex++)
				if (displcTble[displcmentTableIndex].second[entryIndex].blockID
						== i) {
					//      cout<<"start: " << blockStart << " end: " <<blockEnd+1<<endl;
					displcTble[displcmentTableIndex].second[entryIndex].startIndex =
							blockStart
									+ displcTble[displcmentTableIndex].second[entryIndex].displacement
									- 1;
					guangRightShiftCircular(a, blockStart, blockEnd + 1,
							displcTble[displcmentTableIndex].second[entryIndex].displacement);
				}
		//  pirntArray(a);
	}

	//sequential version
	//    for (int i = 0; i < numOfBlocks; i++) {
	//
	//        int blockEnd = blockSize * (i + 1) - 1;
	//        int blockStart = blockSize * (i);
	//
	//        //set block start to the breaker if the block contains a breaker
	//        for (int breakerIndex = 0; breakerIndex < breakers.size(); breakerIndex++)
	//            if (i == breakers[breakerIndex].blockID)
	//                blockStart = breakers[breakerIndex].breakerIndex;
	//        // cout<<"breaker: "<<i<<" : ";
	//        //find displacment
	//        for (int displcmentTableIndex = 0; displcmentTableIndex < displcTble.size(); displcmentTableIndex++)
	//            for (int entryIndex = 0; entryIndex < displcTble[displcmentTableIndex].second.size(); entryIndex++)
	//                if (displcTble[displcmentTableIndex].second[entryIndex].blockID == i) {
	//                    //      cout<<"start: " << blockStart << " end: " <<blockEnd+1<<endl;
	//                    guangRightShiftCircular(a, blockStart, blockEnd + 1, displcTble[displcmentTableIndex].second[entryIndex].displacement);
	//                }
	//        //  pirntArray(a);
	//    }
}

/* This method calculates the displacment table
 * @param displcTble displacment table which shows which element should be moved from which to block to the destension block
 * @param flatDisplacementTable displacment table which shows which element should be moved from which to block to the destension block
 * @param numOfBlocks number of blocks
 */

void getFlatDisplacmentTable(
		std::vector<std::pair<int, std::vector<DisplacmentEntry> > > displcTble,
		vector<DisplacmentEntry>& flatDisplacementTable, int numOfBlocks) {
	for (int blockID = 0; blockID < numOfBlocks; blockID++) {
		DisplacmentEntry entry;
		flatDisplacementTable.push_back(entry);
	}
#pragma omp parallel num_threads(numOfBlocks)
//for( int i =0 ; i < numOfBlocks; i++)
	{
		int blockID = omp_get_thread_num();
		bool found = false;
		for (int breakerID = 0; breakerID < displcTble.size(); breakerID++) {
			for (int entryID = 0; entryID < displcTble[breakerID].second.size();
					entryID++)
				if (displcTble[breakerID].second[entryID].blockID == blockID) {
					found = true;
					flatDisplacementTable[blockID].blockID = blockID;
					flatDisplacementTable[blockID].displacement =
							displcTble[breakerID].second[entryID].displacement;
					flatDisplacementTable[blockID].startIndex =
							displcTble[breakerID].second[entryID].startIndex;
					flatDisplacementTable[blockID].breakerID = breakerID;
					flatDisplacementTable[blockID].breakerBlockID =
							displcTble[breakerID].second[entryID].breakerBlockID;
					break;
				}
			if (found)
				break;
			else {
				flatDisplacementTable[blockID].blockID = blockID;
				flatDisplacementTable[blockID].displacement = -1;
				flatDisplacementTable[blockID].startIndex = -1;
				flatDisplacementTable[blockID].breakerID = -1;
				flatDisplacementTable[blockID].breakerBlockID = -1;
			}
		}
	}

	//seq version
	//    for (int blockID = 0; blockID < numOfBlocks; blockID++) {
	//        bool found = false;
	//        for (int breakerID = 0; breakerID < displcTble.size(); breakerID++) {
	//            for (int entryID = 0; entryID < displcTble[breakerID].second.size(); entryID++)
	//                if (displcTble[breakerID].second[entryID].blockID == blockID) {
	//                    found = true;
	//                    flatDisplacementTable[blockID].blockID = blockID;
	//                    flatDisplacementTable[blockID].displacement = displcTble[breakerID].second[entryID].displacement;
	//                    flatDisplacementTable[blockID].startIndex = displcTble[breakerID].second[entryID].startIndex;
	//                    flatDisplacementTable[blockID].breakerID = breakerID;
	//                    break;
	//                }
	//            if (found)
	//                break;
	//            else {
	//                flatDisplacementTable[blockID].blockID = blockID;
	//                flatDisplacementTable[blockID].displacement = -1;
	//                flatDisplacementTable[blockID].startIndex = -1;
	//                flatDisplacementTable[blockID].breakerID = -1;
	//            }
	//        }
	//    }
}
/* This method does the data movement of Guang's algorithm, it is element synchronized
 * @param a input list of 2-segment
 * @param displcTble displacment table which shows which element should be moved from which to block to the destension block
 * @param flatDisplacementTable displacment table which shows which element should be moved from which to block to the destension block
 * @param numOfBlocks number of blocks
 */

void dataMovement(vector<int>& a,
		std::vector<std::pair<int, std::vector<DisplacmentEntry> > > displcTble,
		vector<DisplacmentEntry>& flatDisplacementTable, int numOfBlocks) {
	getFlatDisplacmentTable(displcTble, flatDisplacementTable, numOfBlocks);

	if (debug) {
		cout
				<< "Flat dispTable(blockid, displacment value,startIndex, breakerID):\n";
		for (int i = 0; i < flatDisplacementTable.size(); i++) {
			cout << "(" << flatDisplacementTable[i].blockID << ","
					<< flatDisplacementTable[i].displacement << ","
					<< flatDisplacementTable[i].startIndex << ","
					<< flatDisplacementTable[i].breakerID << ") ";
		}
		cout << endl;
	}

	vector<int> series1Data(numOfBlocks); //one entry for each block
	vector<int> series2Data(displcTble.size()); //one entry for each breaker
	vector<int> blockIDReceivingRecordsFromSeries2(displcTble.size()); //one entry for each breaker

	for (int breakerID = 0; breakerID < displcTble.size(); breakerID++) {
		for (int entryID = 0; entryID < displcTble[breakerID].second.size();
				entryID++)
			if (displcTble[breakerID].second[entryID].displacement > 0) {
				blockIDReceivingRecordsFromSeries2[breakerID] =
						displcTble[breakerID].second[entryID].blockID;
				break;
			}
	}

	int maxDisplacement = -1;

	for (int breakID = 0; breakID < displcTble.size(); breakID++) {
		if (displcTble[breakID].second[displcTble[breakID].second.size() - 1].displacement
				> maxDisplacement)
			maxDisplacement =
					displcTble[breakID].second[displcTble[breakID].second.size()
							- 1].displacement;
	}

	//seq version

	for (int iteration = 0; iteration < maxDisplacement; iteration++) {

#pragma omp parallel num_threads(numOfBlocks)
		//for( int i =0 ; i < numOfBlocks; i++)
		{
			int blockID = omp_get_thread_num();
			if (flatDisplacementTable[blockID].displacement > iteration) {
				//testing if the next block the break block
				if (flatDisplacementTable[blockID].breakerBlockID - blockID
						== 1) {
					//write current block record to breakers
					series2Data[flatDisplacementTable[blockID].breakerID] =
							a[flatDisplacementTable[blockID].startIndex
									- iteration];
					//write series2 to this breaker first block
					series1Data[blockIDReceivingRecordsFromSeries2[flatDisplacementTable[blockID].breakerID]] =
							a[flatDisplacementTable[blockID].breakerBlockID
									* blockSize + iteration];
				} else {
					series1Data[blockID + 1] =
							a[flatDisplacementTable[blockID].startIndex
									- iteration];
				}
			}
		}

		//        for (int blockID = 0; blockID < numOfBlocks; blockID++) {
		//            if (flatDisplacementTable[blockID].displacement > iteration) {
		//                //testing if the next block the break block
		//                if (flatDisplacementTable[blockID].breakerBlockID - blockID == 1) {
		//                    //write current block record to breakers
		//                    series2Data[flatDisplacementTable[blockID].breakerID] = a[flatDisplacementTable[blockID].startIndex - iteration];
		//                    //write series2 to this breaker first block
		//                    series1Data[blockIDReceivingRecordsFromSeries2[flatDisplacementTable[blockID].breakerID]] = a[flatDisplacementTable[blockID].breakerBlockID * blockSize + iteration];
		//                } else {
		//                    series1Data[blockID + 1] = a[flatDisplacementTable[blockID].startIndex - iteration];
		//                }
		//            }
		//        }

		//update next block to receive from breaker
		for (int breakerID = 0; breakerID < displcTble.size(); breakerID++) {
			for (int entryID = 0; entryID < displcTble[breakerID].second.size();
					entryID++)
				if (displcTble[breakerID].second[entryID].displacement
						> (iteration + 1)) {
					blockIDReceivingRecordsFromSeries2[breakerID] =
							displcTble[breakerID].second[entryID].blockID;
					break;
				}
		}

		if (debug) {
			cout << "ok:\n";
			for (int i = 0; i < series1Data.size(); i++)
				cout << series1Data[i] << ",";
			cout << endl;

			for (int i = 0; i < series2Data.size(); i++)
				cout << series2Data[i] << ",";
			cout << endl;
		}

		//write cycle
#pragma omp parallel num_threads(numOfBlocks)
		//for( int i =0 ; i < numOfBlocks; i++)
		{
			int blockID = omp_get_thread_num();
			if (flatDisplacementTable[blockID].displacement > iteration) {
				if (flatDisplacementTable[blockID].breakerBlockID - blockID
						== 1) {
					a[flatDisplacementTable[blockID].breakerBlockID * blockSize
							+ iteration] =
							series2Data[flatDisplacementTable[blockID].breakerID];
				}
				a[flatDisplacementTable[blockID].startIndex - iteration] =
						series1Data[blockID];
#if IS_DEBUG == 1
#pragma omp critical
				{
					if (noOfGuangSwaps == ULONG_MAX) {
						counter2++;
						noOfGuangSwaps = 1;
					} else
					noOfGuangSwaps++;
				}
#endif
			}
		}

		//        for (int blockID = 0; blockID < numOfBlocks; blockID++) {
		//            if (flatDisplacementTable[blockID].displacement > iteration) {
		//                if (flatDisplacementTable[blockID].breakerBlockID - blockID == 1) {
		//                    a[flatDisplacementTable[blockID].breakerBlockID * blockSize + iteration] = series2Data[flatDisplacementTable[blockID].breakerID];
		//                }
		//                a[flatDisplacementTable[blockID].startIndex - iteration] = series1Data[blockID];
		//            }
		//        }

		for (int i = 0; i < series1Data.size(); i++)
			series1Data[i] = -1;

	}
}

/* This method reverses the specified range of any block
 * @param a input list of 2-segment
 * @param startIndex start index of reverse
 * @param endIndex end index of reverse
 *
 */

void reverse(vector<int>& a, int startIndex, int endIndex) {
	for (int i = startIndex, j = endIndex; i < j; i++, j--) {
		int temp = a[i];
		a[i] = a[j];
		a[j] = temp;
#if IS_DEBUG  == 1
#pragma omp critical
		{
			if (noOfGuangSwaps == ULONG_MAX) {
				counter2++;
				noOfGuangSwaps = 1;
			} else
			noOfGuangSwaps++;
		}
#endif
	}
}

/* This method reverses sub-block elements
 * @param a input list of 2-segment
 * @param breakers list of breakers founded by gunag's algorithm
 * @param flatDisplacementTable displacment table which shows which element should be moved from which to block to the destension block
 * @param numOfBlocks number of blocks
 */

void reverseSubblocks(vector<int>& a, vector<Breaker> breakers,
		vector<DisplacmentEntry>& flatDisplacementTable, int numOfBlocks) {

#pragma omp parallel num_threads(numOfBlocks)
	//for( int i =0 ; i < numOfBlocks; i++)
	{
		int blockID = omp_get_thread_num();
		if (flatDisplacementTable[blockID].displacement > 1) {
			int startReverseIndex = -1;
			int reverseSize = -1;
			if (blockID > 0
					&& flatDisplacementTable[blockID - 1].breakerBlockID
							== flatDisplacementTable[blockID].breakerBlockID) {
				reverseSize = flatDisplacementTable[blockID].displacement
						- flatDisplacementTable[blockID - 1].displacement;
				startReverseIndex = flatDisplacementTable[blockID].startIndex
						- flatDisplacementTable[blockID - 1].displacement
						- reverseSize + 1;
				reverse(a, startReverseIndex,
						startReverseIndex + reverseSize - 1);
			} else
				reverse(a,
						flatDisplacementTable[blockID].startIndex
								- flatDisplacementTable[blockID].displacement
								+ 1, flatDisplacementTable[blockID].startIndex);
		}
	}
//    for (int blockID = 0; blockID < numOfBlocks; blockID++) {
//        if (flatDisplacementTable[blockID].displacement > 1) {
//            int startReverseIndex = -1;
//            int reverseSize = -1;
//            if (blockID > 0 && flatDisplacementTable[blockID - 1].breakerBlockID == flatDisplacementTable[blockID].breakerBlockID) {
//                reverseSize = flatDisplacementTable[blockID].displacement - flatDisplacementTable[blockID - 1].displacement;
//                startReverseIndex = flatDisplacementTable[blockID].startIndex - flatDisplacementTable[blockID - 1].displacement - reverseSize + 1;
//                reverse(a, startReverseIndex, startReverseIndex + reverseSize - 1);
//            } else
//                reverse(a, flatDisplacementTable[blockID].startIndex - flatDisplacementTable[blockID].displacement + 1, flatDisplacementTable[blockID].startIndex);
//        }
//    }

	for (int breakerID = 0; breakerID < breakers.size(); breakerID++) {
		reverse(a, blockSize * breakers[breakerID].blockID,
				breakers[breakerID].breakerIndex - 1);
	}
}

/* This method check if the current block contains breaker or not
 * @param breakers list of breakers founded by gunag's algorithm
 * @param blockID block id (thread id)
 *
 */

int isBreakerBlock(vector<Breaker> breakers, int blockID) {
	int start = 0, end = breakers.size();
	int mid;
	while (start <= end) {
		mid = (start + end) / 2;
		if (breakers[mid].blockID == blockID)
			return mid;
		else if (breakers[mid].blockID < blockID)
			start = mid + 1;
		else
			end = mid - 1;
	}
	return -1;
}
/* This method checks if the array in order or not, it is used after completing the merging process
 * @param a list of input, contains two segments\
 *
 */

bool checkArrayOrder(vector<int>& a) {
	for (int i = 1; i < a.size(); i++)
		if (a[i] < a[i - 1]) {
			cout << "error at: a[" << i - 1 << "]=" << a[i - 1] << " a[" << i
					<< "]=" << a[i] << endl;
			return false;
		}
	cout << "Checking is done successfully!\n";
	return true;
}
/* This method doing the local merge of adjacent elements requiring merging by guang's algorithm
 * @param a list of input, contains two segments
 * @param breakers list of breakers founded by gunag's algorithm
 * @param flatDisplacementTable displacment table which shows which element should be moved from which to block to the destension block
 * @param numOfBlocks number of blocks
 */

void localMerging(vector<int>& a, vector<Breaker> breakers,
		vector<DisplacmentEntry> flatDisplacementTable, int numOfBlocks) {

#pragma omp parallel num_threads(numOfBlocks)
	//for( int i =0 ; i < numOfBlocks; i++)
	{
		int blockID = omp_get_thread_num();
		vector<int> indexes;
		indexes.push_back(blockSize * blockID);
		int index = isBreakerBlock(breakers, blockID);
		if (index > -1)
			indexes.push_back(breakers[index].breakerIndex);
		if (flatDisplacementTable[blockID].displacement > 0
				&& flatDisplacementTable[blockID].startIndex > -1) {
			if (blockID > 0
					&& flatDisplacementTable[blockID - 1].breakerBlockID
							== flatDisplacementTable[blockID].breakerBlockID) {
				indexes.push_back(
						flatDisplacementTable[blockID].startIndex
								- flatDisplacementTable[blockID - 1].displacement
								+ 1);
			} else
				indexes.push_back(
						flatDisplacementTable[blockID].startIndex + 1);
		}
		indexes.push_back(blockSize * (blockID + 1));

		if (debug) {
			if (indexes.size() > 2) {
				for (int i = 0; i < indexes.size() - 1; i++)
					cout << "(" << indexes[i] << "," << indexes[i + 1] - 1
							<< ") ";
				cout << endl;
			} else
				cout << "no merging is required for this block!\n";
		}
		while (indexes.size() > 2) {
			for (int i = 0; i < indexes.size() - 1; i++) {
				iterativeSplitMerge(a, indexes[0], indexes[1], indexes[2]);
				indexes.erase(indexes.begin() + 1);
			}
		}
	}
}

/* This method doing the 2-way guang's algorithm; for input size is multiple of number of blocks(threads)
 * @param a list of input, contains two segments
 * @param limits range of each segment
 * @param numOfBlocks number of threads available for use
 *
 */
//actual indexes
void guanLangston2wayMerging(vector<int>& a,
		vector<std::pair<int, int> > limits, int numOfBlocks) {

	vector<pair<int, int> > blocksSwappingSrcDest = sortByBlocksTails(a, limits,
			numOfBlocks);
	swapBlocks(a, blocksSwappingSrcDest, numOfBlocks);
	if (debug)
		printArray(a, limits[0].first, limits[1].second);
	vector<Breaker> breakers = findBreakers(a, numOfBlocks);
	if (breakers.size() == 0)
		return;

	std::vector<std::pair<int, std::vector<DisplacmentEntry> > > displcTble;
	displacmentTable(a, breakers, displcTble);

	rotateBlocks(a, breakers, displcTble, numOfBlocks);
	if (debug)
		printArray(a, limits[0].first, limits[1].second);

	vector<DisplacmentEntry> flatDisplacementTable;
	dataMovement(a, displcTble, flatDisplacementTable, numOfBlocks);
	if (debug) {
		cout << "after movement:\n";
		printArray(a, limits[0].first, limits[1].second);
	}

	reverseSubblocks(a, breakers, flatDisplacementTable, numOfBlocks);

	localMerging(a, breakers, flatDisplacementTable, numOfBlocks);

}

/* This method doing the 2-way guang's algorithm; for arbitrary size
 * @param a list of input, contains two segments
 * @param limits range of each segment
 * @param numOfBlocks number of threads available for use
 *
 */

void guanLangstonMergingArbitrarySize(vector<int>& a,
		vector<std::pair<int, int> > limits, int numOfBlocks) {

	int L5StartIndex, L6StartIndex;
	int L1Size = limits[0].second - limits[0].first + 1;
	int L2Size = limits[1].second - limits[1].first + 1;

	int totSize = L1Size + L2Size;
	blockSize = (totSize / numOfBlocks);

	int L5Size = L1Size % blockSize;
	int L6Size = L2Size % blockSize;

	if (L5Size != 0)
		L5StartIndex = limits[0].second - L1Size % blockSize + 1;
	if (L6Size != 0)
		L6StartIndex = limits[1].second - L2Size % blockSize + 1;

	if (L5Size > 0) {
		guangRightShiftCircular(a, L5StartIndex, limits[1].second - L6Size + 1,
				L2Size - L6Size);
		limits[0].second = L5StartIndex - 1;
		limits[1].first = limits[0].second + 1;
		limits[1].second -= L5Size;
	}

	if (debug)
		printArray(a, 0, a.size() - 1);

	if (L6Size > 0) {
		limits[1].second -= L6Size;
	} else {
		if (L5Size == 0 && blockSize * numOfBlocks < a.size()) {
			limits[1].second = blockSize * numOfBlocks - 1;
			L6Size = a.size() - blockSize * numOfBlocks;
			L6StartIndex = blockSize * numOfBlocks;
		}
	}

	vector<std::pair<int, int> > L5L6limits;

	if (L5Size > 0)
		L5L6limits.push_back(
				std::make_pair(limits[1].second + 1,
						limits[1].second + L5Size));

	if (L6Size > 0)
		L5L6limits.push_back(
				std::make_pair(L6StartIndex, L6StartIndex + L6Size - 1));

	if (debug) {
		cout << "L1 & L2 limits:\n";
		cout << "(" << limits[0].first << "," << limits[0].second << ") " << "("
				<< limits[1].first << "," << limits[1].second << ") " << endl;
		cout << "L5 & L6 limits:\n";
		if (L5Size > 0)
			cout << "(" << L5L6limits[0].first << "," << L5L6limits[0].second
					<< ") ";
		if (L6Size > 0)
			cout << "(" << L5L6limits[L5L6limits.size() - 1].first << ","
					<< L5L6limits[L5L6limits.size() - 1].second << ") " << endl;
	}

	//sequentially merge L5& L6 >=L4
	if (L5L6limits.size() > 1)
		iterativeSplitMerge(a, L5L6limits[0].first, L5L6limits[1].first,
				L5L6limits[1].second + 1);

	//padding L4 to be of size n/k with last value
	int s = 0, paddingSize = 0;

	if (L5L6limits.size() == 2) {
		s = L5L6limits[0].second - L5L6limits[0].first + 1;
		s += L5L6limits[1].second - L5L6limits[1].first + 1;
	} else if (L5L6limits.size() == 1) {
		s = L5L6limits[0].second - L5L6limits[0].first + 1;
	}

	if (L5L6limits.size() > 0) {
		paddingSize = blockSize
				- ((s % blockSize == 0) ? blockSize : (s % blockSize));
		L5L6limits[L5L6limits.size() - 1].second += paddingSize;
		for (int i = 0; i < paddingSize; i++)
			a.push_back(a[a.size() - 1]);
	}

	//parallely merge L1&L2=>L3
	numOfBlocks = limits[1].second / blockSize + 1;
	if (numOfBlocks > 1)
		guanLangston2wayMerging(a, limits, numOfBlocks);

	if (debug) {
		cout << "l1 l2 done\n";
		printArray(a, limits[0].first, limits[1].second);
	}

	if (L5Size > 0 || L6Size > 0) {
		vector<std::pair<int, int> > L3L4limits;
		L3L4limits.push_back(std::make_pair(limits[0].first, limits[1].second));

		if (L5L6limits.size() > 1) {
			L3L4limits.push_back(
					std::make_pair(L5L6limits[0].first, L5L6limits[1].second));

		} else {
			L3L4limits.push_back(
					std::make_pair(L5L6limits[0].first, L5L6limits[0].second));
		}
		numOfBlocks += (L5Size + L6Size + paddingSize) / blockSize;

		//        cout<<"before mergering\n";
		//        printArray(a, L3L4limits[0].first, L3L4limits[1].second);
		//parallely merge L3 & L4
		guanLangston2wayMerging(a, L3L4limits, numOfBlocks);
		for (int i = 0; i < paddingSize; i++)
			a.erase(a.end() - 1);
	}
}

/* This method doing the k-way guang's algorithm; for each 2-way merging task the availabe threads are used
 * so it is performing a tournament tree of 2-way sequentially where each merging task is done using multithread
 * @param a list of input lists, each two list is an entry after each phase this list size is halved
 * @param segmentsEnds list of list of end segments
 * @param numOfAvailableThreads number of threads available for use
 * @param myfile file stream to write the results
 *
 */

void guanLangstonMultiwayMerging(vector<vector<int> >& a,
		vector<vector<int> > & segmentsEnds, int numOfAvailableThreads) {

	noOfGuangSwaps = 0;

	vector<std::pair<int, int> > limits; //actual index
	int size = a.size();

	for (int i = 0; i < segmentsEnds.size(); i++) {
		limits.push_back(std::make_pair(0, segmentsEnds[i][0]));
		if (segmentsEnds[i].size() == 2)
			limits.push_back(
					std::make_pair(segmentsEnds[i][0] + 1, segmentsEnds[i][1]));
	}

	int numOfSegments = limits.size();

	int noOfMergingTasks, tempStart, tempEnd;
	typedef std::vector<std::pair<int, int> > my_vector;
	my_vector::iterator it;
	int phases = std::ceil(std::log(limits.size()) / std::log(2));

	struct timeval startwtime1, endwtime1;
	double seq_time1;
	//double accumTime = 0;

	for (int i = 0; i < phases; i++) {

		gettimeofday(&startwtime1, NULL);
		noOfMergingTasks = limits.size() / 2;

		int usedThreads = numOfAvailableThreads / noOfMergingTasks;
		resetNoOfSwaps();
		//#pragma omp parallel num_threads(noOfMergingTasks)
		for (int i = 0; i < noOfMergingTasks; i++) {
			//int i = omp_get_thread_num();
			vector<std::pair<int, int> > localLimits;
			localLimits.push_back(limits[2 * i]);
			localLimits.push_back(limits[2 * i + 1]);
			guanLangstonMergingArbitrarySize(a[i], localLimits,
					numOfAvailableThreads);
			gettimeofday(&endwtime1, NULL);

			seq_time1 = (double) ((endwtime1.tv_usec - startwtime1.tv_usec)
					/ 1.0e6 + endwtime1.tv_sec - startwtime1.tv_sec);
			accumTime += seq_time1;

		}

		noOfGuangSwaps += getNumOfSwapsOfSplitMerge();

		for (int threadId = 0; threadId < noOfMergingTasks; threadId++) {
			it = limits.begin();
			tempStart = (it + threadId)->first;
			tempEnd = (it + (threadId + 1))->second;
			limits.erase(it + threadId, it + (threadId + 2));
			limits.insert(it + threadId, std::make_pair(tempStart, tempEnd));
		}
		for (int threadID = 1; threadID < limits.size(); threadID += 2) {
			limits[threadID].second = limits[threadID - 1].second
					+ (limits[threadID].second - limits[threadID].first + 1);
			limits[threadID].first = limits[threadID - 1].second + 1;
		}
		for (int i = 1; i < a.size(); i++) {
			for (int k = 0; k < a[i].size(); k++)
				a[i - 1].push_back(a[i][k]);
			a.erase(a.begin() + i);
		}
	}
}
