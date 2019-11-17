/* 
 * File:   splitMergeOnDistance.h
 * Author: ahmad
 *
 * Created on December 26, 2013, 3:06 PM
 */

#ifndef SPLITMERGEONDISTANCE_H
#define	SPLITMERGEONDISTANCE_H



#endif	/* SPLITMERGEONDISTANCE_H */

#include <vector>
#include "segmentInfo.h"

using namespace std;

#define IS_DEBUG 0


void rightShiftCircularOnDistanceCache(vector<int>& arr, int shiftStartIndex, int shiftEndIndex, int shiftValue);
void SplitMergeOnDistance(vector<int>& a, int f1, int f2, int last);
void multiwaySplitMergeOnDistanceCache(vector<int>& arr, vector<segmentInfo > segmentsInfos, long& threadTime);
void multiwaySplitMergeOnDistanceMap(vector<int>& arr, vector<segmentInfo > segmentsInfos);
void multiwaySplitMergeOnDistance(vector<int>& arr, vector<segmentInfo > segmentsInfos);
long getNumOfSwapsSplitMergeOnDistance();
bool checkArray(vector<int>& indexes,vector<int>& a, vector<segmentInfo > segmentsInfos);
