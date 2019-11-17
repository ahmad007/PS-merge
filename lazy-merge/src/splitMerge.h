/* 
 * File:   splitMerge.h
 * Author: ahmad
 *
 * Created on December 26, 2013, 2:59 PM
 */

#ifndef SPLITMERGE_H
#define	SPLITMERGE_H
#include <vector>
#include <fstream>


using namespace std;

void rightShiftCircular(vector<int>& arr, int shiftStartIndex, int shiftEndIndex, int shiftValue);
void SplitMerge(vector<int>& arr, int f1, int f2, int last);
void testSplitMerge();
void multiWaySplitMerge();
void multiWaySplitMerge(vector<int>& a, vector<int > segmentsEnds, int usedThreads, ofstream* myfile);
unsigned long getNumOfSwapsOfSplitMerge();
bool check(vector<int> a);
void iterativeSplitMerge(vector<int>& a, int f1, int f2, int last);
void resetNoOfSwaps();
unsigned long getNumOfComparisonsSplitMerge();
void resetNumOfComparisonsSplitMerge();


#endif	/* SPLITMERGE_H */



