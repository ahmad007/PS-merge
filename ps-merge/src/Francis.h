/* 
 * File:   Francis.h
 * Author: ahmad
 *
 * Created on December 15, 2013, 4:20 PM
 */

#ifndef FRANCIS_H
#define	FRANCIS_H



#endif	/* FRANCIS_H */

#include <vector>

using namespace std;

void FrancisPartitioning();
void FrancisPartitioning(vector<int > a, const int noOfSegments, const int noOfThreads, const vector<int > segmentsEnds, vector<std::pair<int,vector<int> > >& boundaries);
int * getBoundries(int * a, int noOfSegments, int partitionLen, int * b);
unsigned getNumOfComparisonsFrancis();
