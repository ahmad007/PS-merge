/* 
 * File:   Guan_langston.h
 * Author: ahmad
 *
 * Created on March 16, 2014, 2:50 PM
 */

#ifndef GUAN_LANGSTON_H
#define	GUAN_LANGSTON_H

#include <vector>
#include <fstream>
using namespace std;

extern unsigned long noOfGuangSwaps;
extern double accumTime;

#define IS_DEBUG 0


void guanLangstonMultiwayMerging(vector<vector<int> >& a, vector< vector<int > > & segmentsEnds, int numOfAvailableThreads);
void guanLangston2wayMerging(vector<int>& a, vector<std::pair<int, int> > limits, int numOfBlocks);
void guanLangstonMergingArbitrarySize(vector<int>& a, vector<std::pair<int, int> >  limits,int  numOfBlocks);
bool checkArrayOrder(vector<int>& a);
#endif	/* GUAN_LANGSTON_H */

