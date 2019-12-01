/* 
 * File:   generateData.h
 * Author: ahmad
 *
 * Created on December 28, 2013, 2:30 PM
 */

#ifndef GENERATEDATA_H
#define	GENERATEDATA_H



#endif	/* GENERATEDATA_H */


#include <vector>
using namespace std;
void generateArray(vector<int >& arr, int n, const vector<int > segmentsEnds);
void generateArray(vector<int >& arr, int n, int numOfPortions, const vector<int > segmentsEnds);
void setSegmentsEndEqually(int size, int numOfSegments, vector<int >& segmentsEnds);
void generateDataForGuang(vector<vector<int> >& arr, int n, int numOfSegments, vector<vector<int > >& segmentsEnds);
