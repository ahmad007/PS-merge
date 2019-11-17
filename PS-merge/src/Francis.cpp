#include <iostream>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <ctime>
#include <sys/time.h>
#include "Francis.h"


using namespace std;
unsigned numOfComparisonsFransic = 0;

vector<int> counter1;
vector<double> seq_time1;

unsigned getNumOfComparisonsFrancis() {
	return numOfComparisonsFransic;
}
vector<int> getBoundries(vector<int> a, int noOfSegments, int partitionLen,
		vector<int> b, int id) {
	struct timeval startwtime1, endwtime1;
		gettimeofday(&startwtime1, NULL);
	int deltaMax, deltaMin, delta;
	int rmin = 100, lmax = -10;
	int segContainLmax = -1, segContainRmin = 100;
	int temp;

	int threadid = id;//omp_get_thread_num();


//	 vector<int> s(noOfSegments, 0);
//	 vector<int> l(noOfSegments, 0);
//	 vector<int> r(noOfSegments, 0);

	int* s = new int[noOfSegments];
	int* l = new int[noOfSegments];
	int* r = new int[noOfSegments];

	for (int i = 0; i < noOfSegments; i++) {
		l[i] = b[i];
		r[i] = b[i + 1];
		temp = ((r[i] - l[i]) < (partitionLen / noOfSegments)) ?
				(r[i] - l[i]) : (partitionLen / noOfSegments);
		s[i] = l[i] + temp + ((i + 1 <= partitionLen % noOfSegments) ? 1 : 0);
	}



	while (true) {
		segContainLmax = -1;
		segContainRmin = 100;
		rmin = 100;
		lmax = -100;

		register int temp_s, temp_l, temp_r, temp_a, temp_a1;
		for (int i = 0; i < noOfSegments; i++) {
			temp_s = s[i];
			temp_l = l[i];
			temp_r = r[i];
			temp_a = a[temp_s - 1];
			temp_a1 = a[temp_s];

			if (temp_s > temp_l && (temp_s - 1) >= temp_l
					&& (temp_s - 1) < temp_r && temp_a >= lmax) {
				lmax = temp_a;
				segContainLmax = i;
			}
			if (temp_s < temp_r && temp_s >= temp_l && temp_a1 <= rmin) {
				rmin = temp_a1;
				segContainRmin = i;
			}
		}

		if (s[segContainLmax] == l[segContainLmax]
				|| s[segContainRmin] == r[segContainRmin] || lmax < rmin
				|| (lmax == rmin && segContainLmax < segContainRmin)) {
			break;
		}

		r[segContainLmax] = s[segContainLmax] - 1;
		l[segContainRmin] = s[segContainRmin] + 1;

		deltaMax = ((r[segContainLmax] - l[segContainLmax]) / 2);
		deltaMin = ((r[segContainRmin] - l[segContainRmin]) / 2);

		delta = (deltaMax < deltaMin) ? deltaMax : deltaMin;

		s[segContainLmax] = r[segContainLmax] - delta;
		s[segContainRmin] = l[segContainRmin] + delta;
		counter1[threadid]++;
	}

	gettimeofday(&endwtime1, NULL);
	seq_time1[threadid] = (double) ((endwtime1.tv_usec - startwtime1.tv_usec)
			/ 1.0e6 + endwtime1.tv_sec - startwtime1.tv_sec);

	cout<< "outside"<<(double) ((endwtime1.tv_usec - startwtime1.tv_usec)
			/ 1.0e6 + endwtime1.tv_sec - startwtime1.tv_sec)
<<endl;
	vector<int> values(s, s + noOfSegments);

	delete[] s;
	delete[] l;
	delete[] r;

	return values;
}

void FrancisPartitioning(vector<int> a, const int noOfSegments,
		const int noOfThreads, const vector<int> segmentsEnds,
		vector<std::pair<int, vector<int> > >& boundaries) {

	int totalArrayLength = a.size();
	int partitionLen[noOfThreads]; // D == rank
	vector<int> b; //last index  + 1
	b.push_back(0);
	numOfComparisonsFransic = 0;

	for (int segmentIndex = 0;
			segmentsEnds.begin() + segmentIndex < segmentsEnds.end();
			segmentIndex++)
		b.push_back(*(segmentsEnds.begin() + segmentIndex) + 1);

	for (int i = 0; i < noOfThreads; i++) {
		counter1.push_back(0);
		seq_time1.push_back(0);
		vector<int> x;
		int a = 0;
		boundaries.push_back(std::make_pair(a, x));
	}

	int x = 2;
	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(noOfThreads / x);
//#pragma omp parallel num_threads(noOfThreads/x)
#pragma omp parallel for schedule(static, x) num_threads(noOfThreads/x)
	for (int i = 0; i < noOfThreads; i++)
	{
		int threadid = i; //omp_get_thread_num();

		partitionLen[threadid] = (totalArrayLength / noOfThreads)
				* (threadid + 1);

		vector<int> res;
		{
			int partitionLen1 = partitionLen[threadid];
			int deltaMax, deltaMin, delta;
			int rmin = 100, lmax = -10;
			int segContainLmax = -1, segContainRmin = 100;
			int temp;

			//int threadid = thread_id;//omp_get_thread_num();
			int* s = new int[noOfSegments];
			int* l = new int[noOfSegments];
			int* r = new int[noOfSegments];

			for (int i = 0; i < noOfSegments; i++) {
				l[i] = b[i];
				r[i] = b[i + 1];
				temp = ((r[i] - l[i]) < (partitionLen1 / noOfSegments)) ?
						(r[i] - l[i]) : (partitionLen1 / noOfSegments);
				s[i] = l[i] + temp + ((i + 1 <= partitionLen1 % noOfSegments) ? 1 : 0);
			}

			while (true) {
				segContainLmax = -1;
				segContainRmin = 100;
				rmin = 100;
				lmax = -100;

				register int temp_s, temp_l, temp_r, temp_a, temp_a1;
				for (int i = 0; i < noOfSegments; i++) {
					temp_s = s[i];
					temp_l = l[i];
					temp_r = r[i];
					temp_a = a[temp_s - 1];
					temp_a1 = a[temp_s];

					if (temp_s > temp_l && (temp_s - 1) >= temp_l
							&& (temp_s - 1) < temp_r && temp_a >= lmax) {
						lmax = temp_a;
						segContainLmax = i;
					}
					if (temp_s < temp_r && temp_s >= temp_l && temp_a1 <= rmin) {
						rmin = temp_a1;
						segContainRmin = i;
					}
				}

				if (s[segContainLmax] == l[segContainLmax]
						|| s[segContainRmin] == r[segContainRmin] || lmax < rmin
						|| (lmax == rmin && segContainLmax < segContainRmin)) {
					break;
				}

				r[segContainLmax] = s[segContainLmax] - 1;
				l[segContainRmin] = s[segContainRmin] + 1;

				deltaMax = ((r[segContainLmax] - l[segContainLmax]) / 2);
				deltaMin = ((r[segContainRmin] - l[segContainRmin]) / 2);

				delta = (deltaMax < deltaMin) ? deltaMax : deltaMin;

				s[segContainLmax] = r[segContainLmax] - delta;
				s[segContainRmin] = l[segContainRmin] + delta;
			}
			vector<int> values(s, s + noOfSegments);

			delete[] s;
			delete[] l;
			delete[] r;

			res = values;
		}
		vector<int> temp;
		for (int i = 0; i < noOfSegments; i++) {
			temp.push_back(res[i] - 1);
		}
		boundaries[threadid] = std::make_pair(threadid, temp);
	}
	sort(boundaries.begin(), boundaries.end());
	//adjust the ending boundaries of each segments to make sure the last element of each segment
	//is included
	for (int i = 0; segmentsEnds.begin() + i < segmentsEnds.end(); i++)
		if ((boundaries.end() - 1)->second[i] != *(segmentsEnds.begin() + i))
			(boundaries.end() - 1)->second[i] = *(segmentsEnds.begin() + i);
}

