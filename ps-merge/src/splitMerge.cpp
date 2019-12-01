#include "splitMerge.h"
#include <iostream>
#include <vector>
#include <omp.h>
#include <cmath>
#include <omp.h>
#include <stack>
#include <climits>

using namespace std;
ofstream file;

long *threads;
int shiftTime = 0;
unsigned long noOfSwaps = 0;
int counter = 0;
unsigned long comps = 0;
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
void resetNoOfSwaps() {
	noOfSwaps = 0;
}

unsigned long getNumOfSwapsOfSplitMerge() {
	if (counter > 0)
		cout << "********swaps multiple*******" << counter << endl;
	return noOfSwaps;
}


unsigned long getNumOfComparisonsSplitMerge(){
	return comps;
}

void resetNumOfComparisonsSplitMerge(){
	comps = 0;
}

void rightShiftCircular(vector<int>& arr, int shiftStartIndex,
		int shiftEndIndex, int shiftValue, double & swaps) {
	// const clock_t begin_time = clock();
	// int thread_id = omp_get_thread_num();
	int temp = 0;
	for (int i = shiftStartIndex, j = shiftEndIndex; --j > i; i++) {
		temp = arr[i];
		arr[i] = arr[j];
		arr[j] = temp;
#if IS_DEBUG == 1
#pragma omp critical
        {
         swaps++;
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
         swaps++;
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
         swaps++;
        }
#endif
	}
}

void iterativeSplitMerge(vector<int>& a, int f1, int f2, int last,
		double & swaps) {
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

		if (f1 >= f2 || f2 >= last)
			continue;

		l = f1;
		r = f2;
		ldash = f2;
		rdash = last;
		m = 0;
		mdash = 0;

		while (!(l >= r && ldash >= rdash)) {
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
		rightShiftCircular(a, r, ldash, (ldash - f2), swaps);
		param pFirst(f1, r, r + rdash - f2);
		param pSecond(l + ldash - f2, ldash, last);
		parameters.push(pSecond);
		parameters.push(pFirst);
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
		// cout<<f1 << ", " << f2 << ", " << last <<endl;
		if (f1 >= f2 || f2 >= last)
			continue;
		l = f1;
		r = f2;
		ldash = f2;
		rdash = last;
		m = 0;
		mdash = 0;

		while (!(l >= r && ldash >= rdash)) {
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
		double x = 0;
		rightShiftCircular(a, r, ldash, (ldash - f2), x);
		param pFirst(f1, r, r + rdash - f2);
		param pSecond(l + ldash - f2, ldash, last);
		parameters.push(pSecond);
		parameters.push(pFirst);
	}
}

void print(vector<int> a, int n) {
	for (int i = 0; i < n; i++)
		cout << a[i] << ", ";
	cout << endl;
}

void testSplitMerge() {
	//    int x[] = {1, 2, 3, 4, 5, 1, 2, 3, 4};
	//    SplitMerge(x, 0, 5, 9);
	//    print(x, 9);
}

void multiWaySplitMerge(vector<int>& a, vector<int> segmentsEnds,
		ofstream* myfile) {
	// cout << "\nTournament split merging started...." << endl;
	const clock_t begin_time = clock();

	noOfSwaps = 0;
	vector<std::pair<int, int> > limits(segmentsEnds.size()); //actual index
	int size = a.size();

	int start = 0;
	for (int i = 0; segmentsEnds.begin() + i < segmentsEnds.end(); i++) {
		limits.push_back(std::make_pair(start, *(segmentsEnds.begin() + i)));
		start = *(segmentsEnds.begin() + i) + 1;
	}

		int noOfThreads, tempStart, tempEnd;
	typedef std::vector<std::pair<int, int> > my_vector;
	my_vector::iterator it;
	int phases = std::ceil(std::log(limits.size()) / std::log(2));

	noOfThreads = limits.size() / 2;

	threads = new long[noOfThreads];

		for (int i = 0; i < phases; i++) {
		noOfThreads = limits.size() / 2;

		for (int k = 0; k < noOfThreads; k++) {
			threads[k] = 0;
		}

		const clock_t begin_time1 = clock();

#pragma omp parallel num_threads(noOfThreads)
		{
			//  splitmerge2Portions(a, limits);
		}

		for (int threadId = 0; threadId < noOfThreads; threadId++) {
			it = limits.begin();
			tempStart = (it + threadId)->first;
			tempEnd = (it + (threadId + 1))->second;
			limits.erase(it + threadId, it + (threadId + 2));
			limits.insert(it + threadId, std::make_pair(tempStart, tempEnd));
		}

		long swaps = -1;
		for (int k = 0; k < noOfThreads; k++) {
			if (threads[k] > swaps)
				swaps = threads[k];
		}
	}
}

bool check(vector<int> a) {
	for (int i = 0; i < a.size() - 1; i++)
		if (a[i] > a[i + 1]) {
			cout << "error at: a[" << i << "]=" << a[i] << "[" << i + 1 << "]="
					<< a[i + 1] << endl;
			return false;
		}
	return true;
}
