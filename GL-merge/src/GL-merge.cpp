//============================================================================
// Name        : GL-merge.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
using namespace std;

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
#include "Guan_langston.h"

using namespace std;

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
 * @return
 */
int main() {
	vector<vector<int> > array;

	int noOfSegments = 2; //k value in the k-way merge
	int array_size = std::pow(2.0, 10);
	int numberOfthreads = 2;
	vector<vector<int> > segmentsEnds;
	cout << "\nGenerating data is started...." << flush;
	generateDataForGuang(array, array_size, noOfSegments, segmentsEnds);
	cout << "Done!" << endl;

	cout << "Tournament GuangLangston merging is started...." << flush;
	guanLangstonMultiwayMerging(array, segmentsEnds, numberOfthreads);
	cout << "..merging is done!" << endl;
	cout << "Array size: " << array_size << ", # of segments: "
			<< segmentsEnds.size() << ", # of threads: " << numberOfthreads
			<< ", ";

#if IS_DEBUG == 0
	cout << "wall clock time = " << accumTime << " seconds." << endl;
#endif
#if IS_DEBUG == 1
	cout << "# of swaps: " << noOfGuangSwaps << endl;

#endif
	if (checkArrayOrder(array[0])) {
		cout << "The merged array is checked successfully!" << endl;
	}

	cout << "\nAll tests are done!" << endl;
	return 0;
}

