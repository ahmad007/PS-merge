#include "splitMerge.h"
#include "IndexInfoCache.h"
#include <iostream>
#include <vector>
#include <omp.h>
#include <cmath>
#include <omp.h>
#include <stack>
#include <climits>
#include <sys/time.h>


using namespace std;
ofstream file;

long *threads;
int shiftTime = 0;
unsigned long noOfSwaps = 0;
int counter = 0;
int mergeThreads = 1;
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
void resetNoOfSwaps(){
    noOfSwaps = 0;
}

unsigned long getNumOfComparisonsSplitMerge(){
	return comps;
}

void resetNumOfComparisonsSplitMerge(){
	comps = 0;
}


unsigned long getNumOfSwapsOfSplitMerge() {
    if (counter > 0)
        cout << "********swaps multiple*******" << counter << endl;
    return noOfSwaps;
}

void rightShiftCircular(vector<int>& arr, int shiftStartIndex, int shiftEndIndex, int shiftValue) {
    const clock_t begin_time = clock();
    int thread_id = omp_get_thread_num();
    int temp = 0;
    for (int i = shiftStartIndex, j = shiftEndIndex; --j > i; i++) {
        temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
//#pragma omp critical
//        {
//            if (noOfSwaps == ULONG_MAX) {
//                counter++;
//                noOfSwaps = 1;
//            } else
//                noOfSwaps++;
//        }
    }

    for (int i = shiftStartIndex, j = i + shiftValue; --j > i; i++) {
        temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
//#pragma omp critical
//        {
//            if (noOfSwaps == ULONG_MAX) {
//                counter++;
//                noOfSwaps = 1;
//            } else
//                noOfSwaps++;
//        }
    }

    for (int i = shiftStartIndex + shiftValue, j = shiftEndIndex; --j > i; i++) {
        temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
//#pragma omp critical
//        {
//            if (noOfSwaps == ULONG_MAX) {
//                counter++;
//                noOfSwaps = 1;
//            } else
//                noOfSwaps++;
//        }
    }

//#pragma omp critical
//    shiftTime += float( clock() - begin_time);
}

void par_rightShiftCircular1(vector<int>& arr, int shiftStartIndex,
		int shiftEndIndex, int shiftValue) {

	int temp = 0;
	//for (int i = shiftStartIndex, j = shiftEndIndex; --j > i; i++) {
	int j = shiftEndIndex;
	int endIndex =
			((shiftEndIndex - shiftStartIndex - 1) % 2 == 0) ?
					((shiftEndIndex - shiftStartIndex - 1) / 2) - 1 :
					((shiftEndIndex - shiftStartIndex - 1) / 2);

	int threadSize =
			((endIndex / mergeThreads) == 0) ? 1 : endIndex / mergeThreads;

#pragma omp parallel for private( j,temp) schedule(static ,threadSize) num_threads(mergeThreads)
	for (int i = shiftStartIndex; i <= shiftStartIndex + endIndex; i++) {

		j = shiftEndIndex - (i - shiftStartIndex + 1);
		temp = arr[i];
		arr[i] = arr[j];
		arr[j] = temp;

//#pragma omp critical
//		        {
//		         numOfSwaps++;
//		         shuffleSwaps++;
//		        }
	}

	j = shiftStartIndex + shiftValue;
	endIndex =
			((shiftStartIndex + shiftValue - shiftStartIndex - 1) % 2 == 0) ?
					((shiftStartIndex + shiftValue - shiftStartIndex - 1) / 2)
							- 1 :
					((shiftStartIndex + shiftValue - shiftStartIndex - 1) / 2);
	threadSize = ((endIndex / mergeThreads) == 0) ? 1 : endIndex / mergeThreads;

#pragma omp parallel for private( j,temp) schedule(static ,threadSize) num_threads(mergeThreads)
	for (int i = shiftStartIndex; i <= shiftStartIndex + endIndex; i++) {

		j = shiftStartIndex + shiftValue - (i - shiftStartIndex + 1);

		temp = arr[i];
		arr[i] = arr[j];
		arr[j] = temp;
//#pragma omp critical
//		        {
//		         numOfSwaps++;
//		         shuffleSwaps++;
//		        }
	}

	j = shiftEndIndex;
	endIndex =
			((shiftEndIndex - (shiftStartIndex + shiftValue) - 1) % 2 == 0) ?
					((shiftEndIndex - (shiftStartIndex + shiftValue) - 1) / 2)
							- 1 :
					((shiftEndIndex - (shiftStartIndex + shiftValue) - 1) / 2);
	threadSize = ((endIndex / mergeThreads) == 0) ? 1 : endIndex / mergeThreads;

//		//#pragma omp parallel for private( j,temp) schedule(static ,endIndex/numOfThreads) num_threads(numOfThreads)
//		for (int i = shiftStartIndex + shiftValue; i <= shiftStartIndex + shiftValue + endIndex;
//				i++) {
//
//			j = shiftEndIndex -(i-shiftStartIndex + shiftValue+1);

#pragma omp parallel for private( j,temp) schedule(static ,threadSize) num_threads(mergeThreads)
	for (int i = shiftStartIndex + shiftValue;
			i <= shiftStartIndex + shiftValue + endIndex; i++) {
		j = shiftEndIndex - (i - (shiftStartIndex + shiftValue) + 1);
		temp = arr[i];
		arr[i] = arr[j];
		arr[j] = temp;
//#pragma omp critical
//		        {
//		         numOfSwaps++;
//		         shuffleSwaps++;
//		        }
	}

	//#pragma omp critical
	//    shiftTime += float( clock() - begin_time);
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

//        #pragma omp critical
//        		{
//        			comps++;
//        		}

        if (f1 >= f2 || f2 >= last)
            continue;
        l = f1;
        r = f2;
        ldash = f2;
        rdash = last;
        m = 0;
        mdash = 0;

        while (!(l >= r && ldash >= rdash)) {
//#pragma omp critical
//			{
//				comps += 4;
//			}

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
        //ahmad: this method is sequentil and used in coarse grained parallel binary merge tree
        rightShiftCircular(a, r, ldash, (ldash - f2));
        //ahmad: this method is parallel and used in fine-grained parallel binary merge tree
        //par_rightShiftCircular1(a, r, ldash, (ldash - f2));

        param pFirst(f1, r, r + rdash - f2);
        param pSecond(l + ldash - f2, ldash, last);
        parameters.push(pSecond);
        parameters.push(pFirst);
        // SplitMerge(a, f1, r, r + rdash - f2);
        //SplitMerge(a, l + ldash - f2, ldash, last);
    }
}


void print(vector<int> a, int n) {
    for (int i = 0; i < n; i++)
        cout << a[i] << ", ";
    cout << endl;
}


void splitmerge2Portions(vector<int>& a, vector<std::pair<int, int> > limits,int w) {
    int thread_id = w;//omp_get_thread_num();

    //#pragma omp critical
    //    cout << thread_id << " ok" << endl;
    int size = (limits.begin() + 2 * thread_id)->second - (limits.begin() + 2 * thread_id)->first + 1;
    size += (limits.begin() + 2 * thread_id + 1)->second - (limits.begin() + 2 * thread_id + 1)->first + 1;

//    struct timeval startwtime1, endwtime1;
//    		gettimeofday(&startwtime1, NULL);
//
//    			cout<<"f1 "<< (limits.begin() + 2 * thread_id)->first<<" f2: "<<(limits.begin() + 2 * thread_id + 1)->first<<endl;

    iterativeSplitMerge(a, (limits.begin() + 2 * thread_id)->first, (limits.begin() + 2 * thread_id + 1)->first, (limits.begin() + 2 * thread_id)->first + size);
  //  gettimeofday(&endwtime1, NULL);
//    			cout<< (double) ((endwtime1.tv_usec - startwtime1.tv_usec) / 1.0e6	+ endwtime1.tv_sec - startwtime1.tv_sec)<<endl;

    //SplitMerge(a, (limits.begin() + 2 * thread_id)->first, (limits.begin() + 2 * thread_id + 1)->first, (limits.begin() + 2 * thread_id)->first + size);
}


void multiWaySplitMerge(vector<int>& a, vector<int > segmentsEnds, int usedThreads, ofstream* myfile) {
    // cout << "\nTournament split merging started...." << endl;    
    const clock_t begin_time = clock();

    mergeThreads = usedThreads;

    noOfSwaps = 0;
    vector<std::pair<int, int> > limits(segmentsEnds.size()); //actual index
    int size = a.size();

    int start = 0;
    for (int i = 0; segmentsEnds.begin() + i < segmentsEnds.end(); i++) {
        limits.push_back(std::make_pair(start, *(segmentsEnds.begin() + i)));
        start = *(segmentsEnds.begin() + i) + 1;
    }


    //    for (int i = 0; limits.begin() + i < limits.end(); i++)
    //        cout << "(" << (limits.begin() + i)->first << "," <<
    //            (limits.begin() + i)->second << ")" << endl;


    int noOfThreads, tempStart, tempEnd;
    typedef std::vector<std::pair<int, int> > my_vector;
    my_vector::iterator it;
    int phases = std::ceil(std::log(limits.size()) / std::log(2));

    noOfThreads = limits.size() / 2;


    threads = new long[noOfThreads];


    //    cout << std::log(limits.size()) << endl;
    //    cout << std::ceil(std::log(limits.size())) << endl;
    //    cout << phases << endl << endl;


    for (int i = 0; i < phases; i++) {
        noOfThreads = limits.size() / 2;

        for (int k = 0; k < noOfThreads; k++) {
            threads[k] = 0;
        }


        const clock_t begin_time1 = clock();

//ahmad: uncomment this line to turn to coarse grain parallel Binary merge tree
#pragma omp parallel for schedule(static,1) num_threads(usedThreads)
        for(int w = 0 ; w < noOfThreads; w++)
//not important line
        	//#pragma omp parallel num_threads(noOfThreads)
        {
            splitmerge2Portions(a, limits,w);
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
            if (threads [k] > swaps)
                swaps = threads[k];
        }
        //*myfile << "max swaps phase " << i << ": " << swaps << "\n";
    }
  //  *myfile << "swaps: " << getNumOfSwaps() << "\n\n";


    //cout << "Tournament split merging ended " << endl;
    //std::cout << "Elapsed time:  " << float( clock() - begin_time) / CLOCKS_PER_SEC << " seconds\n";
    //std::cout << "Elapsed time:  " << float( clock() - begin_time) / 1000l << " milli-seconds\n\n";

    // cout << shiftTime / 1000l << endl;
    //    for (int i = 0; i < size; i++)
    //        cout << a[i] << ",";
    //    cout << endl;

}

bool check(vector<int> a) {
    for (int i = 0; i < a.size() - 1; i++)
        if (a[i] > a[i + 1]) {
            cout << "error at: a[" << i << "]=" << a[i] <<
                    "[" << i + 1 << "]=" << a[i + 1] << endl;
            return false;
        }
    cout << "check done!!!!!!!!!\n";
    return true;
}
