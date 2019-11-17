/* 
 * File:   IndexInfoCache.h
 * Author: ahmad
 *
 * Created on December 30, 2013, 2:21 PM
 */

#ifndef INDEXINFOCACHE_H
#define	INDEXINFOCACHE_H

struct IndexCache {
public:
    int size;
    int * virtualStart;
    int * virtualEnd;
    int * actualStart;
    int * actualEnd;

    void setSize(int cacheSzie) {
        size = cacheSzie;
        virtualStart = new int [cacheSzie];
        virtualEnd = new int [cacheSzie];
        actualStart = new int [cacheSzie];
        actualEnd = new int [cacheSzie];
        for (int i = 0; i < cacheSzie; i++)
            virtualStart[i] = virtualEnd [i] = actualStart[i] = actualEnd[i] = 0;
    }

    void finalize() {
        delete [] virtualStart;
        delete [] virtualEnd;
        delete [] actualStart;
        delete [] actualEnd;
    }

};

#endif	/* INDEXINFOCACHE_H */

