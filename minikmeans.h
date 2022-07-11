#ifndef MINIKMEAN_H
#define MINIKMEAN_H

#include "abstractkmeans.h"


/**
Implementation of "Web-scale K-Means Clustering", by D. Schulley

@author:  Cheng-Hao Deng
@modified by Wan-Lei Zhao

@date: 2016-May-2

**/

class MiniKMeans: public AbstractKMeans
{
    unsigned int B0;
    unsigned int T0;
    unsigned int rate;

public:

    MiniKMeans();

    bool   init(const char * srcfn);
    bool   init(float *mat, const int row, const int dim);
    bool   refresh();

    int    clust(const unsigned int clust_num, const char *dstfn, const int verbose);

    bool   config(const char *_seed_, const char *lg_first, const char *crtrn, int verbose);

    void   saveCenters(const char *dstfn, bool append);
    int    fetchCenters(float *centers);

    void   setTO(const unsigned int t);

    virtual ~MiniKMeans();


    static void test();
    static void pctest();
};





#endif // MINKMEAN_H
