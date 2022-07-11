#ifndef KMeans_H
#define KMeans_H

#include "abstractkmeans.h"
#include "nnitem.h"
#include <vector>

using namespace std;

/**
*
*
*
The implementation of Traditional k-means and k-means++
*
@author:  Wan-Lei Zhao
@date:    Mar.05-2016 - Mar./30/2016
*
*
*
*
**/

class KMeans: public AbstractKMeans
{
    private:
        static const unsigned int M_ITER, NTRAILS, NRfn;
        static const float EPS, Err0;
        double *tmparrayD, *bstArrayD, *Cs;
        int    *bstLabels;

    public:
        KMeans();
        virtual ~KMeans();
        bool   init(const char *srcfn);
        bool   init(float *mat, const int row, const int dim);
        bool   refresh();

        bool   (KMeans::*sedFunc)(const int k, int rseeds[], const int bound);

        bool   config(const char *_seed_, const char *lg_first, const char *crtrn, int verbose);

        double *arrayDs;
        int    getL2norms(const unsigned int n, const unsigned int d0);
        double pairwDst(const int nclust);
        int    clust(const unsigned int clust_num, const char *dstfn, const int verbose);
        bool   rndSeeds(const int k, int rseeds[], const int bound);
        bool   kppSeeds(const int k, int rseeds[], const int bound);
        bool   lucSeeds(const int k, int rseeds[], const int bound);
        void   saveCenters(const char *dstfn, bool append);
        int    fetchCenters(float *centers);

        static void test();
};

#endif

