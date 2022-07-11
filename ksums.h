#ifndef KSUMS_H
#define KSUMS_H

/**
implementation of generic k-means
it is designed to work in the generic space

@author: Wan-Lei Zhao, Hui Ye
**/

#include <vector>
#include <set>
#include <map>
#include <random>
#include <time.h>

#include "abstractkmeans.h"
#include "nnitem.h"

using namespace std;

class KSums: public AbstractKMeans
{
    public:
        KSums();
        virtual ~KSums();

        vector<vector<MiniNN> > knnGraph;  ///keeps i's neighbors
        vector<NbHood> nbGraph;
        map<unsigned int, double> topcSums;
        map<unsigned int, unsigned int> topcNums;
        map<unsigned int, set<unsigned int> > clusters;

        static const int k0, BD0, nIter1, cNum, nIter2, nTrail;
        double errSum, dstBound;

        clock_t bg, ed;
        int verbose;

    public:
        void getRndSds(std::mt19937 &rng, const unsigned int k, unsigned int *addr, const unsigned int N);

    public:
        bool   config(const char *_seed_, const char *lg_first, const char *crtrn, int verbose);
        void   saveCenters(const char *dstFn, bool append);
        int    fetchCenters(float *centers);
        bool   refresh();

    public:
        bool   init(const char *srcFn);
        bool   init(float *mat, const int row, const int dim);

        int    initNon(const unsigned int nclust);
        int    clust(unsigned int clust_num, const char *dstFn, const int verbose);

        int    optzIn(const unsigned int nIter, const int k0);
        int    optzI2(const unsigned int nIter, const int k0);
        int    optzIs(const unsigned int nIter, const int k0);
        int    optzIx(const unsigned int nIter, const int k0);
        int    numInKnn(const unsigned int i, const int clabel);
        int    overlap(const unsigned int i, const unsigned int j);
        int    sumShares(const unsigned int i, set<unsigned int> &clust);
        int    evalNeib(const unsigned int k0);
        double calAvgDistort(const unsigned int k0);

        double calDist(double *D, const unsigned int x, const unsigned int n);
        double cosDst(double *D, const unsigned int x, bool _plus_);

        int    loadkNNGraph(const char *graphFn);
        int    saveKNNGraph(const char *dstFn);
        int    initKnnGraph(const unsigned int k);
        int    updateLst(const unsigned int id, const unsigned int nb, float dst);
        int    appndLst(const unsigned int id,  const unsigned int nb, float dst);
        int    appRvNbs(unsigned int n);
        unsigned nnDescent();
        int    getNBGraph(const unsigned int smplNum);
        int    saveRNNGraph(const char *dstFn);

        void    saveDistort(const char *dstFn);

    protected:
        int *bstLabels;
        double *dstSum; ///keep the sum of distance in each cluster
        unsigned int nCmps, tCmps; //number of elements and dimension of data
        unsigned int nnum, nIter; ///number of nearest elements for each data
        bool _INIT_, _REFER_;


    public:
        static void test();
        static void wltest();


};

#endif // GCKMEANS_H
