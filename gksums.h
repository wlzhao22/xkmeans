#ifndef GKSUM_H
#define GKSUM_H

#include <vector>
#include <set>
#include <map>
#include <time.h>
#include <random>

#include "abstractkmeans.h"
#include "nnitem.h"

class GKSums: public AbstractKMeans
{
    public:
        GKSums();
        float entrp0;
        virtual ~GKSums();
        vector<vector<MiniNN> > knnGraph;  ///keeps i's neighbors
        vector<NbHood> nbGraph;

        double *distortions1, *tcs1, *distortions2, *tcs2; ///keeps the average distortions in each iteration
        static const int k0, BD0, nIter1, cNum, nIter2, nTrail;
        float *entropies;
        double errSum;
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
        int    numInKnn(const unsigned int i, const int clabel);
        int    evalNeib(const unsigned int k0);
        double calAvgDistort(const unsigned int k0);

        double calDist(double *D, const unsigned int x, const unsigned int n);
        double cosDst(double *D, const unsigned int x, bool _plus_);

        int    loadkNNGraph(const char *graphFn);
        int    saveKNNGraph(const char *dstFn);
        int    initKnnGraph(const unsigned int k);
        int    updateLst(const unsigned int id, const unsigned int nb, float dst);
        unsigned nnDescent();
        int    getNBGraph(const unsigned int smplNum);

        void    saveDistort(const char *dstFn);

    protected:
        int *bstLabels;
        char *checked;
        double *dstSum; ///keep the sum of distance in each cluster
        unsigned int nCmps, tCmps; //number of elements and dimension of data
        unsigned int nnum, nIter; ///number of nearest elements for each data
        bool _INIT_, _REFER_;


    public:
        static void test();
        static void wltest();

};

#endif // GKSUM_H
