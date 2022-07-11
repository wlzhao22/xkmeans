#ifndef LVECTQUANT_H
#define LVECTQUANT_H

#include "abstractkmeans.h"

#include "nnitem.h"
#include <vector>

/**
Implementation of online Learning Vector Quantization
LVQ1 from paper "The Self-Organizing Map", Tuovo Kohonen

author: Wan-Lei Zhao
date:   2016-Sept.-21

**/
class LVectQuant: public AbstractKMeans
{
    private:
        static const unsigned int M_ITER, NTRAILS, NRfn;
        static const float EPS, Err0, alpha0;

        double *tmparrayD, *bstArrayD, *Cs, *Ds;
        double *tmpCs;
        int    *bstLabels;

    public:
        LVectQuant();
        virtual ~LVectQuant();
        bool   init(const char *srcfn);
        bool   init(float *mat, const int row, const int dim);
        bool   refresh();
        bool   allocMem(const unsigned int clustNum);
        bool   deallocMem();
        bool   (LVectQuant::*sedFunc)(const int k, int rseeds[], const int bound);

        bool   config(const char *_seed_, const char *lg_first, const char *crtrn, int verbose);

        int    clust(const unsigned int clust_num, const char *dstfn, const int verbose);
        bool   rndSeeds(const int k, int rseeds[], const int bound);
        void   saveCenters(const char *dstfn, bool append);
        int    fetchCenters(float *centers);
        static void test();
};

#endif
