#ifndef ONLINEKMEANS_H
#define ONLINEKMEANS_H

#include "abstractkmeans.h"

/***
an implementation of sequential k-means

@author: Wan-Lei Zhao
@date:   Sep-17-2019double SeqKMeans::l2nDst(double *D, const unsigned int x, const unsigned int n)
**/

class SeqKMeans:public AbstractKMeans
{
    public :

    SeqKMeans();

    bool   init(const char * srcfn);
    bool   init(float *mat, const int row, const int dim);
    bool   refresh();

    static const int nTrails;
    double *Ds;
    int    record_num;
    int    getL2norms(const unsigned int n, const unsigned int d0);
    double l2Dst(const unsigned int c, const unsigned int s, const unsigned int n);
    double pairwDst(const int nclust, int seq[], int locount);
    double crsDst(double Dsc, double *D, const unsigned int x, const unsigned int n);
    double l2nDst(double *D, const unsigned int x, const unsigned int n);
    double avgDist(const unsigned int nclust, int seq[], int locount, int lp);
    int    opt1trail(unsigned int clust_num);
    int    clust(const unsigned int clust_num, const char *dstfn, const int verbose);

    bool   config(const char *_seed_, const char *lg_first, const char *crtrn, int verbose);

    void   saveCenters(const char *dstfn, bool append);
    int    fetchCenters(float *centers);

    ~SeqKMeans();

    static void test();



};


#endif // ONLINEKMEANS_H
