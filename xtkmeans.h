#ifndef XTKMEANS_H
#define XTKMEANS_H

#include "abstractkmeans.h"

//#define _WATCH_

/**
*implementation of external k-means, shorten as XT-K-means
*
*
*@quthor: Wanlei Zhao
*@date:   Apr.-4 2016
*
*
**/

class XTKMeans: public AbstractKMeans
{
protected:
    static const int NTRAILS;
    static const int NIter0;
    static const float Err0;
    static const int NRef;
    static const int Round;
    static const int Pubnum;

    double *Ds;
    double *Cs;
    double *tmpCs;
    bool LG_FIRST;
    #ifdef _WATCH_
    int   hits[100000][50];
    int  *hists;
    #endif // _WATCH_
    OPTZ myoptz;
    double (XTKMeans::*optFunc)(const double *Ds, const unsigned int k, const unsigned int dim,
                                const int *tmpns, double *dsts);

public:
    XTKMeans();
    virtual ~XTKMeans();

    bool   init(const char *srcfn);
    bool   init(float *data, const int row, const int dim);
    bool   config(const char *_seed_, const char *lg_first, const char *crtrn, int verbose);
    int    clust(const unsigned int clust_num, const char *dstfn, const int verbose);
    bool   incrOptzI2(double *Es, const unsigned int clustNum);
    bool   incrOptzI2s(double *Es, const unsigned int clustNum);
    bool   incrOptzI4(double *Es, const unsigned int clustNum);
    bool   incrOptzI4s(double *Es, const unsigned int clustNum);
    bool   augOptz(double *Es, map<int, const char*> datafile, const unsigned int clustNum);
    bool   augOptz(double *Es, const unsigned int clustNum);
    bool   tkmOptz(const int clabel,  double &E1, double &E2, int &n1, int &n2, const int nwlbl);
    double initialCenters(const float *data, double *Ds, int *Nss, const unsigned int col, const unsigned int row, int *labels);
    double extraRefine(const float *data, double *Ds, int *Nss, const unsigned int col, const unsigned int row, int *labels);

    bool   allocMem(const unsigned int clustNum);
    bool   deallocMem();
    bool   refresh();

    bool   rndSeeds(const unsigned int k, int rseeds[], const unsigned int bound);
    bool   kppSeeds(const unsigned int k, int rseeds[], const unsigned int bound);
    double I2Func(const double *Ds, const unsigned int k, const unsigned int dim,
                  const int *tmpns, double *dsts);

    int    getL2norms(const unsigned int n, const unsigned int d0);
    double pairwDst(const int nclust, int* tmpNs, int* mylabels);

    void    saveCenters(const char *dstfn, bool append);
    void    saveHits(const char *dstFn, int cN);
    void    saveHist(const char *dstFn, int cN);
    int     fetchCenters(float *centers);

    static void   test();
};

#endif
