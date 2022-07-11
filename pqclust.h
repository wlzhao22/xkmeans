#ifndef PQCLUST_H
#define PQCLUST_H

#include "abstractkmeans.h"



class PQClust : public AbstractKMeans
{

private:
    static const unsigned int M_ITER, NTRAILS, NRfn;
    static const float EPS, Err0;
    static const char * pqfn;
    double *Ds, *bstArrayD, *Cs, *tmpCs;
    int    *bstLabels;

    unsigned int *data;
    unsigned pql, pqm, pqn, cmps;
    double *pqdis;
    map<unsigned, float *> pqmap;

public:
    PQClust();
    virtual ~PQClust();
    bool   init(const char *srcfn);
    bool   init(unsigned *mat, const int row, const int dim);
    bool   init(float *mat, const int row, const int dim);
    bool   refresh();

    bool   config(const char *_seed_, const char *lg_first, const char *crtrn, int verbose);

    int    tkm(const unsigned int clust_num, const char *dstfn, const int verbose);
    int    I2(const unsigned int clust_num, const char *dstfn, const int verbose);
    int    clust(const unsigned int clust_num, const char *dstfn, const int verbose);
    bool   rndSeeds(const int k, int rseeds[], const int bound);
    bool   kppSeeds(const int k, int rseeds[], const int bound);
    void   saveCenters(const char *dstfn, bool append);
    int    fetchCenters(float *centers);

    double PQDist(unsigned * pq1, unsigned *pq2);

    static void test();
};

#endif // PQCLUST_H
