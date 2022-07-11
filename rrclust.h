#ifndef RRCLUST_H
#define RRCLUST_H

#include "abstractkmeans.h"

#include <map>

/**

Round-robin clustering,
However, its performance is only close to xtkm
while it is slower than xbk

Overall, it is in between xtk and xbk in terms of both speed
and clustering quality
**/

class RRClust : public AbstractKMeans
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
    OPTZ myoptz;
    double (RRClust::*optFunc)(const double *Ds, const unsigned int k, const unsigned int dim,
                                const int *tmpns, double *dsts);

public:
    RRClust();
    virtual ~RRClust();

    bool   init(const char *srcfn);
    bool   init(float *data, const int row, const int dim);
    bool   config(const char *_seed_, const char *lg_first, const char *crtrn, int verbose);
    int    clust(const unsigned int clust_num, const char *dstfn, const int verbose);
    bool   incrOptzI2(double *Es, const unsigned int clustNum);
    double initialCenters(const float *data, double *Ds, int *Nss, const unsigned int col, const unsigned int row, int *labels);
    bool   BiMove(vector<int> &cluster1, vector<int > &cluster2, double *D1, double *D2);
    bool   BiOptz(int clbl1, int clbl2, double *D1, double *D2, int &n1, int &n2);
    bool   allocMem(const unsigned int clustNum);
    bool   deallocMem();
    bool   refresh();

    double I2Func(const double *Ds, const unsigned int k, const unsigned int dim,
                  const int *tmpns, double *dsts);

    void    saveCenters(const char *dstfn, bool append);
    int     fetchCenters(float *centers);

    static void   test();

};






#endif // RR_H
