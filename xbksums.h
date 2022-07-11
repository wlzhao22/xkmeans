#ifndef XBKSUMS_H
#define XBKSUMS_H

/**implementation of bisecting ksums

@author: Hui Ye, Wan-Lei Zhao
**/

#include "abstractkmeans.h"

#include <vector>
#include <map>
#include <set>

using namespace std;

class XBKSums: public AbstractKMeans
{
private:
    static const int nTrails;
    vector<int> crntClust;
    int *nwlabels, *bstLabels;
    double *Ds;

    static const double ERR0;

public:
    XBKSums();
    virtual ~XBKSums();

    bool        init(const char *srcfn);
    bool        init(float *data, const int row, const int dim);
    bool        config(const char *_seed_, const char *lg_first, const char *crtrn, int verbose);
    int         getL2norms(const unsigned int n, const unsigned int d0);

    int         clust(const unsigned int nclust, const char *dstFn, const int verbose);
    bool        bisect(const unsigned int odlbl, const unsigned int nwlbl, unsigned int &n1, unsigned int &n2);
    bool        optzI1(const unsigned int odlbl, const unsigned int nwlbl);
    bool        optzI2(const unsigned int odlbl, const unsigned int nwlbl);

    double      refineI2(const int nclust, const unsigned int nIter);

    double      crsDst(double Dsc, double *D, const unsigned int x, const unsigned int n);
    double      calDist(double *D,   const unsigned int x, const unsigned int n);
    double      cosDst(double *D,    const unsigned int x, bool _plus_);
    double      calAvgDistort(const unsigned int nclust);
    double      pairwDst(const int nclust);

    bool        refresh();
    void        saveCenters(const char *dstfn, bool append);
    int         fetchCenters(float *centers);

    static void test();

protected:
    bool _INIT_, _REFER_;




};

#endif
