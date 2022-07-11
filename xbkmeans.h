#ifndef XBKMEANS_H
#define XBKMEANS_H

#include "abstractkmeans.h"

#include <vector>
using namespace std;

class XBKMeans: public AbstractKMeans
{
private:
    static const int NTRAILS;
    static const int NIter0;
    static const int NRfn;
    static const float Err0;

    double *D1, *D2;
    double *tmpD1, *tmpD2;
    double *tmpC1, *tmpC2;
    double *C1, *C2;
    double *innDcts;
    bool LG_FIRST;
   // OPTZ myoptz;
    vector<int> crntClust;
    int    *nwlabels;

    double (XBKMeans::*optFunc)(const double *D1, const double *D2,
                                const unsigned int dim,   const unsigned int tmpn1,
                                const unsigned int tmpn2);
    bool   (XBKMeans::*sedFunc)(unsigned int &r1, unsigned int &r2, vector<int> &vect);

public:
    XBKMeans();
    virtual ~XBKMeans();

    bool          init(const char *srcfn);
    bool          init(float *data, const int row, const int dim);
    bool          config(const char *_seed_, const char *lg_first, const char *crtrn, int verbose);
    int           clust(const unsigned int clust_num, const char *dstfn, const int verbose);
    bool          incrOptzI2(const int clabel, double &E1, double &E2, int &n1, int &n2, const int nwlbl);
    bool          incrOptzI4(const int clabel, double &E1, double &E2, int &n1, int &n2, const int nwlbl);
    bool          incrOptzE1(const int clabel, double &E1, double &E2, int &n1, int &n2, const int nwlbl);
    bool          incrOptzE2(const int clabel, double &E1, double &E2, int &n1, int &n2, const int nwlbl);
    bool          incrOptzI1(const int clabel, double &E1, double &E2, int &n1, int &n2, const int nwlbl);
    bool          incrOptzI2s(const int clabel, double &E1, double &E2, int &n1, int &n2, const int nwlbl);
    bool          incrOptzE1s(const int clabel, double &E1, double &E2, int &n1, int &n2, const int nwlbl);
    bool          incrOptzE2s(const int clabel, double &E1, double &E2, int &n1, int &n2, const int nwlbl);
    bool          incrOptzI4s(const int clabel, double &E1, double &E2, int &n1, int &n2, const int nwlbl);
    bool          incrOptzI1s(const int clabel, double &E1, double &E2, int &n1, int &n2, const int nwlbl);
    bool          tkmOptz(const int clabel,  double &E1, double &E2, int &n1, int &n2, const int nwlbl);

    /** random k seeds with the way of k-means++ **/
    bool          kppTwin(unsigned int &r1, unsigned int &r2, vector<int> &vect);
    bool          rndTwin(unsigned int &r1, unsigned int &r2, vector<int> &vect);

    double        I3Func(const double *D1, const double *D2,
                         const unsigned int dim,   const unsigned int tmpn1,
                         const unsigned int tmpn2);

    double        I1Func(const double *D1, const double *D2,
                         const unsigned int dim,   const unsigned int tmpn1,
                         const unsigned int tmpn2);

    double        I2Func(const double *D1, const double *D2,
                         const unsigned int dim,   const unsigned int tmpn1,
                         const unsigned int tmpn2);

    double        E1Func(const double *D1, const double *D2,
                         const unsigned int dim,   const unsigned int tmpn1,
                         const unsigned int tmpn2);

    double        E2Func(const double *D1, const double *D2,
                         const unsigned int dim,   const unsigned int tmpn1,
                         const unsigned int tmpn2);

    void          saveCenters(const char *dstfn, bool append);
    int           fetchCenters(float *centers);
    bool          refresh();
    static void   test();
    static void   hyeTest();
};

#endif
