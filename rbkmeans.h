#ifndef RBKMEANS_H
#define RBKMEANS_H

#include <iostream>
#include <vector>
#include <map>

#include "abstractkmeans.h"
#include "nnitem.h"

/**************bisecting one cluster into two****************************
Implementation is according to Ying Zhao and G. Karypis's paper

Step1. Randomly choose two initial centers
            a. Assign points to these two initial centers
            b. Calculate overall intra-cluster similarity score

Step2. perform iteration (incremental iteration)
            a. Randomly choose one vector;
            b. Try to move it to another cluster, to see whether that this
               move increases overall score
            c. Loop 2.a, 2.b for several times to find out vector which gives the best move
    T: end if no significant improvement can be made by one best move

    Criterion: I2, which is according to Zhao Ying, Karypis's paper

    Website: http://glaros.dtc.umn.edu/gkhome/index.php

    The procedure in k-means++ that optimizes the initial center has been incoporated

    paper: k-means++: The Advantages of Careful Seeding

*************************************************************************/

using namespace std;

enum OPT_FUN{ I1 = 1, I2, I3};

class RBKMeans: public AbstractKMeans {

    private:
        float   *D1, *D2, *tmpC1, *tmpC2, *tmpD1, *tmpD2;
        bool     LG_FIRST;
        OPT_FUN  myoptz;
        static const unsigned int NTRAILS;

        double (RBKMeans::*crtrn_func)(const float *D1, const float *D2, const int dim, const int n1, const int n2, double &E1, double &E2);
        bool   (RBKMeans::*optim_func)(const int clabel, double &E1, double &E2,
                                       int &n1, int &n2, const int nwlabel);

        static const float Err0, EPS;

    protected:

        bool   incrOptz(const int clabel, double &E1, double &E2, int &n1, int &n2, const int nwlbl);
        double I1(const float *D1, const float *D2, const int dim, const int n1, const int n2, double &E1, double &E2);
        double I2(const float *D1, const float *D2, const int dim, double &E1, double &E2);
        double getI2(const double *Ds, const int dim, const unsigned int clust_num);

    public:

        RBKMeans();
        bool  init(const char *srcfn);
        bool  init(float *data, const int row, const int dim);
        bool  refresh();

        bool  config(const char *lg_first, const char *distfunc, const char *crtrn, int verbose);
        int   clust(const unsigned int clust_num, const char *dstfn, const int verbose);

        static double  l2_norm(const float *vect1, const unsigned int dim);
        float  cos(const float v1[], const unsigned int idx1, const float v2[],
                   const unsigned int idx2, const unsigned int d0);
        float  l2(const float v1[], const unsigned int idx1, const float v2[],
                        const unsigned int idx2, const unsigned int d0);

        bool   randTwo(const int bound,   unsigned int &r1, unsigned int &r2);
        bool   randTwin(unsigned int &r1, unsigned int &r2, vector<int> &vect);
        bool   KppTwin(unsigned int &r1,  unsigned int &r2, vector<int> &vect);
        static unsigned int randOne(const int bound);

        double fvec_norm(const float *v, const long n, const double norm);
        bool   fvec_scale(float *v, const long n, const double sc);
        bool   cpyCenters(double *dstArray, float *srcArray, const unsigned int ndim);
        void   saveCenters(const char *dstfn, bool append);
        int    fetchCenters(float *centers);

        static void test();
        virtual ~RBKMeans();
};

#endif
