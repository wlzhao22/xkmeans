#ifndef PQMATH_H
#define PQMATH_H

/**
Distribute tasks to different classes, that are
responsible. It is designed for the sake of clarity of
the code.

@author: Wanlei Zhao
@email:  stonescx@gmail.com
@date:   Mar-28-2015

**/

#include "sparsematrix.h"

class PQMath
{
private:
    static const float smallVal0;
public:
    static double cos(const float v1[], const unsigned int s1,
                      const float v2[], const unsigned int s2, const unsigned int d0);

    static float  l2f(const float v1[], const unsigned int s1,
                      const float v2[], const unsigned int s2,
                      const unsigned int dim0);
    static float  l2f(const SparseMatrix &v1, const unsigned int s1,
                      const SparseMatrix &v2, const unsigned int s2,
                      const unsigned int dim0);
    static double l2d(const double v1[], const unsigned int s1,
                      const float v2[],  const unsigned int s2,
                      const unsigned int dim0);
    static double l2d(const double v1[], const unsigned int s1,
                      const double v2[], const unsigned int s2,
                      const unsigned int dim0);

    static double l1d(const double v1[], const unsigned int s1,
                     const double v2[], const unsigned int s2, const unsigned int d0);

    static double l2d(const double v1[], const unsigned int s1,
                      const SparseMatrix &b,  const unsigned int s2,
                      const unsigned int dim0);
    static double getI2(const double *Ds, const unsigned int k, const unsigned int dim, const int *tmpns);
    static double getE2(const double *Ds, const unsigned int k, const unsigned int dim, const int *tmpns);
    static double getI1(const double *Ds, const unsigned int k, const unsigned int dim, const int *tmpns);
    static double getE1(const double *Ds, const unsigned int k, const unsigned int dim, const int *tmpns);

    static bool   l2_norm(float *vect1,  const unsigned int dim);
    static bool   l2_norm(double *vect1, const unsigned int dim);
    static float  l2norm(const float *v1, const unsigned long s1,
                    const float *v2, const unsigned long s2, const unsigned int d0);
    static bool   is_l2Norm(const float *vect1, const unsigned int dim);
    static double fvec_norm(const float *v, const long n,  const double norm);
    static double dvec_norm(const double *v, const long n, const double norm);

    static double innerProduct(const double *v1, const double *v2, const long n);
    static double innerProduct(const SparseMatrix &v1,unsigned int s1, const SparseMatrix &v2, unsigned int s2);

    static void   normVects(float *vects, const unsigned int d0,
                            const unsigned int n0, float *lens);
    static void   normVects(float *vects, const unsigned int d0,
                            const unsigned int n0);
    static void   sqrt_SIFT(float *vects, const unsigned int d0, const unsigned int numb);
    static bool   isZeros(const float *feat, const  unsigned int dim);
    static int    sign(const float val);
    static void   powerLaw(float *vect, const float alpha, const unsigned int d0);
    ///static double fvec_norm(const float *v, const long n, const double norm);
    static bool   fvec_scale(float *v,  const long n, const double sc);
    static bool   dvec_scale(double *v, const long n, const double sc);

public:
    PQMath() {}
    virtual ~PQMath() {}
};

#endif
