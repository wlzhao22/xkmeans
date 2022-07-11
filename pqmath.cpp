#include "pqmath.h"

#include <unordered_map>
#include <iostream>
#include <cassert>
#include <cstring>
#include <cmath>
#include <map>

using namespace std;

const float PQMath::smallVal0 = 0.000001f;

double PQMath::cos(const float v1[], const unsigned int s1,
                   const float v2[], const unsigned int s2, const unsigned int d0)
{
    double w1  = 0.0f, w2 = 0.0f, val = 0.0f;
    unsigned int loc1 = d0 * s1;
    unsigned int loc2 = d0 * s2;

    for(unsigned int i = 0; i < d0; i++)
    {
        w1  = w1  + v1[loc1+i]*v1[loc1+i];
        w2  = w2  + v2[loc2+i]*v2[loc2+i];
        val = val + v1[loc1+i]*v2[loc2+i];
    }

    w1 = sqrt(w1);
    w2 = sqrt(w2);

    if(w1 == 0 && w2 == 0)
    {
        return 1;
    }

    if(w1 == 0)
    {
        return fabs(1-w2/2);
    }
    else if(w2 == 0)
    {
        return fabs(1-w1/2);
    }
    else
    {
        val = val/(w1*w2);
        return val;
    }
    return val;
}

float PQMath::l2f(const float v1[], const unsigned int s1,
                  const float v2[], const unsigned int s2, const unsigned int d0)
{
    float w1 = 0.0f, w2 = 0.0f;
    int loc1 = d0 * s1;
    int loc2 = d0 * s2;

    for(unsigned int i = 0; i < d0; i++)
    {
        w1  = v1[loc1+i] - v2[loc2+i];
        w2 += w1*w1;
    }

    ///w2 = sqrt(w2);

    return w2;
}

float PQMath::l2norm(const float *v1, const unsigned long s1,
                    const float *v2, const unsigned long s2, const unsigned int d0)
{
    float res = 0, diff = 0;
    const float *pVect1 = v1 + s1*d0;
    const float *pVect2 = v2 + s2*d0;

    size_t qty4 = d0/4;
    const float* pEnd1 = pVect1 + qty4 * 4;
    const float* pEnd2 = pVect1 + d0;

    ///#pragma omp parallel
    {
    while (pVect1 < pEnd1)
    {
        diff = *pVect1++ - *pVect2++;
        res += diff * diff;
        diff = *pVect1++ - *pVect2++;
        res += diff * diff;
        diff = *pVect1++ - *pVect2++;
        res += diff * diff;
        diff = *pVect1++ - *pVect2++;
        res += diff * diff;
    }
    }

    while (pVect1 < pEnd2)
    {
        diff = *pVect1++ - *pVect2++;
        res += diff * diff;
    }

    return res;
}

float PQMath::l2f(const SparseMatrix &v1, const unsigned int s1,
                  const SparseMatrix &v2, const unsigned int s2,
                  const unsigned int d0)
{
    float w1 = 0.0f, w2 = 0.0f;
    unsigned int i = 0;

    /**/
    unsigned int *mult = NULL;
    unordered_map<unsigned int, unsigned int> s;
    unordered_map<unsigned int, unsigned int>::iterator it;

    if((v1.col[s1+1] - v2.col[s2+1] - v1.col[s1] + v2.col[s2]) > 0)
    {

        mult = new unsigned int[v1.col[s1+1] - v1.col[s1]];
        memset(mult, 0, (v1.col[s1+1] - v1.col[s1])*sizeof(unsigned int));

        for(i = v1.col[s1]; i < v1.col[s1+1]; i++)
        {
            s.insert(make_pair(v1.index[i],i));
        }

        for(i = v2.col[s2]; i < v2.col[s2+1]; i++)
        {
            it = s.find(v2.index[i]);
            if(it != s.end())
            {
                w1  = v1.data[it->second] - v2.data[i];
                mult[it->second - v1.col[s1]] = 1;
            }
            else
                w1 = v2.data[i];
            w2 += w1*w1;
        }
        for(i = v1.col[s1]; i < v1.col[s1+1]; i++)
        {
            if(!mult[i - v1.col[s1]])
            {
                w1 = v1.data[i];
                w2 += w1*w1;
            }
        }
    }
    else
    {
        mult = new unsigned int[v2.col[s2+1] - v2.col[s2]];
        memset(mult, 0, (v2.col[s2+1] - v2.col[s2])*sizeof(unsigned int));

        for(i = v2.col[s2]; i < v2.col[s2+1]; i++)
        {
            s.insert(make_pair(v2.index[i],i));
        }

        for(i = v1.col[s1]; i < v1.col[s1+1]; i++)
        {
            it = s.find(v1.index[i]);
            if(it != s.end())
            {
                w1  = v2.data[it->second] - v1.data[i];
                mult[it->second - v2.col[s2]] = 1;
            }
            else
                w1 = v1.data[i];
            w2 += w1*w1;
        }
        for(i = v2.col[s2]; i < v2.col[s2+1]; i++)
        {
            if(!mult[i - v2.col[s2]])
            {
                w1 = v2.data[i];
                w2 += w1*w1;
            }
        }
    }
    if(mult != NULL)
    {
       delete [] mult;
       mult = NULL;
    }

    w2 = sqrt(w2);
    return w2;
    /**/
    /**
    unsigned int *id = new unsigned[d0];
    for(unsigned int i = 0; i < d0; i++)
    {
        id[i] = -1;
    }

    for(unsigned int i = v1.col[s1]; i < v1.col[s1+1]; i++)
    {
        id[v1.index[i]] = i;
    }

    for(unsigned int i = v2.col[s2]; i < v2.col[s2+1]; i++)
    {
        if(id[v2.index[i]] != -1)
        {
            w1  = v1.data[id[v2.index[i]]] - v2.data[i];
            id[v2.index[i]] = -1;
        }
        else
            w1 = v2.data[i];
        w2 += w1*w1;
    }
    for(unsigned int i = v1.col[s1]; i < v1.col[s1+1]; i++)
    {
        if(id[v1.index[i]] != -1)
        {
            w1 = v1.data[i];
            w2 += w1*w1;
        }
    }
    delete [] id;
    w2 = sqrt(w2);
    return w2;
    /**/

}


double PQMath::l2d(const double v1[], const unsigned int s1,
                   const float v2[], const unsigned int s2, const unsigned int d0)
{
    double w1 = 0.0f, w2 = 0.0f;
    int loc1 = d0*s1;
    int loc2 = d0*s2;
    for(unsigned int i = 0; i < d0; i++)
    {
        w1  = v1[loc1+i] - (double)v2[loc2+i];
        w2 += w1*w1;
    }

    w2 = sqrt(w2);

    return w2;
}


double PQMath::l2d(const double v1[], const unsigned int s1,
                   const double v2[], const unsigned int s2, const unsigned int d0)
{
    double w1 = 0.0f, w2 = 0.0f;
    int loc1 = d0*s1;
    int loc2 = d0*s2;

    for(unsigned int i = 0; i < d0; i++)
    {
        w1  = v1[loc1+i] - v2[loc2+i];
        w2 += w1*w1;
    }

    w2 = sqrt(w2);

    return w2;
}

double PQMath::l1d(const double v1[], const unsigned int s1,
                   const double v2[], const unsigned int s2, const unsigned int d0)
{
    double w1 = 0.0f, w2 = 0.0f;
    int loc1 = d0*s1;
    int loc2 = d0*s2;

    for(unsigned int i = 0; i < d0; i++)
    {
        w1  = abs(v1[loc1+i] - v2[loc2+i]);
        w2 += w1;
    }

    return w2;
}


double PQMath::l2d(const double v1[], const unsigned int s1,
                   const SparseMatrix &b, const unsigned int s2, const unsigned int d0)
{
    double w1 = 0.0f, w2 = 0.0f;
    int loc1 = d0*s1;
    unsigned int *id = new unsigned int[d0];
    memset(id, 0, d0*sizeof(unsigned int));
    for(unsigned int i = b.col[s2]; i < b.col[s2+1]; i++)
    {
        w1 = v1[loc1+b.index[i]] - (double)b.data[i];
        w2 += w1*w1;
        id[b.index[i]] = 1;
    }
    for(unsigned int i = 0; i < d0; i++)
    {
        if(!id[i])
            w2 += v1[i]*v1[i];
    }

    w2 = sqrt(w2);
    //exit(0);
    delete id;
    return w2;
}

double PQMath::getI2(const double *Ds, const unsigned int k, const unsigned int dim, const int *tmpns)
{
    double sumDst = 0, tmpE = 0;
    unsigned int i = 0, j = 0, loc = 0;

    for(j = 0; j < k; j++)
    {
        if(tmpns[j] <= 0)
            continue;

        loc = dim*j;
        tmpE = 0;

        for(i = 0; i < dim; i++)
        {
            tmpE += Ds[loc+i]*Ds[loc+i];
        }
        sumDst += tmpE/tmpns[j];
    }
    return sumDst;
}

double PQMath::getE2(const double *Ds, const unsigned int k, const unsigned int dim, const int *tmpns)
{
    double sumDst = 0, tmpE = 0;
    unsigned int i = 0, j = 0, loc = 0, loc1 = 0;

    for(j = 0; j < k; j++)
    {
        if(tmpns[j] <= 0)
            continue;

        loc1 = dim*j;
        for(unsigned int j2 = j + 1; j2 < k; j2++)
        {
            tmpE = 0;
            if(tmpns[j2] <= 0)
                continue;
            loc = dim*j2;
            for(i = 0; i < dim; i++)
            {
                tmpE += (Ds[loc+i]/tmpns[j2] - Ds[loc1+i]/tmpns[j])*(Ds[loc+i]*tmpns[j] - Ds[loc1+i]*tmpns[j2]);
            }
            sumDst += tmpE;
        }

    }
    return sumDst;
}

double PQMath::getE1(const double *Ds, const unsigned int k, const unsigned int dim, const int *tmpns)
{
    return getE2(Ds, k, dim, tmpns);
}

double PQMath::getI1(const double *Ds, const unsigned int k, const unsigned int dim, const int *tmpns)
{
    return getI2(Ds, k, dim, tmpns);
}


bool PQMath::l2_norm(float *vect1, const unsigned int dim)
{
    float  sz = 0;
    for(unsigned int i = 0; i < dim; i++)
    {
        sz += vect1[i]*vect1[i];
    }

    if(sz == 0.0f)
    {
        return false;
    }
    sz = sqrt(sz);
    for(unsigned int i = 0; i < dim; i++)
    {
        vect1[i] = vect1[i]/sz;
    }

    return true;
}


bool PQMath::l2_norm(double *vect1, const unsigned int dim)
{
    double  sz = 0;
    for(unsigned int i = 0; i < dim; i++)
    {
        sz += vect1[i]*vect1[i];
    }

    if(sz == 0.0f)
    {
        return false;
    }
    sz = sqrt(sz);
    for(unsigned int i = 0; i < dim; i++)
    {
        vect1[i] = vect1[i]/sz;
    }

    return true;
}


bool PQMath::is_l2Norm(const float *vect1, const unsigned int dim)
{
    float  sz = 0;
    for(unsigned int i = 0; i < dim; i++)
    {
        sz += vect1[i]*vect1[i];
    }

    if(sz == 0.0f)
    {
        return true;
    }

    if(sz < 1.0001 && sz > 0.9998)
    {
       return true;
    }else{
       return false;
    }
}

void PQMath::powerLaw(float *vect, const float alpha, const unsigned int d0)
{
    unsigned int i = 0;
    float tmpVal;
    for(i = 0; i < d0; i++)
    {
        tmpVal = fabs(vect[i]);
        vect[i] = PQMath::sign(vect[i])*pow(tmpVal, alpha);
    }
}

void PQMath::normVects(float *vects, const unsigned int d0,
                       const unsigned int n0, float *lens)
{
    assert(vects);
    assert(lens);
    unsigned int i = 0, j = 0, loc = 0;
    float len = 0;
    for(i = 0; i < n0; i++)
    {
        len = 0;
        for(j = 0; j < d0; j++)
        {
            len += vects[loc+j]*vects[loc+j];
        }
        len     = sqrt(len);
        lens[i] = len;
        if(len > PQMath::smallVal0)
        {
            for(j = 0; j < d0; j++)
            {
                vects[loc+j] = vects[loc+j]/len;
            }
        }
        else
        {
            lens[i] = 0;
        }
        loc += d0;
    }
    return ;
}


void PQMath::normVects(float *vects, const unsigned int d0,
                       const unsigned int n0)
{
    assert(vects);
    unsigned int i = 0, j = 0, loc = 0;
    float len = 0;

    for(i = 0; i < n0; i++)
    {
        len = 0;
        for(j = 0; j < d0; j++)
        {
            len += vects[loc+j]*vects[loc+j];
        }
        len  = sqrt(len);
        if(len > PQMath::smallVal0)
        {
            for(j = 0; j < d0; j++)
            {
                vects[loc+j] = vects[loc+j]/len;
            }
        }
        loc += d0;
    }
    return ;
}

void PQMath::sqrt_SIFT(float *feats, const unsigned int dim, const unsigned int numb)
{
    assert(feats);
    unsigned int i, j, loc;
    for(i = 0; i < numb; i++)
    {
        loc = i*dim;
        for(j = 0; j < dim; j++)
        {
            feats[loc+j] = sqrt(feats[loc+j]);
        }
    }
    return ;
}

double PQMath::fvec_norm(const float *v, const long n, const double norm)
{
    if(norm == 0)
        return n;

    long i;
    double s = 0;

    if(norm == 1)
    {
        for(i = 0 ; i < n ; i++)
            s += fabs(v[i]);
        return s;
    }

    if(norm == 2)
    {
        for(i = 0 ; i < n ; i++)
        {
            s += v[i]*v[i];
        }
        return sqrt(s);
    }

    if(norm == -1)
    {
        for(i = 0 ; i < n ; i++)
            if(fabs(v[i]) > s)
                s = fabs(v[i]);
        return s;
    }

    for (i = 0 ; i < n ; i++)
    {
        s += pow (v[i], norm);
    }

    return pow (s, 1 / norm);
}


double PQMath::dvec_norm(const double *v, const long n, const double norm)
{
    if(norm == 0)
        return n;

    long i;
    double s = 0;

    if(norm == 1)
    {
        for(i = 0 ; i < n ; i++)
            s += fabs(v[i]);
        return s;
    }

    if(norm == 2)
    {
        for(i = 0 ; i < n ; i++)
        {
            s += v[i]*v[i];
        }
        return sqrt(s);
    }

    if(norm == -1)
    {
        for(i = 0 ; i < n ; i++)
            if(fabs(v[i]) > s)
                s = fabs(v[i]);
        return s;
    }

    for (i = 0 ; i < n ; i++)
    {
        s += pow (v[i], norm);
    }

    return pow (s, 1 / norm);
}

double PQMath::innerProduct(const double *v1, const double *v2, const long n)
{
    double s = 0;
    for(int i = 0 ; i < n ; i++)
    {
        s += v1[i]*v2[i];
    }
    return s;
}

double PQMath::innerProduct(const SparseMatrix &v1, unsigned int s1, const SparseMatrix &v2, unsigned int s2)
{
    unsigned int i = 0;
    double w = 0.0f;

    /**/
    unordered_map<unsigned int, unsigned int> s;
    unordered_map<unsigned int, unsigned int>::iterator it;

    if((v1.col[s1+1] - v2.col[s2+1] - v1.col[s1] + v2.col[s2]) > 0)
    {
        for(i = v1.col[s1]; i < v1.col[s1+1]; i++)
        {
            s.insert(make_pair(v1.index[i],i));
        }

        for(i = v2.col[s2]; i < v2.col[s2+1]; i++)
        {
            it = s.find(v2.index[i]);
            if(it != s.end())
            {
                w += v1.data[it->second] * v2.data[i];
            }
        }

    }
    else
    {
        for(i = v2.col[s2]; i < v2.col[s2+1]; i++)
        {
            s.insert(make_pair(v2.index[i],i));
        }

        for(i = v1.col[s1]; i < v1.col[s1+1]; i++)
        {
            it = s.find(v1.index[i]);
            if(it != s.end())
            {
                w  += v2.data[it->second] * v1.data[i];
            }
        }
    }

    return w;
}

bool PQMath::isZeros(const float *feat, const  unsigned int dim)
{
    unsigned int j;
    float len = 0;
    for(j = 0; j < dim; j++)
    {
        len += feat[j]*feat[j];
    }

    len = sqrt(len);
    if(len <= 0.0000)
    {
        return true;
    }
    else
    {
        cout<<len<<endl;
        return false;
    }
}

int PQMath::sign(const float val)
{
    if(val >=0 )
    {
        return 1;
    }
    else
    {
        return -1;
    }
}


bool PQMath::fvec_scale(float *v, const long n, const double sc)
{
    assert(sc!=0.0f);

    for(long i = 0; i < n; i++)
    {
        v[i] = v[i]*sc;
    }
    return true;
}

bool PQMath::dvec_scale(double *v, const long n, const double sc)
{
    assert(sc!=0.0f);

    for(long i = 0; i < n; i++)
    {
        v[i] = v[i]*sc;
    }
    return true;
}
