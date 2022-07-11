#ifndef DATACOVERT_H
#define DATACOVERT_H

#include "sparsematrix.h"
#include <vector>
#include <map>
#include <set>

using namespace std;

/***
 * Author: Wan-Lei Zhao
 * 
 * since 2014 Sep.
 * 
 * */

class DataConvert
{
public :

    static void getIdf(float *mat, unsigned int row, unsigned col, float *idf);
    static void tfidfMat(const float *mat, const float *idf, unsigned int row, unsigned col, float *nmat);
    static void normMat(float *mat, unsigned int row, unsigned col, float *nmat);
    static void yingzhaoForm(const char *dataFn, const char *srcDir, const char *dstDir);
    static void iggg2ii(const char *igfn, const char * ggfn, const char *iifn);

    static void label2set(vector<int> labels, map<int, set<int> *> &sets);
    static void label2set(vector<pair<int, int> >labels, map<int, set<int > *> &sets);
    static void label2label(const char  *fn, const char* igfn, const char *dstfn);
    static void rawcsv2mat(const char *srcFn, const char *dstFn);
    static void xcsv2mat(const char *srcFn, const char *dstFn);
    static void similarityBoost(SparseMatrix &sdata, unsigned int row, unsigned int ndim, float boostRate);

    static void test();

};



#endif
