#ifndef IOAGENT_H
#define IOAGENT_H

#include "sparsematrix.h"

#include <string>
#include <vector>
#include <map>
#include <set>

using std::set;
using std::map;
using std::vector;
using std::string;
class IOAgent
{
public:
    IOAgent() {}
    virtual ~IOAgent() {}
    static float *load_fvecs(const char *fvecfn, unsigned int &d, unsigned int &r);
    static float *load_fvecs(vector<string> &filenames, unsigned int &d, unsigned int &r);
    static float *loadMatrix(const char *fn,     unsigned int &row, unsigned int &col);
    static float *loadMatrix(vector<string> &filenames,     unsigned int &row, unsigned int &col);
    static float *loadDat(const char *fn,     unsigned int &row, unsigned int &col);
    static float *loadDat(vector<string> &filenames,     unsigned int &row, unsigned int &col);
    static unsigned int *loadPQ(const char *srcfn, unsigned int &row, unsigned int&col);
    static float* loadPQInfo(const char *srcfn, unsigned int &pqn, unsigned int &pql, unsigned &pqm);

    static float *loadItms(const char *fn,  const char *idxKey,  unsigned int &row, unsigned int &col);
    static float *loadItms(const char *fn,  const char *idxKey, const unsigned int line,  unsigned int &row, unsigned int &col);
    static SparseMatrix loadBOVW(const char *fn, unsigned int &row, unsigned int &col);
    static SparseMatrix loadSparse(const char *fn, unsigned int &row, unsigned int &col);
    static float *loadCSV(const char *srcFn, const char delm, unsigned int &row, unsigned int &col);

    static unsigned int getMatrixInfo(const char *matfn, unsigned int &col, unsigned int &row);

    static void saveDat(const char *fn, const unsigned int &row, const unsigned int &col,float *Dat);
    static void saveSparse(const char* fn, const unsigned int& row, const unsigned int &col, const SparseMatrix &sdata);

    bool static readFileName(const char* srcdir, vector<string >& filename);
    static void saveEx(vector<string> &filenames, float *e, float *p, const char *result);
    static void addRslt(const char *filename, float *e, float *p, const char *result);
    static void saveAsFullMat(const float *dat, unsigned int &row, unsigned int &col, const char *file);

/// colnum1 data11 data12 .....datacolnum1
/// colnum2 data21 data22....
    static void loadClust(const char *fn, map<unsigned int, set<int> *> &cluster);
///label to clust
///label1
///label2
///.....
///labeln
    static void loadmClust(const char *srcFn, map<unsigned int, set<int> *> &clusts);
    static void saveClust(const char* fn, map<unsigned int, set<int> *>&cluster);

    static void test();
};

#endif // IOAGENT_H
