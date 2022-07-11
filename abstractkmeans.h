#ifndef ABSTRACTKMEANS_H
#define ABSTRACTKMEANS_H

#include "sparsematrix.h"
#include "clsinfo.h"
#include "nnitem.h"

#include <vector>
#include <map>

struct CLSIndex
{
public:
    CLSIndex()
    {
        index = 0;
        val   = 0.0;
    }
    CLSIndex(const unsigned int idx, const float val0)
    {
        index = idx;
        val   = val0;
    }
    unsigned int index;  /// number of members
    float val;           /// sum of intra-sim

    static bool LtComp (const CLSIndex *a,const CLSIndex *b)
    {
        return (a->val < b->val);
    }
    static bool LgComp (const CLSIndex *a,const CLSIndex *b)
    {
        return (a->val < b->val);
    }
};


enum KMN_OPT {_tkmn_ = 0, _kppkmn_ = 1, _xbkmn_ = 2, _xtkmn_ = 3, _mnkmn_ = 4, _olkmn_ = 5, _rbkmn_, _rrkmn_, _lvq_, _imk_};
enum OPTZ {_I1_ = 0, _I2_ = 1, _I3_, _I4_, _E1_, _E2_,_T1_};
enum SEED {_rnd_ = 0, _kpp_ = 1, _non_, _luc_};

class AbstractKMeans
{
protected:
    AbstractKMeans();

protected:
    static const unsigned int paraBound;
    unsigned int count, ndim, clnumb;
    unsigned int nITER;
    vector<NNItem*>  sorted_stack;
    static const float smallVal0;
    double  *lens, *Es, *arrayD, allEg, *record;
    bool    _INIT_,  _REFER_, _REFINE_;
    char    srcmatfn[1024];
    char    dataFn[32];
    char    mthStr[8];
    unsigned long mvs[2];
    double  dstort;
    SEED    seed;
    float   *data;
    double *kmLogs[128];
    unsigned nLogs;


    SparseMatrix sdata;

    CLSInfo *infoMap;
    CLData *dataMap;
    int     *labels;
    int     *Ns;
    KMN_OPT kmMtd;
    bool dataType;
    bool SHUFFLE;
    OPTZ  myoptz;
public:
    unsigned int realSz;
    unsigned int buildcluster(const char *srcfn, const char *dstfn,  const char *_seed_, const char *lg_first,
                              const char *crtrn, const int num, bool _refine_);
    unsigned int buildcluster(float *mat, const int row, const int dim, const char *dstfn,
                              const char *_seed_, const char  *lg_first, const char *crtrn, const int num, bool _refine_);

    bool  initMemry(const unsigned int dNum, const unsigned int clustNum);
    bool  relsMemry();
    virtual bool  init(const char *srcfn) = 0;
    virtual bool  init(float *data, const int row, const int dim) = 0;
    virtual bool  config(const char *_seed_, const char *lg_first, const char *crtrn, int verbose) = 0;
    virtual int   clust(const unsigned int clust_num, const char *dstfn, const int verbose) = 0;
    unsigned int  nvclust(const unsigned int clnumb, const char *dstfn, const int verbose);
    double  I2FastM(const double *D1, const float *v1, const unsigned int dim, const unsigned int n0);
    double  I2FastP(const double *D1, const float *v1, const unsigned int dim, const unsigned int n0);

    double  getI1(const double *Ds, const unsigned int k, const unsigned int dim, const int *tmpns);
    double  getI2(const double *Ds, const unsigned int k, const unsigned int dim, const int *tmpns);
    double  getI4(const double *Ds, const unsigned int k, const unsigned int dim, const int *tmpns);
    double  getE1(const double *Ds, const unsigned int k, const unsigned int dim, const int *tmpns);
    double  getE2(const double *Ds, const unsigned int k, const unsigned int dim, const int *tmpns);

    double  refine(const unsigned int clustNum, const unsigned int NRef);
    double  refine(int label1, int label2, int label3, const unsigned int NRef);
    double  r2refine(const unsigned int clustNum, const unsigned int NRef);
    bool    BiOptz(int clbl1, int clbl2, double *D1, double *D2, int &n1, int &n2);
    bool    SBiOptz(int clbl1, int clbl2, double *D1, double *D2, int &n1, int &n2);
    double  refineSparse(const unsigned int clustNum, const unsigned int NRef);
    void    setIter(unsigned int iter0)
    {
        nITER = iter0;
    }

    double  calAVGDist(const double *Ds, const unsigned int clustNum, CLSInfo *infos);
    double  calAvgDistort(const unsigned int nclust, int lp);

    float  *getCluster(const unsigned int clabel0, unsigned int &row, unsigned int &dim);
    virtual void saveCenters(const char *dstfn, bool append) = 0;
    void    saveCenters(const char *dstfn, bool append, bool norm);
    virtual int  fetchCenters(float *centers) = 0;
    virtual bool refresh() = 0;

    unsigned int getndim()
    {
        return this->ndim;
    }
    int    idx2clabel(const int i);
    bool   setLogOn(const unsigned logSize);
    void   saveLogs(const char *dstFn, const unsigned nIt);

    void   printClusters(const char *dstdir);
    void   printCluster(const char  *dstfn);
    bool   save_clust(const char    *dstfn);

    static double normVects(float *vects, double *lens, const unsigned int d0, const unsigned int n0);
    static void normVects(const SparseMatrix &svects, double *lens, const unsigned int d0, const unsigned int n0);
    static void normVects(float *vects, const unsigned int d0, const unsigned int n0, double *lens);
    static void normVects(SparseMatrix &vects, const unsigned int d0, const unsigned int n0, double *lens);

    void resetlabels(unsigned int &clustnum);

    virtual ~AbstractKMeans();

};



#endif
