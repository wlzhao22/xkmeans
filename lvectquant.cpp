#include "lvectquant.h"
#include "vstring.h"
#include "ioagent.h"
#include "cleaner.h"
#include "pqmath.h"


#include <iostream>
#include <fstream>
#include <cstring>
#include <cassert>
#include <cmath>

using namespace std;

const unsigned int LVectQuant::NTRAILS = 1;
const unsigned int LVectQuant::M_ITER  = 30;
const float LVectQuant::alpha0   = 0.01;

LVectQuant::LVectQuant()
{
    Ds = Cs = tmpCs   = NULL;
    kmMtd   = _lvq_;
    strcpy(mthStr, "_lvq_");
    this->_REFER_ = false;
    cout<<"Method ........................... Learning Vector Quantization\n";
}

bool LVectQuant::init(const char *srcfn)
{
    refresh();
    cout<<"Loading matrix ................... ";
    assert(srcfn);
    strcpy(srcmatfn, srcfn);
    if(VString::endWith(srcfn, ".txt"))
    {
        this->data = IOAgent::loadMatrix(srcfn, this->count, this->ndim);
    }
    else if(VString::endWith(srcfn, ".fvecs"))
    {
        this->data = IOAgent::load_fvecs(srcfn, this->ndim, this->count);
    }
    else if(VString::endWith(srcfn, ".mat"))
    {
        this->data = IOAgent::loadDat(srcfn, this->count, this->ndim);
    }
    else if(VString::endWith(srcfn, ".itm"))
    {
        this->data = IOAgent::loadItms(srcfn, "fsvtab", this->count, this->ndim);
    }
    else
    {
        cout<<"Unrecognizable input file format!!!\n";
        this->data = NULL;
        exit(0);
    }

    cout<<this->count<<"x"<<this->ndim<<endl;

    if(this->data == NULL)
    {
        cout<<"Exceptions ....................... Loading matrix failed!\n";
        exit(0);
    }

    this->_REFER_= false;
    return true;
}

bool LVectQuant::init(float *mat, const int row, const int dim)
{
    assert(mat);
    refresh();

    this->data   = mat;
    this->count  = row;
    this->ndim   = dim;

    this->_INIT_  = true;
    this->_REFER_ = true;
    return true;
}

bool LVectQuant::refresh()
{
    this->count  = 1;
    this->ndim   = 0;
    this->_INIT_ = false;

    if(!this->_REFER_)
    {
        if(this->data != NULL)
        {
            delete [] this->data;
            this->data = NULL;
        }
    }
    else
    {
        if(this->data != NULL)
        {
            this->data = NULL;
        }
    }

    return true;
}

bool LVectQuant::config(const char *_seed_, const char *lg_first, const char *crtrn, int verbose)
{
    if(verbose)
        cout<<"Seeds ............................ ";

    if(verbose)
        cout<<"rand\n";
    sedFunc = &LVectQuant::rndSeeds;
    return true;
}

int LVectQuant::clust(const unsigned int clust_num, const char *dstfn, const int verbose)
{
    double dst = 0, minDst = 0.0f, sumDst0 = 0;
    double sumDst = 0, bstSumDst = RAND_MAX + 0.0f, optEg = 0;
    unsigned int loc = 0, d = 0, di = 0, niter = 0;
    unsigned int i = 0, j = 0, clabel = 0;
    unsigned int k = 0, ri = 0;
    unsigned long cmps = 0;

    if(!this->_INIT_)
    {
        cout<<"Error ........................... ahead of init!\n";
        return 0;
    }


    int cpysize = this->ndim*clust_num*sizeof(double);
    Cs        = new double[this->ndim*clust_num];
    tmparrayD = new double[this->ndim*clust_num];
    bstArrayD = new double[this->ndim*clust_num];
    memset(tmparrayD, 0, cpysize);
    memset(bstArrayD, 0, cpysize);
    memset(Cs,        0, cpysize);

    int    *rseeds = new int[clust_num];
    double *infos  = new double[clust_num];
    int    *counts = new int[clust_num];
    double *tmpDat = new double[ndim];
    bstLabels = new int[count];
    float alpha =  alpha0;
    //ofstream *foutDist = new ofstream("/home/wlzhao/datasets/clustrslt/lvq_10k_001.txt", ios::out);

    for(ri = 0; ri < LVectQuant::NTRAILS; ri++)
    {
        cmps = 0;
        (this->*sedFunc)(clust_num, rseeds, count);
        for(j = 0; j < clust_num; j++)
        {
            cmps++;
            loc = this->ndim*rseeds[j];
            di  = this->ndim*j;
            for(d = 0; d < ndim; d++)
            {
                Cs[di+d] = this->data[loc+d];
            }
        }

        niter  = 0;
        do
        {
            memset(tmparrayD, 0, cpysize);
            memset(infos,  0, clust_num*sizeof(double));
            memset(counts, 0, clust_num*sizeof(int));
            ///compute data to arrayD
            for(i = 0; i < count; i++)
            {
                minDst = RAND_MAX+0.0f;
                loc = this->ndim*i;
                for(d = 0; d < ndim; d++)
                {
                    tmpDat[d] = this->data[loc+d];
                }
                for(k = 0; k < clust_num; k++)
                {
                    if(counts[k] == -1)
                        continue;

                    dst  = PQMath::l2d(Cs, k, tmpDat, 0, this->ndim);
                    if(dst < minDst)
                    {
                        clabel = k;
                        minDst = dst;
                    }
                    cmps++;
                }

                this->labels[i] = clabel;
                counts[clabel] += 1;
                infos[clabel]  += minDst;
                di  = clabel*this->ndim;

                /**incremental update on the clustering center**/
                for(d = 0; d < ndim; d++)
                {
                    this->Cs[di+d] += alpha*(tmpDat[d] - this->Cs[di+d]);
                }
                /**/

                for(d = 0; d < ndim; d++)
                {
                    this->tmparrayD[di+d] += tmpDat[d];
                }
            }///for(i)
            alpha = alpha - 0.00004;

            ///adjust the centroid for next round by taking average of a cluster
            sumDst = 0;
            for(j = 0; j < clust_num; j++)
            {
                if(counts[j] == 0)
                {
                    di = j*ndim;
                    memset(this->arrayD+di, 0, sizeof(float)*ndim);
                    memset(this->Cs+di,     0, sizeof(float)*ndim);
                    infos[j]  = 0;
                    counts[j] = -1;
                }
                else
                {
                    di = j*ndim;
                    for(d = 0; d < ndim; d++)
                    {
                        this->Cs[di+d] = this->tmparrayD[di+d]/counts[j];
                    }

                    sumDst += infos[j];
                }
            }///for(j)
            niter++;

            optEg = getI2(this->tmparrayD, clust_num, ndim, counts);
            sumDst0 = sumDst;
            //(*foutDist)<<niter<<"\t"<<(allEg - optEg)/this->count<<endl;
        }while(niter < LVectQuant::M_ITER);
        if(sumDst0 < bstSumDst)
        {
            memcpy(arrayD, tmparrayD, this->ndim*clust_num*sizeof(double));
            memcpy(bstLabels, labels, count*sizeof(int));
            bstSumDst = sumDst0;
            for(j = 0; j < clust_num; j++)
            {
                this->Ns[j] = (counts[j]==-1)?0:counts[j];
            }
        }
    }///for(ri)
    //foutDist->close();

    memcpy(labels, bstLabels, count*sizeof(int));
    cout<<"Comparisons: "<<cmps<<endl;

    double sumEg = calAVGDist(this->arrayD, clust_num, this->infoMap);
    cout<<"before refine: "<<sumEg<<endl;

    if(false)
    {
        refine(clust_num, 1);
        double sumEg = calAVGDist(this->arrayD, clust_num, this->infoMap);
        cout<<"after refine: "<<sumEg<<endl;
    }

    for(j = 0; j < clust_num; j++)
    {
        NNItem *crnt_item = new NNItem(j, this->infoMap[j].E, this->infoMap[j].n);
        sorted_stack.push_back(crnt_item);
    }

    if(strlen(dstfn) > 0)
    {
        char infofn[512];
        VString::parsePath(infofn, dstfn);

        strcat(infofn, "_km_info.txt");
        LVectQuant::save_clust(dstfn);
    }

    Cleaner::clearVector(sorted_stack);

    delete [] bstLabels;
    delete [] tmparrayD;
    delete [] bstArrayD;
    delete [] tmpDat;
    delete [] counts;
    delete [] infos;
    delete [] rseeds;

    bstLabels = NULL;
    tmparrayD = NULL;
    bstArrayD = NULL;
    rseeds = NULL;
    infos  = NULL;
    counts = NULL;
    tmpDat = NULL;

    return clust_num;
}

bool LVectQuant::rndSeeds(const int k, int rseeds[], const int bound)
{
    unsigned int NITER, it = 0, sed0 = (unsigned int)time(NULL);
    float v  = rand_r(&sed0)/(RAND_MAX+1.0f);
    bool FAILED = false;
    int i = 1, j = 0;
    int r = 0;

    rseeds[0] = (int)floor(v*bound);
    rseeds[0] = (rseeds[0]==bound)?(rseeds[0] - 1):rseeds[0];
    NITER = 20*k;

    while(i < k && it < NITER)
    {
        v      = rand_r(&sed0)/(RAND_MAX+1.0f);
        r      = (int)floor(v*bound);
        r      = (r==bound)?(r-1):r;
        FAILED = false;
        for(j = 0; j < i; j++)
        {
            if(r == rseeds[j])
            {
                FAILED = true;
                break;
            }
        }
        if(!FAILED)
        {
            rseeds[i] = r;
            i++;
        }
        it++;
    }

    if(FAILED && it == NITER)
    {
        return false;
    }
    else
    {
        return true;
    }
}

void LVectQuant::saveCenters(const char *dstfn, bool append)
{
    unsigned int clabel = 0, j, loc, rCNum = 0;
    ofstream *outStrm  = NULL;

    for(clabel = 0; clabel < this->clnumb; clabel++)
    {
        if(this->infoMap[clabel].n > 0)
        {
            rCNum++;
        }
    }

    if(!this->_INIT_||rCNum == 0)
    {
        return ;
    }

    if(append)
    {
        outStrm = new ofstream(dstfn, ios::app);
    }
    else
    {
        outStrm = new ofstream(dstfn, ios::out);
        (*outStrm)<<rCNum<<" "<<this->ndim<<endl;;
    }
    for(clabel = 0; clabel < this->clnumb; clabel++)
    {
        loc  = clabel*this->ndim;

        if(this->infoMap[clabel].n <= 0)
            continue;

        for(j = 0; j < this->ndim; j++)
        {
            (*outStrm)<<this->arrayD[loc+j]/this->infoMap[clabel].n<<" ";
        }
        (*outStrm)<<endl;
    }

    outStrm->close();
    cout<<"done\n";
    return ;
}

int  LVectQuant::fetchCenters(float *centers)
{
    unsigned int clabel = 0, j = 0, loc = 0, idxi = 0, rCNum = 0;
    assert(centers);

    for(clabel = 0; clabel < this->clnumb; clabel++)
    {
        if(this->infoMap[clabel].n > 0)
        {
            rCNum++;
        }
    }

    if(!this->_INIT_||rCNum == 0)
    {
        memset(centers, 0, this->clnumb*this->ndim*sizeof(float));
        return 0;
    }

    for(clabel = 0; clabel < this->clnumb; clabel++)
    {
        loc   = clabel*this->ndim;

        if(this->infoMap[clabel].n <= 0)
            continue;

        for(j = 0; j < this->ndim; j++)
        {
            centers[idxi + j] = this->arrayD[loc+j]/this->infoMap[clabel].n;
        }
        idxi += this->ndim;
    }
    return rCNum;
}

LVectQuant::~LVectQuant()
{
    if(this->bstLabels != NULL)
    {
        delete [] this->bstLabels;
        this->bstLabels = NULL;
    }

    if(this->infoMap != NULL)
    {
        delete [] this->infoMap;
        this->infoMap = NULL;
    }

    if(this->arrayD != NULL)
    {
        delete [] this->arrayD;
        this->arrayD = NULL;
    }

    if(this->tmparrayD != NULL)
    {
        delete [] this->tmparrayD;
        this->tmparrayD = NULL;
    }

    if(this->Cs != NULL)
    {
        delete [] this->Cs;
        this->Cs = NULL;
    }
}

void LVectQuant::test()
{
   const char *srcMatFn1 = "/home/wlzhao/datasets/bignn/mat/sift1m/sift_base.txt";
   const char *srcMatFn2 = "/home/wlzhao/datasets/bignn/mat/sift_learn.txt";
   const char *srcMatFn3 = "itm_vlad_hesaff_flickr10m.itm";
   const char *srcfn = "sift_learn.txt";
   const char *dstfn = "itm_cvlad_dense_tst.result10";
   int i = 1024;
   //do{
       LVectQuant *mykm = new LVectQuant();
       mykm->buildcluster(srcfn, dstfn, "rnd", "large", "i2", 1024, false);
       delete mykm;
      i = i*2;
   //}while(i <= 8192);
}

