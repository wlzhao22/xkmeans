#include "kmeans.h"

#include "print2scrn.h"
#include "ioagent.h"
#include "cleaner.h"
#include "vstring.h"
#include "nnitem.h"
#include "pqmath.h"
#include "timer.h"

#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <cstdlib>
#include <thread>
#include <vector>
#include <cmath>
#include <ctime>

using namespace std;

const unsigned int KMeans::M_ITER = 40;
const unsigned int KMeans::NTRAILS= 128;
const unsigned int KMeans::NRfn   = 10;
const float KMeans::EPS  = 0.5f;
const float KMeans::Err0 = 0;
const int interval = 1;
KMeans::KMeans()
{
    arrayD = tmparrayD = bstArrayD = Cs = NULL;
    bstLabels = NULL;
    _INIT_ = false;
    kmMtd  = _tkmn_;
    strcpy(mthStr, "_tkm_");
    cout<<"Method ........................... Traditional K-means\n";
}

bool KMeans::init(const char *srcfn)
{
    cout<<"Loading matrix ................... ";
    assert(srcfn);
    strcpy(srcmatfn, srcfn);
    refresh();

    if(VString::endWith(srcfn, ".txt"))
    {
        this->data = IOAgent::loadMatrix(srcfn,  this->count, this->ndim);
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
    else if(VString::endWith(srcfn, ".csv"))
    {
        this->data = IOAgent::loadCSV(srcfn, ',', this->count, this->ndim);
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

    this->bstLabels = new int[this->count];
    memset(this->bstLabels, 0 , sizeof(int)*this->count);

    this->_INIT_  = true;
    this->_REFER_ = false;

    return true;
}

bool KMeans::init(float *mat, const int row, const int dim)
{
    this->refresh();

    this->data  = mat;
    this->count = row;
    this->ndim  = dim;

    this->bstLabels = new int[this->count];
    memset(this->bstLabels, 0 , sizeof(int)*this->count);

    this->_INIT_  = true;
    this->_REFER_ = true;
    return true;
}

bool KMeans::refresh()
{
    this->count  = 1;
    this->ndim   = 0;
    this->_INIT_ = false;

    if(this->_REFER_)
    {
        if(this->data != NULL)
        {
            this->data = NULL;
        }
    }
    else
    {
        if(this->data != NULL)
        {
            delete [] this->data;
            this->data = NULL;
        }
    }

    if(this->bstLabels != NULL)
    {
        delete [] this->bstLabels;
        this->bstLabels = NULL;
    }

    if(this->arrayD != NULL)
    {
        delete [] this->arrayD;
        this->arrayD = NULL;
    }

    if(this->arrayDs != NULL)
    {
        delete [] this->arrayDs;
        this->arrayDs = NULL;
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

    if(this->Ns != NULL)
    {
        delete [] this->Ns;
        this->Ns = NULL;
    }

    return true;
}

bool KMeans::config(const char *_seed_, const char *crtrn, const char *lg_first, int verbose)
{
    if(verbose)
        cout<<"Distance function ................ l2\n";

    if(verbose)
        cout<<"Seeds ............................ ";
    if(!strcmp(_seed_, "rnd"))
    {
        if(verbose)
            cout<<"rand\n";
        sedFunc = &KMeans::rndSeeds;
        seed = _rnd_;
    }
    else if(!strcmp(_seed_, "kpp"))
    {
        if(verbose)
            cout<<"kpp\n";
        sedFunc = &KMeans::kppSeeds;
        seed = _kpp_;
    }else if(!strcmp(_seed_, "luc"))
    {
        if(verbose)
            cout<<"kpp\n";
        sedFunc = &KMeans::lucSeeds;
        seed = _luc_;
    }
    else
    {
        if(verbose)
            cout<<"kpp\n";
        sedFunc = &KMeans::kppSeeds;
        seed = _kpp_;
    }

    return true;
}

bool KMeans::lucSeeds(const int k, int rseeds[], const int bound)
{
    unsigned int i = 0, j = 0, row = 0, col = 0;
    const char *srcFn0 = "/home/wlzhao/datasets/clust/sift1m/seed";
    char srcFn[512];

    sprintf(srcFn, "%s%d.txt", srcFn0, bound);
    float *mat = IOAgent::loadMatrix(srcFn, row, col);
    cout<<"seeding is done\t"<<row<<"\t"<<col<<endl;
    for(i = 0; i < row; i++)
    {
       for(j = 0; j < col; j++)
       {
          this->Cs[i*col+j] = mat[i*col+j];
       }
    }
    delete [] mat;
    mat = NULL;
    return true;
}

int KMeans::getL2norms(const unsigned int n, const unsigned int d0)
{
   unsigned int i = 0, j = 0,loc = 0;
   assert(this->data != NULL);
   assert(this->lens != NULL);

   for(i = 0; i < n; i++)
   {
        this->lens[i] = 0;
        loc = i*d0;
        for(j = 0; j < d0; j++)
        {
            this->lens[i] += this->data[loc+j]*this->data[loc+j];
        }
   }
   return 0;
}

double KMeans::pairwDst(const int nclust)
{
    unsigned int i = 0, label = 0, cloc = 0, c0 = 0, nl = 0;
    double l2nsum = 0, pd = 0, *D;
    for(label = 0; label < nclust; label++)
    {
        cloc = label * this->ndim;
        D    = this->arrayDs+cloc;
        for(i = 0; i < this->ndim; i++)
        {
            pd += D[i]*D[i];
        }
    }
    for(i = 0; i < this->count; i++)
    {
        c0 = this->labels[i];
        nl = this->Ns[c0];
        l2nsum += nl * this->lens[i];
    }
    return (l2nsum-pd)/this->count;
}

int KMeans::clust(const unsigned int nclust, const char *dstfn, const int verbose)
{
    double dst = 0, minDst = 0.0f, sumDst0 = 0, err = 0.0f;
    double sumDst = 0, bstSumDst = RAND_MAX + 0.0f, optEg = 0, prvEg = 0, avgDst = 0;
    unsigned int loc = 0, d = 0, di = 0, niter = 0, c = 0;
    unsigned int i = 0, j = 0, clabel = 0;
    unsigned int k = 0, ri = 0;
    unsigned long cmps = 0;

    if(!this->_INIT_)
    {
        cout<<"Error ........................... ahead of init!\n";
        return 0;
    }

    int cpysize = this->ndim*nclust*sizeof(double);
    Cs        = new double[this->ndim*nclust];
    bstArrayD = new double[this->ndim*nclust];
    memset(bstArrayD, 0, cpysize);
    memset(Cs,        0, cpysize);

    if(this->record == NULL)
    {
        this->record = new double[M_ITER*6];
    }
    memset(this->record, 0, sizeof(double)*M_ITER*6);
    cout<<"Initialize record ................ ";
    for(unsigned i = 0; i < M_ITER; i++)
    {
        this->record[6*i]   = numeric_limits<double>::max();
        this->record[6*i+1] = -1;
        this->record[6*i+3] = numeric_limits<double>::max();
        this->record[6*i+4] = -1;
    }
    cout<<"done"<<endl;

    int    *rseeds = new int[nclust];
    double *infos  = new double[nclust];
    this->Ns = new int[nclust];
    double *tmpDat = new double[ndim];
    cout<<"cluster num. ..................... "<<nclust<<endl;
    cout<<"Start clustering ................. \n";
    clock_t t1 = clock();

    for(ri = 0; ri < NTRAILS; ri++)
    {
        (this->*sedFunc)(nclust, rseeds, count);
        this->getL2norms(this->count, this->ndim);
        for(j = 0; j < nclust; j++)
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
            memset(arrayD,   0, cpysize);
            memset(infos,    0, nclust*sizeof(double));
            memset(this->Ns, 0, nclust*sizeof(int));
            ///compute data to arrayD
            for(i = 0; i < count; i++)
            {
                minDst = RAND_MAX+0.0f;
                loc = this->ndim*i;
                for(d = 0; d < ndim; d++)
                {
                    tmpDat[d] = this->data[loc+d];
                }
                for(k = 0; k < nclust; k++)
                {
                    if(this->Ns[k] == -1)
                        continue;

                    dst  = PQMath::l2d(Cs, k, tmpDat, 0, this->ndim);
                    ///dst = PQMath::l1d(Cs, k, tmpDat, 0, this->ndim);
                    ///cout<<k<<"\t"<<dst<<endl;
                    if(dst < minDst)
                    {
                        clabel = k;
                        minDst = dst;
                    }
                    cmps++;
                }

                this->labels[i]   = clabel;
                this->Ns[clabel] += 1;
                infos[clabel]    += minDst;
                di  = clabel*this->ndim;

                for(d = 0; d < ndim; d++)
                {
                    this->arrayD[di+d] += tmpDat[d];
                }
            }///for(i)

            ///update center: arrayD
            sumDst = 0;
            for(j = 0; j < nclust; j++)
            {
                if(this->Ns[j] == 0)
                {
                    di = j*ndim;
                    memset(this->arrayD+di, 0, sizeof(float)*ndim);
                    memset(this->Cs+di,     0, sizeof(float)*ndim);
                    infos[j]  = 0;
                    this->Ns[j] = -1;
                }
                else
                {
                    di = j*ndim;
                    for(d = 0; d < ndim; d++)
                    {
                        this->Cs[di+d] = this->arrayD[di+d]/this->Ns[j];
                    }

                    sumDst += infos[j];
                }
            }///for(j)
            niter++;

            if(this->arrayDs == NULL)
            {
               this->arrayDs = new double[nclust*this->ndim];
            }
            memset(this->arrayDs, 0, sizeof(double)*nclust*this->ndim);
            for(i = 0; i < this->count; i++)
            {
                 c = this->labels[i];
                 loc = c*this->ndim;
                 for(j = 0; j < this->ndim; j++)
                 {
                     this->arrayDs[loc+j] += this->data[i*this->ndim+j];
                 }
            }

            avgDst = pairwDst(nclust);
            sumDst = this->calAvgDistort(nclust, 2);
            sumDst0 = sumDst;
            double t2 = (clock() - t1 + 0.0)/CLOCKS_PER_SEC;
            cout<<"iter "<< niter <<" avgDist :\t"<<sumDst<<"\t"<<avgDst<<"\t"<<t2<<endl;

            int iter = niter - 1;

            if(this->record[6*iter] > avgDst)
                this->record[6*iter] = avgDst;
            if(this->record[6*iter+1] < avgDst)
                this->record[6*iter+1] = avgDst;
            this->record[6*iter+2] += avgDst;

            if(this->record[6*iter+3] > sumDst)
                this->record[6*iter+3] = sumDst;
            if(this->record[6*iter+4] < sumDst)
                this->record[6*iter+4] = sumDst;
            this->record[6*iter+5] += sumDst;

        }while(niter < KMeans::M_ITER);

        if(sumDst0 < bstSumDst)
        {
            memcpy(bstArrayD, arrayD, this->ndim*nclust*sizeof(double));
            memcpy(bstLabels, labels, count*sizeof(int));
            bstSumDst = sumDst0;
        }
    }///for(ri)

    for(unsigned i = 0; i < M_ITER; i++)
    {
        cout<<"avgI2: "<<i<<"\t"<<this->record[6*i]<<"\t"<<this->record[6*i+1]<<"\t"<<this->record[6*i+2]/NTRAILS<<endl;
        cout<<"avgI1: "<<i<<"\t"<<this->record[6*i+3]<<"\t"<<this->record[6*i+4]<<"\t"<<this->record[6*i+5]/NTRAILS<<endl;
    }

    memcpy(labels, bstLabels, count*sizeof(int));
    cout<<"Comparisons: "<<cmps<<endl;

    memcpy(arrayD, bstArrayD, this->ndim*nclust*sizeof(double));

    double sumEg = this->calAvgDistort(nclust, 2);
    cout<<"best: "<<sumEg<<endl;

    if(false)
    {
        refine(nclust, 1);
        double sumEg = this->calAvgDistort(nclust, 2);
    }

    for(j = 0; j < nclust; j++)
    {
        NNItem *crnt_item = new NNItem(j, this->infoMap[j].E, this->infoMap[j].n);
        sorted_stack.push_back(crnt_item);
    }

    if(strlen(dstfn) > 0)
    {
        char infofn[512];
        VString::parsePath(infofn, dstfn);

        strcat(infofn, "_km_info.txt");
        Print2Scrn::printvect(sorted_stack, infofn);
        KMeans::save_clust(dstfn);
    }

    Cleaner::clearVector(sorted_stack);

    delete [] bstLabels;
    delete [] tmparrayD;
    delete [] bstArrayD;
    delete [] infos;
    delete [] rseeds;

    bstLabels = NULL;
    tmparrayD = NULL;
    bstArrayD = NULL;
    rseeds = NULL;
    infos  = NULL;

    return nclust;
}

bool KMeans::rndSeeds(const int k, int rseeds[], const int bound)
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

/** random k seeds with the way of k-means++ **/
bool KMeans::kppSeeds(const int k, int rseeds[], const int bound)
{
    int i = 0, j = 0;
    double sum = 0.0f, rd = 0;
    int sel;


    for(i = 0; i < k; i++)
    {
        rseeds[i] = i;
    }

    double *disbest = new double[bound];
    for(i = 0; i < bound; i++)
    {
        disbest[i] = RAND_MAX + 0.0f;
    }
    double * distmp    = new double[bound];
    unsigned int sed0 = (unsigned int)time(NULL);
    float tmp = 0;

    rseeds[0] = rand_r(&sed0) % bound;
    for (i = 1 ; i < k; i++)
    {
        sel = rseeds[i - 1];
        for (j = 0 ; j < bound; j++)
        {
            tmp = PQMath::l2f(this->data, j, this->data, sel, ndim);
            if(tmp < disbest[j])
            {
                disbest[j] = tmp;
            }
        }

        /** convert the best distances to probabilities **/
        memcpy (distmp, disbest, bound * sizeof (double));

        sum = PQMath::dvec_norm(distmp, bound, 1);
        PQMath::dvec_scale(distmp, bound, 1.0/sum);
        rd  = rand_r(&sed0)/(RAND_MAX+1.0f);

        for (j = 0 ; j < bound-1; j++)
        {
            rd -= distmp[j];
            if (rd < 0)
                break;
        }
        rseeds[i] = j;
    }

    delete [] disbest;
    delete [] distmp;
    disbest = distmp = NULL;
    return true;
}

void KMeans::saveCenters(const char *dstfn, bool append)
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

int KMeans::fetchCenters(float *centers)
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

KMeans::~KMeans()
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

void KMeans::test()
{
   const char *srcMatFn1 = "/home/wlzhao/datasets/bignn/sift1m/sift_base.txt";
   ///const char *srcMatFn2 = "/home/wlzhao/datasets/bignn/sift1m/sift_learn.txt";
   const char *srcMatFn2 = "/home/rqchen/dataset/sift_learn.txt";
   const char *srcMatFn3 = "/home/wlzhao/datasets/bignn/glove/glove1m_norm_base.txt";
   const char *srcfn = "sift_learn.txt";
   ///const char *dstfn = "/home/wlzhao/datasets/hi.txt";
   const char *dstfn = "/home/rqchen/dataset/rslt/hi.txt";
   //int i = 1024;
  // do{
       KMeans *mykm = new KMeans();
       ///mykm->buildcluster(srcMatFn2, dstfn, "rnd", "large", "i2", 1024, false);
       mykm->buildcluster(srcMatFn2, dstfn, "kpp", "large", "i2", 1024, false);
       delete mykm;
   //    i = i*2;
  // }while(i <= 8192);
}
