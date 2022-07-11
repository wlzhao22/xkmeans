#include "xtksums.h"

#include <iostream>
#include <cassert>
#include <cstring>
#include <fstream>
#include <queue>
#include <cmath>

#include "vstring.h"
#include "ioagent.h"
#include "cleaner.h"
#include "pqmath.h"
#include "timer.h"

//#define _DB_

using namespace std;

const int XTKSums::nTrails = 128;
const int XTKSums::nIter   = 41;    ///set iteration times

XTKSums::XTKSums()
{
    this->_INIT_    = false;
    this->bstLabels = NULL;
    this->l2norms   = NULL;
    this->Ds        = NULL;
    this->verbose   = 1;
}

bool XTKSums::init(const char *srcfn)
{
    cout<<"Method ........................... xtksum"<<endl;
    assert(srcfn);
    strcpy(srcmatfn, srcfn);
    VString::parseFName(dataFn, srcmatfn);
    cout<<"Dataset .......................... "<<dataFn<<endl;
    cout<<"Loading matrix ................... ";
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
        this->data = IOAgent::loadItms(srcfn, "fsvtab", 1000000, this->count, this->ndim);
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
    cout<<"nIter is set to .................. "<<this->nIter<<endl;

    if(this->data == NULL)
    {
        cout<<"Exceptions ....................... Loading matrix failed!\n";
        exit(0);
    }

    if(this->labels == NULL)
    {
        this->labels = new int[this->count];
    }
    memset(this->labels, 0, sizeof(int)*this->count);

    if(this->bstLabels == NULL)
    {
       this->bstLabels = new int[this->count];
    }
    memset(this->bstLabels, 0, sizeof(int)*this->count);

    if(this->l2norms == NULL)
    {
        this->l2norms = new double[this->count];
    }
    memset(this->l2norms, 0, sizeof(double)*this->count);

    this->getL2norms(this->count, this->ndim);
    this->_INIT_  = true;
    this->_REFER_ = false;
    return true;
}

bool XTKSums::init(const char *srcfn, unsigned int fixn)
{
    cout<<"Method ........................... xtksum"<<endl;
    refresh();
    assert(srcfn);
    strcpy(srcmatfn, srcfn);
    VString::parseFName(dataFn, srcmatfn);
    cout<<"Dataset .......................... "<<dataFn<<endl;
    cout<<"Loading matrix ................... ";

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
        this->data = IOAgent::loadItms(srcfn, "fsvtab", 1000000, this->count, this->ndim);
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
    cout<<"nIter is set to .................. "<<this->nIter<<endl;

    if(this->data == NULL)
    {
        cout<<"Exceptions ....................... Loading matrix failed!\n";
        exit(0);
    }
    this->count = fixn;

    memset(this->labels, 0, sizeof(int)*this->count);

    if(this->bstLabels == NULL)
    {
       this->bstLabels = new int[this->count];
    }
    memset(this->bstLabels, 0, sizeof(int)*this->count);

    if(this->l2norms == NULL)
    {
        this->l2norms = new double[this->count];
    }
    memset(this->l2norms, 0, sizeof(double)*this->count);

    this->getL2norms(this->count, this->ndim);

    /**
    if(this->dists[0] == NULL)
    {
        this->dists[0] = new double[this->count];
    }
    if(this->dists[1] == NULL)
    {
        this->dists[1] = new double[this->count];
    }
    if(this->dists[2] == NULL)
    {
       this->dists[2] = new double[this->count];
    }
    if(this->dists[3] == NULL)
    {
        this->dists[3] = new double[this->count];
    }/**/
    /**
    memset(this->dists[0], 0, sizeof(double)*this->count);
    memset(this->dists[1], 0, sizeof(double)*this->count);
    memset(this->dists[2], 0, sizeof(double)*this->count);
    memset(this->dists[3], 0, sizeof(double)*this->count);
    /**/

    this->_INIT_  = true;
    this->_REFER_ = false;
    return true;
}


bool XTKSums::init(float *mat, const int row, const int dim)
{
    this->refresh();

    this->data  = mat;
    this->count = row;
    this->ndim  = dim;

    if(this->labels == NULL)
    {
       this->labels = new int[this->count];
    }
    memset(this->labels, 0, sizeof(int)*this->count);

    if(this->bstLabels == NULL)
    {
       this->bstLabels = new int[this->count];
    }
    memset(this->bstLabels, 0, sizeof(int)*this->count);

    this->l2norms = new double[this->count];
    memset(this->l2norms, 0, sizeof(double)*this->count);
    this->getL2norms(this->count, this->ndim);

    this->_INIT_  = true;
    this->_REFER_ = true;
    return true;
}

bool XTKSums::config(const char *_seed_, const char *crtrn, const char *lg_first, int verbose)
{
    if(verbose)
        cout<<"Distance function ................ l2\n";

    if(verbose)
        cout<<"Seeds ............................ ";

    {
        if(verbose)
            cout<<"non\n";
        seed = _non_;
    }

    return true;
}

bool XTKSums::refresh()
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

    if(this->l2norms != NULL)
    {
        delete [] this->l2norms;
        this->l2norms = NULL;
    }

    if(this->Ds != NULL)
    {
        delete [] this->Ds;
        this->Ds = NULL;
    }

    return true;
}

int XTKSums::initNon(const unsigned int nclust)
{
    unsigned int i = 0, j = 0, c = 0, loc = 0;
    memset(this->arrayD, 0, sizeof(double)*nclust*this->ndim);
    memset(this->labels, 0, sizeof(int)*this->count);
    memset(this->Ns, 0, sizeof(int)*nclust);

    if(this->Ds == NULL)
    {
        this->Ds = new double[nclust];
    }
    memset(this->Ds, 0, sizeof(double)*nclust);

    cout<<"Initialize cluster ............... ";
    for(i = 0; i < this->count; i++)
    {
         c = rand()%nclust;
         this->labels[i] = c;
         this->Ns[c] += 1;
    }

    for(i = 0; i < this->count; i++)
    {
         c   = this->labels[i];
         loc = c*this->ndim;
         for(j = 0; j < this->ndim; j++)
         {
             this->arrayD[loc+j] += this->data[i*this->ndim+j];
         }
         this->Ds[c] += this->l2norms[i];
    }

    cout<<"done\n";
    return 0;
}

int XTKSums::initI2(const unsigned int nclust)
{
    unsigned int i = 0, j = 0, c = 0, loc = 0;
    memset(this->arrayD, 0, sizeof(double)*nclust*this->ndim);
    memset(this->labels, 0, sizeof(int)*this->count);

    if(this->Ds == NULL)
    {
        this->Ds = new double[nclust];
    }

    memset(this->Ds, 0, sizeof(double)*nclust);
    memset(this->Ns, 0, sizeof(int)*nclust);

    cout<<"Initialize cluster ............... ";
    for(i = 0; i < this->count; i++)
    {
         c = rand()%nclust;
         this->labels[i] = c;
         this->Ns[c] += 1;
    }

    for(i = 0; i < this->count; i++)
    {
         c = this->labels[i];
         loc = c*this->ndim;
         for(j = 0; j < this->ndim; j++)
         {
             this->arrayD[loc+j] += this->data[i*this->ndim+j];

         }
         this->Ds[c] += this->l2norms[i];
    }

    cout<<"done\n";
    return 0;
}


int XTKSums::getL2norms(const unsigned int n, const unsigned int d0)
{
   unsigned int i = 0, j = 0,loc = 0;
   assert(this->data    != NULL);
   assert(this->l2norms != NULL);

   for(i = 0; i < n; i++)
   {
        this->l2norms[i] = 0;
        loc = i*d0;
        for(j = 0; j < d0; j++)
        {
            this->l2norms[i] += this->data[loc+j]*this->data[loc+j];
        }
   }
   return 0;
}

int XTKSums::clust(const unsigned int nclust, const char *dstFn, const int verbose)
{
    double avgDst = 0, minDst = RAND_MAX + 0.0;

    cout<<"Cluster numb ..................... "<<nclust<<endl;
    cout<<"Optimize ......................... \n";

    double t2 = 0;
    clock_t t1;
    ofstream outStrm;
    this->bg = clock();

    if(this->record == NULL)
    {
        this->record = new double[nIter*6];
    }
    memset(this->record, 0, sizeof(double)*nIter*6);

    cout<<"Initialize record ................ ";
    for(unsigned i = 0; i < nIter; i++)
    {
        this->record[6*i]   = numeric_limits<double>::max();
        this->record[6*i+1] = -1;
        this->record[6*i+2] = 0;
        this->record[6*i+3] = numeric_limits<double>::max();
        this->record[6*i+4] = -1;
        this->record[6*i+5] = 0;
    }
    cout<<"done"<<endl;

    for(unsigned i = 0; i < nTrails; i++)
    {
       this->initI2(nclust);
       t1 = clock();
       avgDst = optzI1(nIter, nclust);
       t2 = (clock() - this->bg+0.0)/CLOCKS_PER_SEC;
       outStrm.open(dstFn, ios::app);
       outStrm<<nclust <<"\t"<<avgDst<<"\t"<<t2<<endl;
       cout<<i<<"\t" << nclust <<" ........."<<avgDst<<"\t"<<t2<<endl;
       outStrm.close();
       if(avgDst < minDst)
       {
           minDst = avgDst;
           memcpy(this->bstLabels, this->labels, this->count*sizeof(int));
       }
    }

    cout<<"\n\t\tdone\n";

    /**/
    ofstream logStrm;
    char destLog[512];
    sprintf(destLog, "./%s_ksum_Im_%d.txt", dataFn, nclust);
    logStrm.open(destLog, ios::out);
    if(!logStrm.is_open())
    {
        cout<<"Log file cannot open for write!\n";
        exit(0);
    }
    logStrm<<"nIter\tI1\t\t"<<"I2\n";
    for(unsigned i = 0; i < nIter; i++)
    {
        logStrm<<i+1<<"\t"<<this->record[6*i]<<"\t\t"<<this->record[6*i+1]<<"\t\t"<<this->record[6*i+2]<<endl;
    }
    cout<<"-------------------------------\n";
    for(unsigned i = 0; i < nIter; i++)
    {
        logStrm<<i+1<<"\t"<<this->record[6*i+3]<<"\t\t"<<this->record[6*i+4]<<"\t\t"<<this->record[6*i+5]<<endl;
    }
    logStrm.close();
    delete [] this->record;
    this->record = NULL;
    /**/

    memcpy(this->labels, this->bstLabels, this->count*sizeof(int));

    cout<<"Averge distortion ................ "<<minDst<<endl;

    return 0;
}

void XTKSums::saveCenters(const char *dstfn, bool append)
{
    unsigned int clabel = 0, j = 0, rCNum = 0, cloc = 0;
    ofstream *outStrm  = NULL;

    for(clabel = 0; clabel < this->clnumb; clabel++)
    {
        if(this->Ns[clabel] > 0)
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
        if(this->Ns[clabel] <= 0)
            continue;

        cloc = clabel*this->ndim;
        for(j = 0; j < this->ndim; j++)
        {
            (*outStrm)<<this->arrayD[cloc+j]/this->Ns[clabel]<<" ";
        }
        (*outStrm)<<endl;
    }

    outStrm->close();
    cout<<"done\n";
    return ;
}

int XTKSums::fetchCenters(float *centers)
{
    unsigned int clabel = 0, j = 0, idxi = 0, rCNum = 0, cloc = 0;
    assert(centers);

    for(clabel = 0; clabel < this->clnumb; clabel++)
    {
        if(this->Ns[clabel] > 0)
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
        if(this->Ns[clabel] <= 0)
            continue;

        cloc = clabel*this->ndim;
        for(j = 0; j < this->ndim; j++)
        {
            centers[idxi + j] = this->arrayD[cloc+j]/this->Ns[clabel];
        }
        idxi += this->ndim;
    }
    return rCNum;
}

double XTKSums::l2nDst(double *D, const unsigned int x, const unsigned int n)
{
    unsigned int i = 0, loc = 0;
    double dist = 0, delta = 0;
    loc = x*ndim;
    for(i = 0; i < this->ndim; i++)
    {
        delta = (D[i] - n*this->data[loc + i]);
        dist += delta*delta;
    }

    return dist;
}

double XTKSums::crsDst(double Dsc, double *D, const unsigned int x, const unsigned int n)
{
    unsigned int i = 0, loc = 0;
    double dist = 0, xd = 0;
    loc = x*this->ndim;
    for(i = 0; i < this->ndim; i++)
    {
        xd += D[i]*this->data[loc + i];
    }
    dist = (Dsc - 2.0*xd + n*this->l2norms[x]);
    return dist;
}

double XTKSums::innDct(double *D, const unsigned int x, bool _plus_)
{
    unsigned int i = 0, loc = 0;
    double dist = 0, pd = 0, xd = 0;
    loc = x*this->ndim;
    for(i = 0; i < this->ndim; i++)
    {
        pd += D[i]*D[i];
        xd += D[i]*this->data[loc + i];
    }
    if(_plus_)
    {
        dist = pd + 2.0*xd + this->l2norms[x];
    }else{
        dist = pd - 2.0*xd + this->l2norms[x];
    }
    return dist;
}

double XTKSums::cosDst(double *D, const unsigned int x, bool _plus_)
{
    unsigned int i = 0, loc = x*this->ndim;
    double sim = 0, sdct1 = 0, sdct2 = 0, dct = 0, dist = 0;
    for(i = 0; i < this->ndim; i++)
    {
        dct   += D[i]*this->data[loc+i];
        sdct1 += D[i]*D[i];
        sdct2 += this->data[loc+i]*this->data[loc+i];
    }
    if(_plus_)
    {
       sim = (dct+sdct2)/(sqrt(sdct1+2*dct+sdct2)*sqrt(sdct2));
    }else{
       sim = dct/(sqrt(sdct1)*sqrt(sdct2));
    }

    dist = 2.0 - 2.0*sim;
    return dist;
}

double XTKSums::optzI1(const unsigned int nIter, const int nclust)
{
    unsigned int i = 0, iter = 0, j = 0, cloc = 0, sloc = 0, dloc = 0;
    unsigned int s = 0, c = 0, c0 = 0, nwC = 0, n1 = 0, n2 = 0;
    double distA = 0, distB = 0, profit = 0, maxProf = 0;
    double avgDistort = 0, avgDst = 0;
    clock_t tcost = 0;
    int label = 0;

    vector<MiniNN>::iterator vit;
    int *seq = NULL;
    this->distortions = new double[nIter+1];
    this->tcs = new double[nIter+1];
    memset(this->distortions, 0, sizeof(double)*(nIter+1));
    memset(this->tcs, 0, sizeof(double)*(nIter+1));
    memset(this->Ns, 0, sizeof(int)*nclust);

    seq = new int[this->count];
    for(i = 0; i < this->count; i++)
    {
        c = this->labels[i];
        this->Ns[c] += 1;
        seq[i] = i;
    }

    avgDistort = calAvgDistort(nclust, 2);
    this->distortions[0] = avgDistort;
    cout<<"iter\t0\t"<< avgDistort<<endl;
    //this->ed = clock();
    //tcost = (double)(this->ed - this->bg)/(0.0+CLOCKS_PER_SEC);
    for(iter = 0; iter < nIter; iter++)///set nIter to 30 normally
    {
        random_shuffle(seq, seq+this->count);
        for(i = 0; i < this->count; i++)///one round
        {
            s = seq[i];

            c0 = this->labels[s]; ///only one element, no movement
            if(this->Ns[c0] == 1)
                continue;

            cloc = c0*ndim;
            n1 = this->Ns[c0];

            /**/
            distA = l2nDst(this->arrayD+cloc, s, n1);
            distA = distA/(n1*n1);
            ///distA = cosDst(this->arrayD+cloc, s, false);

            maxProf = -1;
            ///search for the best movement
            for(label = 0; label < nclust; label++)
            {
                if(label != (int)c0)
                {
                    cloc  = label*ndim;
                    n2    = this->Ns[label]+1;
                    /**/
                    distB = l2nDst(arrayD+cloc, s, n2-1);
                    distB = distB/(n2*n2);/**/
                    ///distB = cosDst(this->arrayD+cloc, s, true);
                    profit = distA - distB;

                    if(profit > maxProf)
                    {
                        maxProf = profit;
                        nwC     = label;
                    }
                }
            }

            if(maxProf < 0.0f)
            {
                continue;
            }

             ///carry out movement
            this->labels[s] = nwC;
            this->Ns[nwC]  += 1;
            this->Ns[c0]   -= 1;
            this->Ds[nwC] += this->l2norms[s];
            this->Ds[c0]  -= this->l2norms[s];
            cloc = c0*ndim;
            dloc = nwC*ndim;
            sloc = s*ndim;

            for(j = 0; j < this->ndim; j++)
            {
                this->arrayD[cloc+j] -= this->data[sloc + j];
                this->arrayD[dloc+j] += this->data[sloc + j];
            }

        }///end-for(i)
        /**/
        this->ed = clock();
        tcost = (double)(this->ed - this->bg)/(0.0 + CLOCKS_PER_SEC);

        avgDistort = calAvgDistort(nclust, 2);
        avgDst     = pairwDst(nclust);
        this->distortions[iter+1] = avgDistort;
        this->tcs[iter+1] = tcost;
        /**/

        if(1)
          cout<<"iter"<<iter+1<<": "<<avgDistort<<"\t"<<avgDst<<"\t"<<tcost<<"\tsec(s)"<<endl;
        if(this->record[6*iter] > avgDst)
            this->record[6*iter] = avgDst;
        if(this->record[6*iter+1] < avgDst)
            this->record[6*iter+1] = avgDst;
        this->record[6*iter+2] += avgDst;

        if(this->record[6*iter+3] > avgDistort)
            this->record[6*iter+3] = avgDistort;
        if(this->record[6*iter+4] < avgDistort)
            this->record[6*iter+4] = avgDistort;
        this->record[6*iter+5] += avgDistort;

    }///end-for(iter)
    avgDistort = calAvgDistort(nclust, 2);

    delete [] seq;
    seq = NULL;

    return avgDistort;
}


double XTKSums::optzI2(const unsigned int nIter, const int nclust)
{
    unsigned int i = 0, iter = 0, j = 0, cloc = 0, sloc = 0, dloc = 0;
    unsigned int s = 0, c = 0, c0 = 0, nwC = 0, n1 = 0, n2 = 0, nUps = 20;
    double distA = 0, distB = 0, profit = 0, maxProf = 0;
    double avgDistort = 0, avgDst = 0;
    clock_t tcost = 0;
    int label = 0;

    vector<MiniNN>::iterator vit;
    int *seq = NULL;
    this->distortions = new double[nIter+1];
    this->tcs = new double[nIter+1];
    memset(this->distortions, 0, sizeof(double)*(nIter+1));
    memset(this->tcs, 0, sizeof(double)*(nIter+1));
    memset(this->Ns, 0, sizeof(int)*nclust);


    seq = new int[this->count];
    for(i = 0; i < this->count; i++)
    {
        c = this->labels[i];
        this->Ns[c] += 1;
        seq[i] = i;
    }

    avgDistort = calAvgDistort(nclust, 2);

    this->distortions[0] = avgDistort;
    cout<<"iter\t0\t"<< avgDistort<<endl;
    this->ed = clock();
    tcost = (double)(this->ed - this->bg)/(0.0+CLOCKS_PER_SEC);

    for(iter = 0; iter < nIter && nUps > 4; iter++)///set nIter to 30 normally
    {
        random_shuffle(seq, seq+this->count);
        nUps = 0;
        for(i = 0; i < this->count; i++)///one round
        {
            s = seq[i];

            c0 = this->labels[s]; ///only one element, no movement
            if(this->Ns[c0] == 1)
                continue;

            cloc = c0*ndim;
            n1 = this->Ns[c0];

            distA = crsDst(this->Ds[c0], this->arrayD+cloc, s, n1);
            maxProf = -1;

            ///search for the best movement
            for(label = 0; label < nclust; label++)
            {
                if(label != (int)c0)
                {
                    cloc  = label*ndim;
                    n2    = this->Ns[label];

                    distB = crsDst(Ds[label], this->arrayD+cloc, s, n2);
                    profit = distA - distB;

                    if(profit > maxProf)
                    {
                        maxProf = profit;
                        nwC     = label;
                    }
                }
            }

            if(maxProf < 0.0f)
            {
                continue;
            }

             ///carry out movement
            this->labels[s] = nwC;
            this->Ns[nwC]  += 1;
            this->Ns[c0]   -= 1;
            this->Ds[nwC] += this->l2norms[s];
            this->Ds[c0]  -= this->l2norms[s];
            cloc = c0*ndim;
            dloc = nwC*ndim;
            sloc = s*ndim;

            for(j = 0; j < this->ndim; j++)
            {
                this->arrayD[cloc+j] -= this->data[sloc + j];
                this->arrayD[dloc+j] += this->data[sloc + j];
            }
            nUps++;

        }///end-for(i)
        this->ed = clock();
        tcost = (double)(this->ed - this->bg)/(0.0 + CLOCKS_PER_SEC);

        avgDistort = calAvgDistort(nclust, 2);
        avgDst     = pairwDst(nclust);
        this->distortions[iter+1] = avgDistort;
        this->tcs[iter+1] = tcost;

        if(1)
          cout<<"iter\t"<<iter+1<<": "<<avgDistort<<"\t"<<avgDst<<"\t"<<nUps<<"\t"<<tcost<<"\tsec(s)"<<endl;


        if(this->record[6*iter] > avgDst)
            this->record[6*iter] = avgDst;
        if(this->record[6*iter+1] < avgDst)
            this->record[6*iter+1] = avgDst;
        this->record[6*iter+2] += avgDst;

        if(this->record[6*iter+3] > avgDistort)
            this->record[6*iter+3] = avgDistort;
        if(this->record[6*iter+4] < avgDistort)
            this->record[6*iter+4] = avgDistort;
        this->record[6*iter+5] += avgDistort;



    }///end-for(iter)

    delete [] seq;
    seq = NULL;

    return avgDistort;
}


void XTKSums::saveDistort(const char* dstFn)
{
    ofstream *outStrm = new ofstream(dstFn, ios::out);
    unsigned int i = 0;

    for(i = 0; i <= this->nIter; i++)
    {
        (*outStrm)<<i<<"\t"<<this->distortions[i]<<"\t"<<this->tcs[i]<<endl;
    }

    outStrm->close();
}

double XTKSums::pairwDst(const int nclust)
{
    unsigned int i = 0, label = 0, cloc = 0, c0 = 0, nl = 0;
    double l2nsum = 0, pd = 0, *D;
    for(label = 0; label < nclust; label++)
    {
        cloc = label * ndim;
        D    = this->arrayD+cloc;
        for(i = 0; i < this->ndim; i++)
        {
            pd += D[i]*D[i];
        }
    }
    for(i = 0; i < this->count; i++)
    {
        c0 = this->labels[i];
        nl = this->Ns[c0];
        l2nsum += nl * this->l2norms[i];
    }
    return (l2nsum-pd)/this->count;
}

void XTKSums::saveCandlDistort(const char* dstFn)
{
    /**
    ofstream *outStrm = new ofstream(dstFn, ios::out);
    unsigned int i = 0;

    for(i = 0; i <= this->nIter; i++)
    {
        (*outStrm)<<i<<"\t"<<this->candlDsts[0][i]<<"\t";
        (*outStrm)<<this->candlDsts[1][i]<<"\t";
        (*outStrm)<<this->candlDsts[2][i]<<"\n";
    }

    outStrm->close();
    /**/
}

void XTKSums::saveMircoDistort(const char* dstFn)
{
    ofstream *outStrm = new ofstream(dstFn, ios::out);
    unsigned int i = 0, cunt = 0;
    double val1 = 0, val2 = 0, val3 = 0, val4 = 0;
    /**
    val1 = this->dists[0][0];
    val2 = this->dists[1][0];
    val3 = this->dists[2][0];
    val4 = this->dists[3][0];

    (*outStrm)<<i<<"\t"<<val1<<"\t";
    (*outStrm)<<val2<<"\t";
    (*outStrm)<<val3<<"\t";
    (*outStrm)<<val4<<"\n";

    for(i = 1; i < this->count; i++)
    {
       if(this->dists[0][i] > 0)
       {
           val1 = this->dists[0][i];
           (*outStrm)<<i<<"\t"<<this->dists[0][i]<<"\t";
       }else{

           (*outStrm)<<i<<"\t"<<val1<<"\t";
       }
       if(this->dists[1][i] > 0)
       {
           val2 = this->dists[1][i];
           (*outStrm)<<this->dists[1][i]<<"\t";
       }else{

           (*outStrm)<<val2<<"\t";
       }
       if(this->dists[2][i] > 0)
       {
           val3 = this->dists[2][i];
           (*outStrm)<<this->dists[2][i]<<"\t";
       }else{

           (*outStrm)<<val3<<"\t";
       }
       if(this->dists[3][i] > 0)
       {
           val4 = this->dists[3][i];
           (*outStrm)<<this->dists[3][i]<<"\n";
       }else{

           (*outStrm)<<val4<<"\n";
       }
    }
    /**/

    outStrm->close();
}

XTKSums::~XTKSums()
{
    refresh();
}

void XTKSums::test()
{
    const char *srcFn1 = "/home/hye/datasets/clust/sift_learn.txt";
    const char *distortFn1 = "/home/hye/result/sift100k_60t_org_distort.txt";

    const char *srcFn2 = "/home/hye/datasets/clust/sift1m/sift_base.txt";
    const char *distortFn2 = "/home/hye/result/sift1m_60t_org_distort.txt";
    const char *dstFn2 = "";

    const char *srcFn3 = "/home/hye/datasets/glove1m/glove1m_base.txt";
    const char *distortFn3 = "/home/hye/datasets/glove1m/glove1m_60t_org_distort.txt";
    int cNum = 1024;

    XTKSums *mykm = new XTKSums();
    Timer *mytm = new Timer();
    mytm->start();
    mykm->init(srcFn3);
    mykm->clust(cNum, dstFn2, 0);
    mytm->end(true);

    //mykm->saveDistort(distortFn3);
}

void XTKSums::wltest()
{
    const char *srcFn1 = "/home/wlzhao/datasets/bignn/sift1m/sift_learn.txt";
    //const char *srcFn1 = "/home/rqchen/dataset/sift_learn.txt";
    const char *distortFn1 = "/home/wlzhao/datasets/sift100k_60t_org_distort.txt";
    const char *microDstFn = "/home/wlzhao/datasets/clust/distort/micro_distort_xtksum_iter=5.txt";
    ///const char *dstFn1 = "/home/wlzhao/datasets/clust/rslt/sift100k_clust.txt";
    const char *dstFn1 = "/home/rqchen/dataset/rslt/sift100k_clust.txt";

    const char *srcFn2 = "/home/wlzhao/datasets/bignn/sift1m/sift_base.txt";
    ///const char *srcFn2 = "/home/rqchen/dataset/sift_base.txt";
    const char *distortFn2 = "/home/wlzhao/datasets/clust/distort/sift1m_80_xtksums_distort_v7.txt";
    const char *dstFn2 = "/home/wlzhao/datasets/clust/distort/sift1m_ksumsI2_varyk.txt";

    const char *srcFn3 = "/home/wlzhao/datasets/bignn/glove/glove1m_base.txt";
    ///const char *srcFn3 = "/home/rqchen/dataset/glove1m_base.txt";
    const char *distortFn3 = "/home/wlzhao/datasets/clust/distort/glove1m_norm_80_xtksums_distort_v7.txt";
    const char *dstFn3 = "/home/wlzhao/datasets/clust/distort/glove1m_ksumsI2_varyk.txt";

    const char *srcFn4 = "/home/wlzhao/datasets/clust/msd-rh.txt";
    ///const char *srcFn4 = "/home/rqchen/dataset/msd-rh.txt";
    const char *distortFn4 = "/home/wlzhao/datasets/clust/distort/msd-rh_80_xtksums_distort_v7.txt";
    const char *dstFn4 = "/home/wlzhao/datasets/clust/distort/msd1m_ksumsI2_varyk.txt";

    const char *srcFn5 = "/home/wlzhao/datasets/clust/SUSY.csv";
    ///const char *srcFn5 = "/home/rqchen/dataset/SUSY.csv";
    const char *distortFn5 = "/home/wlzhao/datasets/clust/distort/susy_80_xtksums_distort_v7.txt";
    const char *dstFn5 = "/home/wlzhao/datasets/clust/distort/susy1m_ksumsI2_varyk.txt";

    int cnum = 128, i = 0, j = 0, n = 0;
    const char *srcs[5] = {srcFn1, srcFn2, srcFn3, srcFn4, srcFn5};
    const char *dsts[5] = {dstFn1, dstFn2, dstFn3, dstFn4, dstFn5};
    for(j = 0; j < 1; j++)
    {
        cnum = 1024;
        for(i = 0; i < 1; i++)
        {
            XTKSums *mykm = new XTKSums();
            mykm->buildcluster(srcs[j], dsts[j], "rnd", "large", "", cnum, false);
            delete mykm;
            cnum = 2*cnum;
        }
    }

    ///mykm->saveDistort(distortFn3);
    ///mykm->saveMircoDistort(microDstFn);
}

