#include "xbksums.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <cassert>
#include <float.h>
#include <cmath>

#include "vstring.h"
#include "ioagent.h"
#include "cleaner.h"
#include "pqmath.h"
#include "timer.h"
#include "nnitem.h"

const int XBKSums::nTrails   = 10;
const double XBKSums::ERR0   = -0.001f;

XBKSums::XBKSums()
{
    this->_INIT_     = false;
    this->bstLabels  = NULL;
    this->nwlabels   = NULL;
    this->Ds         = NULL;
}

XBKSums::~XBKSums()
{

    this->refresh();
}

bool XBKSums::init(const char *srcfn)
{
    cout<<"Method ........................... bisecting ksum"<<endl;
    cout<<"Loading matrix ................... ";
    this->refresh();
    assert(srcfn);
    strcpy(srcmatfn, srcfn);

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
    if(this->bstLabels == NULL)
    {
        this->bstLabels = new int[this->count];
        memset(this->bstLabels, 0, this->count*sizeof(int));
    }

    if(this->data == NULL)
    {
        cout<<"Exceptions ....................... Loading matrix failed!\n";
        exit(0);
    }

    this->_INIT_ = true;
    this->_REFER_ = false;

    return true;
}

bool XBKSums::init(float *data, const int row, const int dim)
{
    this->refresh();

    this->data  = data;
    this->count = row;
    this->ndim  = dim;

    if(this->bstLabels == NULL)
    {
        this->bstLabels = new int[this->count];
    }
    memset(this->bstLabels, 0, this->count*sizeof(int));

    getL2norms(this->count, this->ndim);

    this->_INIT_  = true;
    this->_REFER_ = true;
    return true;
}

bool XBKSums::config(const char *_seed_, const char *lg_first, const char *crtrn, int verbose)
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


int XBKSums::getL2norms(const unsigned int n, const unsigned int d0)
{
   unsigned int i = 0, j = 0,loc = 0;
   assert(this->data != NULL);
   assert(this->lens != NULL);
   memset(this->lens,   0, sizeof(double)*n);

   for(i = 0; i < n; i++)
   {
        loc = i*d0;
        for(j = 0; j < d0; j++)
        {
            this->lens[i] += this->data[loc+j]*this->data[loc+j];
        }
   }
   return 0;
}

int XBKSums::clust(const unsigned int nclust, const char *dstFn, const int verbose)
{
    double distort = 0, minDst = RAND_MAX + 0.0f;
    unsigned int i = 0, j = 0, k = 0, m = 0;
    unsigned int _loc1 = 0, _loc2 = 0;
    double avgD1 = 0, avgD2 = 0;
    unsigned int n1 = 0, n2 = 0;
    int clabel = 0;
    NNItem *crntItm = NULL;
    vector<NNItem*>::iterator vit;
    bool OPTM = false;

    if(verbose)
    {
        cout<<"Clustering ....................... on progress"<<endl;
        cout<<"k is set to ...................... "<<nclust<<endl;
    }

    memset(this->arrayD, 0 ,sizeof(double)*nclust*this->ndim);

    if(this->Ds == NULL)
    {
        this->Ds = new double[nclust];
    }else{

        delete [] this->Ds;
        this->Ds = new double[nclust];
    }
    memset(this->Ds, 0 ,sizeof(double)*nclust);
    getL2norms(this->count, this->ndim);

    for(i = 0; i < this->nTrails; i++)
    {
        clabel = 0;

        memset(this->labels, 0, sizeof(unsigned int)*this->count);
        memset(this->Ns,     0, sizeof(unsigned int)*nclust);
        memset(this->arrayD, 0, sizeof(double)*nclust*this->ndim);

        _loc1 = clabel*this->ndim;
        for(j = 0; j < this->count; j++)
        {
            _loc2 = j*this->ndim;
            for(k = 0; k < this->ndim; k++)
            {
               this->arrayD[_loc1+k] += this->data[_loc2+k];
            }
        }

        this->Ns[clabel] = this->count;

        crntItm = new NNItem(clabel, 0, this->count);
        this->sorted_stack.push_back(crntItm);

        ///produce the other cluster
        j = 1;
        while(j < nclust)
        {
            if(verbose)
            {
                cout<<"begin to incrOptz "<<j<<"..............."<<endl;
            }

            OPTM = bisect(clabel, j, n1, n2);

            if(OPTM)
            {
                crntItm        = sorted_stack[0];
                crntItm->index = clabel;
                crntItm->size  = this->Ns[clabel];
                crntItm->val   = avgD1;

                crntItm = new NNItem(j, avgD2, this->Ns[j]);
                sorted_stack.push_back(crntItm);
                j++;
            }
            else
            {
                crntItm->dvd    = false;
                cout<<"unseparable!"<<endl;
            }

            stable_sort(sorted_stack.begin(), sorted_stack.end(), NNItem::LGSZcomparer);

            ///find the next cluster ready to be separated
            m = 0;
            for(vit = sorted_stack.begin(); vit != sorted_stack.end(); vit++,m++)
            {
                crntItm = *vit;
                if(crntItm->size > 1 && crntItm->dvd)
                    break;
            }

            if(crntItm->size == 1 || !crntItm->dvd)
            {
                cout<<"Failed to reach expected cluster number!\n";
                break;   ///break-while
            }
            else
            {
                clabel = crntItm->index;
            }
        } ///end-while
        refineI2(nclust, 3);

        distort = calAvgDistort(nclust);
        ///distort = pairwDst(nclust);
        if(minDst > distort)
        {
             minDst = distort;
             memcpy(this->bstLabels, this->labels, this->count*sizeof(unsigned int));
        }

        if(1)
        {
            cout<<"nTrails: "<<i<<": "<<distort<<" minDst: "<<minDst<<endl;
        }

        Cleaner::clearVector(sorted_stack);
    } ///end-for

    memcpy(this->labels, this->bstLabels, this->count*sizeof(unsigned int));

    this->save_clust(dstFn);

    return 0;
}

bool XBKSums::bisect(const unsigned int odlbl, const unsigned int nwlbl, unsigned int &n1, unsigned int &n2)
{
    unsigned int x = 0, j = 0, i = 0, csize = 0;
    unsigned int _loc1 = 0, _loc2 = 0, _loc3 = 0;
    double dst1 = 0, dst2 = 0;

    this->crntClust.reserve(this->count);
    for(i = 0; i < this->count; i++)
    {
        if((int)odlbl == this->labels[i])
        {
            this->crntClust.push_back(i);
        }
    }

    csize = this->crntClust.size();
    if(csize == 2)
    {
        n1 = n2 = 0;
        this->Ns[odlbl] = 1;
        this->Ns[nwlbl] = 1;

        i = crntClust[1];
        this->labels[i] = nwlbl;

        _loc1 = this->crntClust[0]*this->ndim;
        _loc2 = this->crntClust[1]*this->ndim;
        for(i = 0; i < this->ndim; i++)
        {
            this->arrayD[odlbl*this->ndim+i] = this->data[_loc1+i];
            this->arrayD[nwlbl*this->ndim+i] = this->data[_loc2+i];
        }

        this->crntClust.clear();
        return true;
    }

    /****** Initial assignment ******/
    if(this->nwlabels != NULL)
        delete [] this->nwlabels;

    this->nwlabels = new int[csize];
    memset(this->nwlabels, -1, sizeof(int)*csize);

    _loc1 = odlbl*this->ndim;
    _loc2 = nwlbl*this->ndim;
    csize = this->crntClust.size();

    this->Ns[nwlbl] = 0;
    this->Ds[nwlbl] = 0;
    this->Ds[odlbl] = 0;
    for(i = 0; i < csize ; i++)
    {
        x = this->crntClust[i];
        dst1 = rand();
        dst2 = rand();
        if(dst1 < dst2)
        {

            this->Ns[nwlbl]++;
            this->Ns[odlbl]--;
            this->nwlabels[i] = x;
            _loc3 = x*this->ndim;
            this->Ds[nwlbl] += this->lens[x];
            for(j = 0; j < this->ndim; j++)
            {
               this->arrayD[_loc1+j] -= this->data[_loc3+j];
               this->arrayD[_loc2+j] += this->data[_loc3+j];
            }
        }
        else
        {
            this->Ds[odlbl]  += this->lens[x];
            this->nwlabels[i] = -1;
        }
    }

    /****** Incremental optimization ******/
    optzI2(odlbl, nwlbl);

    /**
    double distort = 0;
    distort = calAvgDistort(nwlbl+1);
    cout<<"cluster:\t"<<nwlbl+1<<"\t\t"<<distort<<endl;
    /**/

    this->crntClust.clear();
    delete [] this->nwlabels;
    this->nwlabels = NULL;
    return true;
}

bool XBKSums::optzI1(const unsigned int odlbl, const unsigned int nwlbl)
{
    double distA = 0, distB = 0;
    unsigned int sz = 0, _loc1 = 0, _loc2 = 0, _loc3 = 0;
    unsigned int x = 0, c0 = 0, c1 = 0, index = 0;
    unsigned int iter = 0, i = 0, j = 0, nUps = 0;
    double profit = 0, delta = 0;
    int *seq = NULL;

    sz  = this->crntClust.size();
    seq = new int[sz];
    memset(seq, 0, sizeof(int)*sz);

    for(i = 0; i < sz; i++)
    {
        x = this->crntClust[i];
    }

    do
    {
        for(i = 0; i < sz; i++)
        {
            seq[i] = i;
        }
        random_shuffle(seq, seq+sz);
        nUps = 0;
        ///sumDst = 0;
        delta = 0;
        for(i = 0; i < sz; i++)
        {
            index = seq[i];
            x     = this->crntClust[index];

            if(this->nwlabels[index] == -1)
            {
                c0 = odlbl;
                c1 = nwlbl;
            }
            else
            {
                c0 = nwlbl;
                c1 = odlbl;
            }

            if(this->Ns[c0] <= 1)
                continue;

            _loc1 = c0 * this->ndim;
            _loc2 = c1 * this->ndim;
            _loc3 = x  * this->ndim;
            /**
            distA = calDist(this->arrayD+_loc1, x, nA+1);
            distB = calDist(this->arrayD+_loc2, x, nB-1);
            distA = distA/nA; distA = distA/nA; //otherwise, overflow happens
            distB = distB/nB; distB = distB/nB; //otherwise, overflow happens
            /**/

            /**/
            distA = cosDst(this->arrayD+_loc1, x, false);
            distB = cosDst(this->arrayD+_loc2, x, true);
            /**/

            if(distB < distA)
            {
                if(c0 == odlbl)
                {
                    this->nwlabels[index] = x;
                }
                else
                {
                    this->nwlabels[index] = -1;
                }

                this->Ns[c0]--;
                this->Ns[c1]++;

                for(j = 0; j < this->ndim; j++)
                {
                    this->arrayD[_loc1 + j] -= this->data[_loc3 + j];
                    this->arrayD[_loc2 + j] += this->data[_loc3 + j];
                }
                nUps++;
                profit = distA - distB;
                delta += profit;
            }
        }///end-for(i)

        iter++;
    }while(nUps > 8 && iter < 64); ///end-while(_UPDATE)
    cout<<"nUps:\t"<<nUps<<"\t"<<iter<<endl;

    ///save the result in labels and members
    for(i = 0; i < sz; i++)
    {
        x = this->crntClust[i];
        if(this->nwlabels[i] != -1)
        {
            this->labels[x] = nwlbl;
        }
    }

    delete [] seq;
    seq = NULL;

    return true;
}

bool XBKSums::optzI2(const unsigned int odlbl, const unsigned int nwlbl)
{
    double distA = 0, distB = 0;
    unsigned int sz = 0, _loc1 = 0, _loc2 = 0, _loc3 = 0;
    unsigned int x = 0, c0 = 0, c1 = 0, index = 0;
    unsigned int iter = 0, nA = 0, nB = 0;
    unsigned int i = 0, j = 0, nUps = 0;
    int *seq = NULL;
    sz  = this->crntClust.size();
    seq = new int[sz];
    memset(seq, 0, sizeof(int)*sz);

    for(i = 0; i < sz; i++)
    {
        x = this->crntClust[i];
    }

    do
    {
        for(i = 0; i < sz; i++)
        {
            seq[i] = i;
        }
        random_shuffle(seq, seq+sz);
        nUps = 0;
        for(i = 0; i < sz; i++)
        {
            index = seq[i];
            x     = this->crntClust[index];

            if(this->nwlabels[index] == -1)
            {
                c0 = odlbl;
                c1 = nwlbl;
            }
            else
            {
                c0 = nwlbl;
                c1 = odlbl;
            }

            if(this->Ns[c0] <= 1)
                continue;

            _loc1 = c0 * this->ndim;
            _loc2 = c1 * this->ndim;
            _loc3 = x  * this->ndim;
            nA    = this->Ns[c0];
            nB    = this->Ns[c1];

            distA = this->crsDst(this->Ds[c0], this->arrayD+_loc1, x, nA);
            distB = this->crsDst(this->Ds[c1], this->arrayD+_loc2, x, nB);

            if(distB < distA)
            {
                if(c0 == odlbl)
                {
                    this->nwlabels[index] = x;
                }
                else
                {
                    this->nwlabels[index] = -1;
                }

                this->Ns[c0]--;
                this->Ns[c1]++;
                this->Ds[c0] -= this->lens[x];
                this->Ds[c1] += this->lens[x];

                for(j = 0; j < this->ndim; j++)
                {
                    this->arrayD[_loc1 + j] -= this->data[_loc3 + j];
                    this->arrayD[_loc2 + j] += this->data[_loc3 + j];
                }

                nUps++;
            }
        }///end-for(i)
        iter++;
    }while(nUps > 1 && iter < 64);
    cout<<"nUps:\t"<<nUps<<"\t"<<iter<<endl;
    ///save the result in labels and members
    for(i = 0; i < sz; i++)
    {
        x = this->crntClust[i];
        if(this->nwlabels[i] != -1)
        {
            this->labels[x] = nwlbl;
        }
    }

    delete [] seq;
    seq = NULL;

    return true;
}

double XBKSums::refineI2(const int nclust, const unsigned int nIter)
{
   double distA = 0, distB = 0, maxProf = 0, avgDistort = 0, profit = 0, avgDst = 0;
   unsigned int cloc = 0, sloc = 0, n1 = 0, n2 = 0, dloc = 0;
   unsigned int i = 0, nUps = 0, iter = 0, s = 0, j = 0;
   int c = 0, c0 = 0, nwC = 0, label = 0;

   memset(this->Ds, 0, sizeof(double)*nclust);
   memset(this->Ns, 0, sizeof(unsigned)*nclust);
   unsigned int *seq = new unsigned int[this->count];
   for(i = 0; i < this->count; i++)
   {
        c = this->labels[i];
        this->Ds[c] += this->lens[i];
        seq[i] = i;
        this->Ns[c] += 1;
   }
   cout<<"Perform refinement ................. \n";

   for(iter = 0; iter < nIter; iter++)
    {
        random_shuffle(seq, seq+this->count);
        nUps = 0;
        for(i = 0; i < this->count; i++)///one round
        {
            s = seq[i];

            c0 = this->labels[s];
            if(this->Ns[c0] == 1)
                continue;

            cloc = c0*ndim;
            n1 = this->Ns[c0];

            /**/
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
            this->Ds[nwC] += this->lens[s];
            this->Ds[c0]  -= this->lens[s];
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

        avgDistort = calAvgDistort(nclust);
        avgDst     = pairwDst(nclust);

        if(1)
          cout<<"iter\t"<<iter+1<<": "<<avgDistort<<"\t"<<avgDst<<"\t"<<nUps<<endl;
    }
    delete [] seq;

    return avgDistort;
}

double XBKSums::calDist(double *D, const unsigned int x, const unsigned int n)
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

double XBKSums::cosDst(double *D, const unsigned int x, bool _plus_)
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
       sim = (dct)/(sqrt(sdct1)*sqrt(sdct2));
    }

    dist = 2.0 - 2.0*sim;
    return dist;
}

double XBKSums::crsDst(double Dsc, double *D, const unsigned int x, const unsigned int n)
{
    unsigned int i = 0, loc = 0;
    double dist = 0, xd = 0;
    loc = x*this->ndim;
    for(i = 0; i < this->ndim; i++)
    {
        xd += D[i]*this->data[loc + i];
    }
    dist = (Dsc - 2.0*xd + n*this->lens[x]);
    return dist;
}


double XBKSums::calAvgDistort(const unsigned int nclust)
{
    double *centers = NULL;
    double distort = 0;
    centers = new double[nclust * this->ndim];
    memset(centers, 0, sizeof(double)* nclust * this->ndim);

    unsigned int i = 0, j = 0, c = 0, dloc = 0, cloc = 0;
    double avgDistort = 0;

    for(c = 0; c < nclust; c++)
    {
        cloc = this->ndim*c;
        for(j = 0; j < this->ndim; j++)
        {
           centers[cloc + j] = this->arrayD[cloc+j]/this->Ns[c];
        }
    }

    for(i = 0; i < this->count; i++)
    {
        c = this->labels[i];
        cloc = c*this->ndim;
        dloc = i*this->ndim;
        for(j = 0; j < this->ndim; j++)
        {
            distort = abs(centers[cloc + j] - this->data[dloc + j]);
            avgDistort += distort * distort;
        }
    }
    ///caculate avg distortion
    //cout<<"before: "<<avgDistort<<endl;
    avgDistort /= this->count;
    //cout<<"hi: "<<avgDistort<<"\t"<<this->count<<endl;

    delete[] centers;
    centers = NULL;

    return avgDistort;
}

bool XBKSums::refresh()
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

    if(this->nwlabels != NULL)
    {
        delete [] this->nwlabels;
        this->nwlabels = NULL;
    }

    if(this->bstLabels != NULL)
    {
        delete [] this->bstLabels;
        this->bstLabels = NULL;
    }


    if(this->Ds != NULL)
    {
        delete [] this->Ds;
        this->Ds = NULL;
    }

    return true;
}

void XBKSums::saveCenters(const char *dstfn, bool append)
{
    unsigned int clabel = 0, j = 0, rCNum = 0;
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

        for(j = 0; j < this->ndim; j++)
        {
            (*outStrm)<<this->arrayD[clabel*this->ndim+j]/this->Ns[clabel]<<" ";
        }
        (*outStrm)<<endl;
    }

    outStrm->close();
    cout<<"done\n";
    return ;
}

int XBKSums::fetchCenters(float *centers)
{
    unsigned int clabel = 0, j = 0, idxi = 0, rCNum = 0;
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

        for(j = 0; j < this->ndim; j++)
        {
            centers[idxi + j] = this->arrayD[idxi + j]/this->Ns[clabel];
        }

        idxi += this->ndim;
    }
    return rCNum;
}

double XBKSums::pairwDst(const int nclust)
{
    unsigned int i = 0, label = 0, cloc = 0, c0 = 0, nl = 0;
    double l2nsum = 0, pd = 0, *D = NULL;
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
        l2nsum += nl * this->lens[i];
    }
    return (l2nsum-pd)/this->count;
}


void XBKSums::test()
{
    const char *srcFn1 = "/home/wlzhao/datasets/bignn/sift1m/sift_learn.txt";
    const char *srcFn2 = "/home/wlzhao/datasets/bignn/sift1m/sift_base.txt";
    const char *dstFn1 = "/home/wlzhao/datasets/clust/sift100k_xbksums.txt";
    const char *dstFn2 = "/home/wlzhao/datasets/clust/sift1m_xbksums.txt";

    XBKSums *mykm = new XBKSums();
    Timer *mytm = new Timer();
    mytm->start();
    mykm->buildcluster(srcFn2, dstFn2, "non", "large", "i2", 10000, false);
    mytm->end(true);
    return;
}
