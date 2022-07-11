#include "pqclust.h"

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
#include <limits>
#include <cassert>
#include <cstdlib>
#include <thread>
#include <bitset>
#include <vector>
#include <cmath>
#include <ctime>

using namespace std;

const unsigned int PQClust::M_ITER = 20;
const unsigned int PQClust::NTRAILS= 1;
const unsigned int PQClust::NRfn   = 10;
const float PQClust::EPS  = 0.5f;
const float PQClust::Err0 = 0;
const int interval = 1;
const char *PQClust::pqfn = "/home/chenghaod/dataset/sift1m_pq/sift_pq256seg8.txt";


PQClust::PQClust()
{
    arrayD = Ds = bstArrayD = Cs = NULL;
    bstLabels = NULL;
    _INIT_ = false;
    ///kmMtd  = _pqc_;
    strcpy(mthStr, "_pqc_");
    cout<<"Method ........................... PQ clusting K-means\n";
}

bool PQClust::init(const char *srcfn)
{
    cout<<"Loading matrix ................... ";
    assert(srcfn);
    strcpy(srcmatfn, srcfn);
    refresh();

    data = IOAgent::loadPQ(srcfn, this->count, this->pqm);
    float *PQInfo = IOAgent::loadPQInfo(pqfn, pqn, pql, pqm);
    cout<<this->count<<"x"<<this->pqm<<endl;
    this->ndim = pqm*pql;
    this->pqdis = new double[pqn*pqn*pqm];
    memset(pqdis, 0, sizeof(double)*pqn*pqn*pqm);

    if(this->data == NULL)
    {
        cout<<"Exceptions ....................... Loading matrix failed!\n";
        exit(0);
    }

    unsigned int i = 0, j = 0, k = 0, l = 0, m = 0;
    double dis = 0;

    ///get l2 distence
    for(i = 0; i < pqm; i++)
    {
        for(j = 0; j < pqn; j++)
        {
            for(k = 0; k < pqn; k++)
            {
                if(j == k)
                    continue;
                dis = 0;
                for(m = 0; m < pql; m++)
                {
                    dis += (PQInfo[i*pqn*pql + j*pql + m]-PQInfo[i*pqn*pql + k*pql + m])*(PQInfo[i*pqn*pql + j*pql + m]-PQInfo[i*pqn*pql + k*pql + m]);
                }
                /// dis = sqrt(dis);
                this->pqdis[pqn*pqn*i+j*pqn+k] = dis;
            }
        }
    }

    for(i = 0; i < pqm; i++)
    {
        for(j = 0; j < pqn; j++)
        {
            this->pqmap.insert(pair<unsigned, float *>(l, PQInfo+l*pql));
            l++;
        }
    }


    this->Ds = new double[this->clnumb*this->ndim];
    this->Cs = new double[this->clnumb*this->ndim];
    this->arrayD = new double[this->clnumb*this->ndim];
    this->labels = new int[this->count];


    this->bstLabels = new int[this->count];
    memset(this->bstLabels, 0 , sizeof(int)*this->count);

    this->_INIT_  = true;
    this->_REFER_ = false;

    cout<<"clust num\t"<<this->clnumb<<endl;
    return true;
}

bool PQClust::init(float *mat, const int row, const int dim)
{
}

bool PQClust::init(unsigned *mat, const int row, const int dim)
{
    this->refresh();

    this->data  = mat;
    this->count = row;
    this->ndim  = dim;

    this->bstLabels = new int[this->count];
    memset(this->bstLabels, 0 , sizeof(int)*this->count);

    ///memset(this->visited, 0, this->nThrd*sizeof(unsigned char));

    this->_INIT_  = true;
    this->_REFER_ = true;
    return true;
}

bool PQClust::refresh()
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

    if(this->Ds != NULL)
    {
        delete [] this->Ds;
        this->Ds = NULL;
    }

    if(this->Cs != NULL)
    {
        delete [] this->Cs;
        this->Cs = NULL;
    }
    if(this->tmpCs != NULL)
    {
        delete [] this->tmpCs;
        this->tmpCs = NULL;
    }


    ///memset(this->visited, 0, nThrd*sizeof(unsigned char));
    return true;
}

bool PQClust::config(const char *_seed_, const char *crtrn, const char *lg_first, int verbose)
{
    if(verbose)
        cout<<"Distance function ................ l2\n";

    if(verbose)
        cout<<"Seeds ............................ ";
    if(!strcmp(_seed_, "rnd"))
    {
        if(verbose)
            cout<<"rand\n";
        seed = _rnd_;
    }
    else if(!strcmp(_seed_, "kpp"))
    {

        if(verbose)
            cout<<"kpp\n";
        seed = _kpp_;
    }
    else
    {
        if(verbose)
            cout<<"kpp\n";
        seed = _kpp_;
    }

    if(verbose)
        cout<<"Optimization function ............ ";
    if(!strcmp(crtrn, "i1"))
    {
        myoptz     = _I1_;
        if(verbose)
            cout<<"I1\n";
    }
    else if(!strcmp(crtrn, "i2"))
    {
        myoptz  = _I2_;
        if(verbose)
            cout<<"I2\n";
    }
    else if(!strcmp(crtrn, "i3"))
    {
        myoptz  = _I3_;
        if(verbose)
            cout<<"I3\n";
    }
    else if(!strcmp(crtrn, "i4") )
    {
        myoptz  = _I4_;
        if(verbose)
            cout<<"e1\n";
    }
    else if(!strcmp(crtrn, "e1") )
    {
        myoptz  = _E1_;
        if(verbose)
            cout<<"e1\n";
    }
    else if(!strcmp(crtrn, "e2") )
    {
        myoptz  = _E2_;
        if(verbose)
            cout<<"e2\n";
    }
    else if(!strcmp(crtrn, "t1"))
    {
        myoptz = _T1_;
        if(verbose)
            cout<<"t1\n";
    }
    else if(!strcmp(crtrn, "i4"))
    {
        myoptz = _I4_;
        if(verbose)
            cout<<"i4\n";
    }
    else
    {
        cout<<"Unkown optimize option '"<<crtrn<<"'!\n";
        exit(0);
    }

    return true;
}

double PQClust::PQDist(unsigned * pq1, unsigned *pq2)
{
    double dis = 0;
    unsigned int i = 0, loc, n;
    n = pqn*pqn;
    loc = 0;
    for(i = 0; i < pqm; i++)
    {
        dis += pqdis[loc + pqn*pq1[i] + pq2[i]];
        loc += n;
    }
    return dis;
}

int PQClust::clust(const unsigned int clust_num, const char *dstfn, const int verbose)
{
///lens
    int unsigned i, j, k, n;
    memset(lens, 0, sizeof(double)*count);
    memset(labels, 0, sizeof(int)*count);
    for(i = 0; i < count; i++)
    {
        for(j = 0; j < pqm; j++)
        {
            for(k = 0; k < pql; k++)
            {
                lens[i] += pqmap[data[i*pqm+j] + j*pqn][k]*pqmap[data[i*pqm+j] + j*pqn][k];
            }
        }
    }
    if(myoptz == _I2_)
    {
        n = this->I2(clnumb, dstfn, verbose);
    }
    else if(myoptz == _T1_)
        n = this->tkm(clnumb, dstfn, verbose);

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
        Print2Scrn::printvect(sorted_stack, infofn);
        PQClust::save_clust(dstfn);
    }

    Cleaner::clearVector(sorted_stack);
/**
    for(i = 0; i < clnumb; i++)
    {
        cout<<"cluster num "<<i<<" "<<Ns[i]<<endl;
    }
/**/
    return true;
}

int PQClust::I2(const unsigned int clust_num, const char *dstfn, const int verbose)
{
    double dst = 0, minDst = 0.0f, sumDst0 = 0, err = 0.0f, Eg1 = 0, Eg2 = 0, delta = 0;
    double sumDst = 0, bstSumDst = numeric_limits<double>::max(), optEg = 0, prvEg = 0;
    unsigned int loc = 0, d = 0, niter = 0, loci;
    unsigned int i = 0, j = 0, clabel = 0, id = 0, k = 0, ri = 0;
    unsigned int *PQC;
    PQC = new unsigned[clnumb*pqm];
    double *PD = new double[pqn*pqn*pqm];
    double *Es = new double[clnumb];
    id = 0;
    for(i = 0; i < pqm; i++)
    {
        for(j = 0; j < pqn; j++)
        {
            for(k = 0; k < pqn; k++)
            {
                dst = 0;
                for(d = 0; d < pql; d++)
                {
                    dst += pqmap[i*pqn+j][d]*pqmap[i*pqn+k][d];
                }
                PD[id] = dst;
                id++;
            }
        }
    }

    double en = 0;

    for(i = 0; i < count; i++)
    {
        en += this->lens[i];
    }

    if(!this->_INIT_)
    {
        cout<<"Error ........................... ahead of init!\n";
        return 0;
    }

    int cpysize = this->ndim*clnumb*sizeof(double);

    bstArrayD = new double[this->ndim*clnumb];
    memset(Ds, 0, cpysize);
    memset(bstArrayD, 0, cpysize);

    int    *rseeds = new int[clust_num];
    int    *counts = new int[clust_num];

    for(ri = 0; ri < NTRAILS; ri++)
    {
        memset(counts, 0, sizeof(int)*clust_num);
        memset(Ds, 0, sizeof(double)*clnumb*ndim);
        if(this->seed != _non_)
        {
            if(this->seed == _rnd_)
                this->rndSeeds(clnumb, rseeds, count);
            else if(this->seed == _kpp_)
                this->kppSeeds(clnumb, rseeds, count);
            for(j = 0; j < clust_num; j++)
            {
                cmps++;
                loc = rseeds[j];
                memcpy(PQC + j*pqm, data + loc*pqm, sizeof(unsigned)*pqm);
            }
            for(i = 0; i < count; i++)
            {
                minDst = numeric_limits<double>::max();
                for(k = 0; k < clnumb; k++)
                {
                    dst  = this->PQDist(PQC+pqm*k, data + pqm*i);
                    if(dst < minDst)
                    {
                        clabel = k;
                        minDst = dst;
                    }
                    cmps++;
                }

                labels[i] = clabel;
                counts[clabel] += 1;
                loci = clabel*ndim;
                for(j = 0; j < pqm; j++)
                {
                    loc = j*pql;
                    for(k = 0; k < pql; k++)
                    {
                        Ds[loci + loc + k] += pqmap[j*pqn + data[i*pqm+j]][k];
                    }
                }
            }///for(i)

            for(i = 0; i < clnumb; i++)
            {
                id = rseeds[i];
                if(labels[id] != i)
                {
                    clabel = labels[id];
                    counts[clabel]--;
                    loci = clabel*ndim;
                    for(j = 0; j < pqm; j++)
                    {
                        loc = j*pql;
                        for(k = 0; k < pql; k++)
                        {
                            Ds[loci + loc + k] -= pqmap[j*pqn + data[id*pqm+j]][k];
                        }
                    }
                    clabel = i;
                    counts[clabel]++;
                    loci = clabel*ndim;
                    for(j = 0; j < pqm; j++)
                    {
                        loc = j*pql;
                        for(k = 0; k < pql; k++)
                        {
                            Ds[loci + loc + k] += pqmap[j*pqn + data[id*pqm+j]][k];
                        }
                    }
                }
            }
        }///if !non
        else
        {
            for(i = 0; i <  count; i++)
            {
                clabel = rand()%clnumb;
                this->labels[i] = clabel;
                counts[clabel] += 1;
                loci = clabel*ndim;
                for(j = 0; j < pqm; j++)
                {
                    loc = j*pql;
                    for(k = 0; k < pql; k++)
                    {
                        Ds[loci + loc + k] += pqmap[j*pqn + data[i*pqm+j]][k];
                    }
                }
            }
        }///if else seeds

        ///get cs
        for(i = 0; i < clnumb; i++)
        {
            loc = i*ndim;
            if(counts[i] != 0)
                for(j = 0; j < ndim; j++)
                {
                    Cs[loc + j] = Ds[loc + j]/counts[i];
                }
            else
                for(j = 0; j < ndim; j++)
                {
                    Cs[loc + j] = 0;
                }
        }
/// cs to pqcode
        for(unsigned cluster = 0; cluster < clnumb; cluster++)
        {
            loci = cluster*ndim;
            for(i = 0; i < pqm; i++)
            {
                minDst = numeric_limits<double>::max();
                loc = i*pql;
                for(k = 0; k < pqn; k++)
                {
                    dst = 0;
                    for(j = 0; j < pql; j++)
                    {
                        dst += (Cs[loci + loc + j] - pqmap[i*pqn + k][j])*(Cs[loci + loc + j] - pqmap[i*pqn + k][j]);
                    }
                    if(dst < minDst)
                    {
                        minDst = dst;
                        id = k;
                    }
                }
                PQC[cluster*pqm+i] = id;
            }
        }
///updata es
        for(i = 0; i < clnumb; i++)
        {
            Es[i] = 0;
            loc = i*ndim;
            for(j = 0; j < ndim; j++)
            {
                Es[i] += Ds[loc+j]*Ds[loc+j];
            }
        }

        optEg = getI2(this->Ds, clust_num, ndim, counts);
        prvEg = optEg;

///update center: arrayD
        niter  = 0;
        bool UPDATED = false;
        memcpy(Ns, counts, sizeof(unsigned)*clnumb);
        do
        {
            UPDATED = false;
            ///compute data to arrayD
            memset(Ds, 0, sizeof(double)*ndim*clnumb);
            for(i = 0; i < count; i++)
            {
                // cout<<"bug3\n";
                minDst = 0;
                loc = i*pqm;
                id = labels[i];
                //  cout<<"id\t"<<id<<endl;
                if(Ns[id] <= 1)
                    continue;
                Eg1 = 0;
                for(j = 0; j < pqm; j++)
                {
                    Eg1 += PD[pqn*pqn*j + data[loc + j]*pqn + PQC[id*pqm+j]];
                }
                Eg1 = (Es[id] -2*Eg1*counts[id] + lens[i])/(counts[id]-1) - Es[id]/counts[id];
                //  cout<<"bug3.1\n";
                clabel = id;
                for(k = 0; k < clust_num; k++)
                {
                    if(k == id)
                        continue;
                    Eg2 = 0;

                    for(j = 0; j < pqm; j++)
                    {
                        Eg2 += PD[pqn*pqn*j + data[loc + j]*pqn + PQC[k*pqm+j]];
                    }

                    if(Ns[k] == 0)
                    {
                        Eg2 = 0;
                        cout<<"iter "<<niter<<" counts "<<k<<" = 0\n";
                    }
                    else
                        Eg2 = (Es[k] + Eg2*2*counts[k] + lens[i])/(counts[k]+1) - Es[k]/counts[k];
                    delta = Eg2 + Eg1;
                    cmps++;
                    if(delta > minDst)
                    {
                        clabel = k;
                        minDst = delta;
                        UPDATED = true;
                    }
                }
                // cout<<"bug3.2\n";
                if(labels[i] != clabel)
                {
                    Ns[labels[i]] --;
                    Ns[clabel] ++;
                    labels[i] = clabel;
                }
                loci = clabel*ndim;
                for(j = 0; j < pqm; j++)
                {
                    loc = j*pql;
                    for(k = 0; k < pql; k++)
                    {
                        Ds[loci + loc + k] += pqmap[j*pqn + data[i*pqm+j]][k];
                    }
                }
            }
            memcpy(counts, Ns, sizeof(unsigned)*clnumb);
            ///get cs
            for(i = 0; i < clnumb; i++)
            {
                loc = i*ndim;
                if(counts[i] != 0)
                    for(j = 0; j < ndim; j++)
                    {
                        Cs[loc + j] = Ds[loc + j]/counts[i];
                    }
                else
                    for(j = 0; j < ndim; j++)
                    {
                        Cs[loc + j] = 0;
                    }
            }
/// cs to pqcode
            for(unsigned cluster = 0; cluster < clnumb; cluster++)
            {
                loci = cluster*ndim;
                for(i = 0; i < pqm; i++)
                {
                    minDst = numeric_limits<double>::max();
                    loc = i*pql;
                    for(k = 0; k < pqn; k++)
                    {
                        dst = 0;
                        for(j = 0; j < pql; j++)
                        {
                            dst += (Cs[loci + loc + j] - pqmap[i*pqn + k][j])*(Cs[loci + loc + j] - pqmap[i*pqn + k][j]);
                        }
                        if(dst < minDst)
                        {
                            minDst = dst;
                            id = k;
                        }
                    }
                    PQC[cluster*pqm+i] = id;
                }
            }
///updata es
            for(i = 0; i < clnumb; i++)
            {
                Es[i] = 0;
                loc = i*ndim;
                for(j = 0; j < ndim; j++)
                {
                    Es[i] += Ds[loc+j]*Ds[loc+j];
                }
            }
            niter++;

            optEg = getI2(this->Ds, clust_num, ndim, counts);
            double tmperr = optEg - prvEg;
            prvEg = optEg;
            err = fabs(sumDst - sumDst0);
            sumDst0 = sumDst;
            cout<<niter<<"\t"<<optEg<<"\tdistortion\t"<<(en - optEg)/this->count<<endl;
        }
        while(UPDATED == true && niter < PQClust::M_ITER);
        if(sumDst0 < bstSumDst)
        {
            memcpy(arrayD, Ds, this->ndim*clust_num*sizeof(double));
            memcpy(bstLabels, labels, count*sizeof(int));
            bstSumDst = sumDst0;
            for(j = 0; j < clust_num; j++)
            {
                this->Ns[j] = counts[j];
            }
        }
    }///for(ri)

    memcpy(labels, bstLabels, count*sizeof(int));

    delete [] Es;
    delete [] PD;
    delete [] PQC;
    delete [] bstLabels;
    delete [] bstArrayD;
    delete [] counts;
    delete [] rseeds;

    bstLabels = NULL;
    bstArrayD = NULL;
    rseeds = NULL;
    counts = NULL;

    return clust_num;
}



int PQClust::tkm(const unsigned int clust_num, const char *dstfn, const int verbose)
{
    double dst = 0, minDst = 0.0f, sumDst0 = 0, err = 0.0f;
    double sumDst = 0, bstSumDst = numeric_limits<double>::max() + 0.0f, optEg = 0, prvEg = 0;
    unsigned int loc = 0, d = 0, di = 0, niter = 0, loci;
    unsigned int i = 0, j = 0, clabel = 0, id = 0;
    unsigned int k = 0, ri = 0;
    unsigned n = pqn*pqn;
    unsigned int *PQC;
    PQC = new unsigned[clnumb*pqm];
    //unsigned int *Cs;
    // Cs = new unsigned[pqn*pqm];


    double en = 0;

    for(i = 0; i < count; i++)
    {
        en += this->lens[i];
    }

    if(!this->_INIT_)
    {
        cout<<"Error ........................... ahead of init!\n";
        return 0;
    }

    int cpysize = this->ndim*clust_num*sizeof(double);

    bstArrayD = new double[this->ndim*clust_num];
    memset(Ds, 0, cpysize);
    memset(bstArrayD, 0, cpysize);

    int    *rseeds = new int[clust_num];
    int    *counts = new int[clust_num];


    for(ri = 0; ri < NTRAILS; ri++)
    {
        memset(counts, 0, sizeof(int)*clust_num);

        if(this->seed == _rnd_)
            this->rndSeeds(clnumb, rseeds, count);
        else if(this->seed == _kpp_)
            this->kppSeeds(clnumb, rseeds, count);
        for(j = 0; j < clust_num; j++)
        {
            cmps++;
            loc = rseeds[j];
            memcpy(PQC + j*pqm, data + loc*pqm, sizeof(unsigned)*pqm);
        }


///update center and assignment

        niter  = 0;
        bool UPDATED = false;

        do
        {
            UPDATED = false;
            memset(counts, 0, sizeof(unsigned)*clnumb);
            memset(Ds, 0, sizeof(double)*clnumb*ndim);
            for(i = 0; i < count; i++)
            {
                minDst = numeric_limits<double>::max();
                for(k = 0; k < clnumb; k++)
                {
                    ///dst = this->PQDist(PQC+k*pqm, data+pqm*i);
                    /**/
                    id = 0;
                    dst = 0;
                    loc = k*pqm;
                    loci = pqm*i;
                    for(j = 0; j < pqm; j++)
                    {
                        dst += pqdis[id + pqn*PQC[loc+j] + data[loci+j]];
                        id += n;
                    }
                    /**/
                    if(dst < minDst)
                    {
                        clabel = k;
                        minDst = dst;
                    }
                    cmps++;
                }
                if(this->labels[i] != clabel)
                {
                    UPDATED = true;
                }
                labels[i] = clabel;
                counts[clabel] += 1;
                loci = clabel*ndim;
                for(j = 0; j < pqm; j++)
                {
                    loc = j*pql;
                    for(k = 0; k < pql; k++)
                    {
                        Ds[loci + loc + k] += pqmap[j*pqn + data[i*pqm+j]][k];
                    }
                }
            }///for(i)

            for(i = 0; i < clnumb; i++)
            {
                if(counts[i] != 0)
                    for(j = 0; j < ndim; j++)
                    {
                        Cs[i*ndim + j] = Ds[i*ndim + j]/counts[i];
                    }
                else
                    for(j = 0; j < ndim; j++)
                    {
                        Cs[i*ndim+j] = 0;
                    }
            }

/// cs to pqcode
            for(unsigned cluster = 0; cluster < clnumb; cluster++)
            {
                loci = cluster*ndim;
                for(i = 0; i < pqm; i++)
                {
                    minDst = numeric_limits<double>::max();
                    loc = i*pql;
                    for(k = 0; k < pqn; k++)
                    {
                        dst = 0;
                        for(j = 0; j < pql; j++)
                        {
                            dst += (Cs[loci + loc + j] - pqmap[i*pqn + k][j])*(Cs[loci + loc + j] - pqmap[i*pqn + k][j]);
                        }
                        if(dst < minDst)
                        {
                            minDst = dst;
                            id = k;
                        }
                    }
                    PQC[cluster*pqm+i] = id;
                }
            }
            niter++;

            optEg = getI2(this->Ds, clust_num, ndim, counts);
            double tmperr = optEg - prvEg;
            prvEg = optEg;
            err = fabs(sumDst - sumDst0);
            sumDst0 = sumDst;
            cout<<niter<<"\t"<<optEg<<"\tdistortion\t"<<(en - optEg)/this->count<<endl;
        }
        while(UPDATED == true && niter < PQClust::M_ITER);
        if(sumDst0 < bstSumDst)
        {
            memcpy(arrayD, Ds, this->ndim*clust_num*sizeof(double));
            memcpy(bstLabels, labels, count*sizeof(int));
            bstSumDst = sumDst0;
            for(j = 0; j < clust_num; j++)
            {
                this->Ns[j] = counts[j];
            }
        }
    }///for(ri)



    memcpy(labels, bstLabels, count*sizeof(int));

    delete [] PQC;
    delete [] bstLabels;
    delete [] bstArrayD;
    delete [] counts;
    delete [] rseeds;
    delete [] Ds;
    Ds = NULL;
    bstLabels = NULL;
    bstArrayD = NULL;
    rseeds = NULL;
    counts = NULL;

    return clust_num;
}


bool PQClust::rndSeeds(const int k, int rseeds[], const int bound)
{
    unsigned int *ls = new unsigned[bound+1];
    unsigned int i = 0;
    for(i = 0; i <= bound; i++)
    {
        ls[i] = i;
    }
    random_shuffle(ls, ls+bound+1);
    memcpy(rseeds, ls, sizeof(unsigned)*(k+1));
    delete [] ls;
    return true;
}

/** random k seeds with the way of k-means++ **/
bool PQClust::kppSeeds(const int k, int rseeds[], const int bound)
{
    int i = 0, j = 0;
    double sum = 0.0f, rd = 0;
    int sel = 0;

    for(i = 0; i < k; i++)
    {
        rseeds[i] = i;
    }

    double *disbest = new double[bound];
    for(i = 0; i < bound; i++)
    {
        disbest[i] = numeric_limits<double>::max() + 0.0f;
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
            tmp = this->PQDist(data+i*pqm, data+j*pqm);
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

void PQClust::saveCenters(const char *dstfn, bool append)
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

int PQClust::fetchCenters(float *centers)
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

PQClust::~PQClust()
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

    if(this->tmpCs != NULL)
    {
        delete [] this->tmpCs;
        this->tmpCs = NULL;
    }

    if(this->Cs != NULL)
    {
        delete [] this->Cs;
        this->Cs = NULL;
    }
    if(this->pqdis != NULL)
    {
        delete [] this->pqdis;
        this->pqdis = NULL;
    }
    if(pqmap[0] != NULL)
    {
        delete [] pqmap[0];
    }
    pqmap.clear();
}

void PQClust::test()
{
    const char *srcMatFn1 = "/home/wlzhao/datasets/bignn/mat/sift1m/sift_base.txt";
    const char *srcMatFn2 = "/home/wlzhao/datasets/bignn/mat/sift_learn.txt";
    const char *srcMatFn3 = "itm_vlad_hesaff_flickr10m.itm";
    const char *dstfn = "/home/chenghaod/dataset/sift1m_pq/sift1m_pq8.out";
    const char *srcfn = "/home/chenghaod/dataset/sift1m_pq/sift1m_pq256seg8.txt";

    //int i = 1024;
    // do{
    ///PQClust *mykm = new PQClust();
    ///mykm->buildcluster(srcfn, dstfn, "rnd", "t1", "large", 1024, false);
    ///delete mykm;
    //    i = i*2;
    // }while(i <= 8192);
}

