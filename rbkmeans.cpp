#include "rbkmeans.h"

#include "ioagent.h"
#include "cleaner.h"
#include "vstring.h"
#include "pqmath.h"

#include <algorithm>
#include <iostream>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <cmath>

/**************************************************************************
@author: Wan-Lei Zhao
@email:  stonescx@gmail.com
@date:   23/Oct/2010
@Copyright: All rights are reserved by the author
**************************************************************************/

using namespace std;

const unsigned int RBKMeans::NTRAILS = 1;
const float RBKMeans::EPS    = 0.5f;
const float RBKMeans::Err0   = 1.0;

RBKMeans::RBKMeans()
{
    D1 = D2 = NULL;
    tmpD1   = tmpD2 = NULL;
    tmpC1   = tmpC2 = NULL;
    arrayD  = NULL;
    lens    = NULL;
    kmMtd   = _rbkmn_;
    strcpy(mthStr, "_rbk_");
    cout<<"Method ........................... RB K-Means\n";
}

bool RBKMeans::init(const char *srcfn)
{
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
        //this->sdata    = IOAgent::loadSparse(srcfn, this->count, this->ndim);
        //this->dataType = 1;
    }
    else
    {
        this->data = IOAgent::loadItms(srcfn, "fsvtab", this->count, this->ndim);
        //cout<<"Unrecognizable input file format!!!\n";
        //this->data = NULL;
        //exit(0);
    }

    cout<<this->count<<"x"<<this->ndim<<endl;

    if(this->data == NULL)
    {
        cout<<"Exceptions ....................... Loading matrix failed!\n";
        exit(0);
    }

    this->D1     = new float[this->ndim];
    this->D2     = new float[this->ndim];
    this->tmpD1  = new float[this->ndim];
    this->tmpD2  = new float[this->ndim];
    this->tmpC1  = new float[this->ndim];
    this->tmpC2  = new float[this->ndim];
    this->_REFER_= false;
    return true;
}

bool RBKMeans::init(float *mat, const int row, const int dim)
{
    assert(mat);
    refresh();
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
        this->data = NULL;
    }

    this->data   = mat;
    this->count  = row;
    this->ndim   = dim;
    this->D1     = new float[this->ndim];
    this->D2     = new float[this->ndim];
    this->tmpD1  = new float[this->ndim];
    this->tmpD2  = new float[this->ndim];
    this->tmpC1  = new float[this->ndim];
    this->tmpC2  = new float[this->ndim];

    this->_INIT_  = true;
    this->_REFER_ = true;
    return true;
}

bool RBKMeans::refresh()
{
    this->count  = 1;
    this->ndim   = 0;
    this->_INIT_ = false;

    if(this->D1 != NULL)
    {
        delete [] this->D1;
        this->D1 = NULL;
    }

    if(this->D2 != NULL)
    {
        delete [] this->D2;
        this->D2 = NULL;
    }

    if(this->tmpC1 != NULL)
    {
        delete [] this->tmpC1;
        this->tmpC1 = NULL;
    }

    if(this->tmpC2 != NULL)
    {
        delete [] this->tmpC2;
        this->tmpC2 = NULL;
    }

    if(this->tmpD1 != NULL)
    {
        delete [] this->tmpD1;
        this->tmpD1 = NULL;
    }

    if(this->tmpD2 != NULL)
    {
        delete [] this->tmpD2;
        this->tmpD2 = NULL;
    }

    return true;
}

bool RBKMeans::config(const char *_seed_, const char *lg_first, const char *crtrn, int verbose)
{
    if(!strcmp(lg_first, "best"))
    {
        this->LG_FIRST = false;
    }
    else if(!strcmp(lg_first, "large"))
    {
        this->LG_FIRST = true;
    }
    else
    {
        this->LG_FIRST = true;
    }

    if(verbose)
        cout<<"Distance function ................ cos\n";

    if(verbose)
        cout<<"Optimization function ............ ";
    if(!strcmp(crtrn, "i1"))
    {
        optim_func = &RBKMeans::incrOptz;
        //myoptz     = I1;
        if(verbose)
            cout<<"I1\n";
    }
    else if(!strcmp(crtrn, "i2"))
    {
        optim_func = &RBKMeans::incrOptz;
        if(verbose)
            cout<<"I2\n";
    }
    else
    {
        cout<<"Unkown optimize option '"<<crtrn<<"'!\n";
        exit(0);
    }

    this->seed = _rnd_;

    return true;
}

/***implementation of criterion I1***/
double RBKMeans::I1(const float *D1, const float *D2, const int dim,
                          const int n1, const int n2, double &E1, double &E2)
{
    double E = 0.0;
    E1 = 0;
    E2 = 0;
    if(n1 > 0)
        E1 = RBKMeans::l2_norm(D1, this->ndim)/n1;
    if(n2 > 0)
        E2 = RBKMeans::l2_norm(D2, this->ndim)/n2;
    E  = E1 + E2;
    return E;
}

/***implementation of criterion I2***/
double RBKMeans::I2(const float *D1, const float *D2, const int dim, double &E1, double &E2)
{
    double E = 0;
    E1 = RBKMeans::l2_norm(D1, this->ndim);
    E2 = RBKMeans::l2_norm(D2, this->ndim);
    E  = sqrt(E1) + sqrt(E2);
    return E;
}


double RBKMeans::getI2(const double *Ds, const int dim, const unsigned int clust_num)
{
    double sumE = 0, E = 0;
    for(unsigned int i = 0; i < clust_num; i++)
    {
        E = PQMath::dvec_norm(Ds+i*this->ndim, this->ndim, 2);
        E = E*E;
        sumE += E;
    }

    return sumE;
}

int RBKMeans::clust(const unsigned int clust_num, const char *dstfn, const int verbose)
{
    int clabel = 0, n1 = 0, n2 = 0;
    unsigned int i = 0, j = 0;
    double avgD1 = 0, avgD2 = 0, optEg = 0, maxEg = 0;
    bool OPTM = false;
    vector<NNItem*>::iterator vit;
    NNItem *crntItm = NULL;

    double *tmpDs  = new double[clust_num*this->ndim];
    int *tmpLabels = new int[this->count];
    int *tmpNs     = new int[clust_num];

    i = 1;
    for(j = 0; j < clust_num; j++)
    {
        this->infoMap[j].E = 0.0f;
        this->infoMap[j].n = 0;
    }
    cout<<"cluster number: "<<clust_num<<endl;
    if(verbose)
        cout<<"Clustering ....................... on progress\n";

    for(j = 0; j < NTRAILS; j++)
    {
        clabel = 0;
        memset(this->arrayD, 0, sizeof(double)*clust_num*this->ndim);
        memset(this->labels, 0, sizeof(int)*this->count);
        memset(this->Ns,     0, sizeof(int)*clust_num);
        mvs[0] = mvs[1] = 0;
        if(this->LG_FIRST)
        {
            crntItm = new NNItem(clabel, this->count);
        }
        else
        {
            crntItm = new NNItem(clabel, 0.0f);
        }
        i = 1;
        sorted_stack.push_back(crntItm);

        while(i < clust_num)
        {

            OPTM = incrOptz(clabel, avgD1, avgD2, n1, n2, i);

            if(OPTM)
            {
                crntItm           = sorted_stack[0];
                crntItm->index    = clabel;
                crntItm->size     = n1;
                crntItm->val      = avgD1;
                crntItm = new NNItem(i, avgD2, n2);
                crntItm->val      = avgD2;
                sorted_stack.push_back(crntItm);
                Ns[clabel]        = n1;
                Ns[i]             = n2;
                i++;   ///for next cluster
            }
            else
            {
                crntItm->dvd      = false;
                cout<<"unseparable!\n";
            }

            if(this->LG_FIRST)
            {
                stable_sort(sorted_stack.begin(), sorted_stack.end(), NNItem::LGSZcomparer);
            }
            else
            {
                stable_sort(sorted_stack.begin(), sorted_stack.end(), NNItem::LGVALcomparer);
            }

            for(vit = sorted_stack.begin(); vit != sorted_stack.end(); vit++)
            {
                crntItm = *vit;
                if(crntItm->size > 1 && crntItm->dvd)
                    break;
            }

            if(crntItm->size == 1 || !crntItm->dvd)
            {
                cout<<"Failed to reach expected cluster number!\n";
                break;
            }
            else
            {
                clabel = crntItm->index;
            }
        }///while(i)
        //cout<<"bug 1.14\n";
        optEg = getI2(this->arrayD, this->ndim, clust_num);

        if(optEg > maxEg)
        {
            memcpy(tmpDs,     this->arrayD, sizeof(double)*clust_num*this->ndim);
            memcpy(tmpLabels, this->labels, sizeof(int)*this->count);
            memcpy(tmpNs,     this->Ns,     sizeof(int)*clust_num);
            maxEg = optEg;
        }

        Cleaner::clearVector(sorted_stack);
    }

    memcpy(this->arrayD, tmpDs,     sizeof(double)*clust_num*this->ndim);
    memcpy(this->labels, tmpLabels, sizeof(int)*this->count);
    memcpy(this->Ns,     tmpNs,     sizeof(int)*clust_num);


    for(j = 0; j < clust_num; j++)
    {
        this->infoMap[j].n = this->Ns[j];
    }

    /**/
    int k = 0;
    memset(this->arrayD, 0, sizeof(double)*clust_num*this->ndim);

    for(i = 0; i < this->count; i++)
    {
        k = this->labels[i];
        if(k < 0)
        continue;

        for(j = 0; j < ndim; j++)
        {
            this->arrayD[k*ndim+j] += this->data[i*ndim+j]*this->lens[i];
            ///this->data[i*ndim+j] = this->data[i*ndim+j]*this->lens[i];
        }
    }
    /**/

    this->calAVGDist(this->arrayD, clust_num, this->infoMap);
    //this->refine(clust_num, 1);
    //this->calAVGDist(this->arrayD, clust_num, this->infoMap);

    if(strlen(dstfn) > 0)
    {
        stable_sort(sorted_stack.begin(), sorted_stack.end(), NNItem::LLIDXcomparer);
        char infofn[512];
        VString::parsePath(infofn, dstfn);
        strcat(infofn, "_clust_info.txt");
        ///RBKMeans::printvect(sorted_stack, infofn);
        RBKMeans::save_clust(dstfn);
    }
    Cleaner::clearVector(sorted_stack);
    return clust_num;
}

bool RBKMeans::incrOptz(const int clabel, double &E1, double &E2, int &n1, int &n2, const int  nwlbl)
{
    vector<int> crntClust;
    unsigned int i, j, loc1, loc2, loc, r;
    unsigned int cpysize = sizeof(float)*this->ndim;

    for(i = 0; i < count; i++)
    {
        if(clabel == this->labels[i])
        {
            crntClust.push_back(i);
        }
    }
    unsigned int numb = crntClust.size();
    if(numb == 2)
    {
        E2 = E1 = 1.0f;
        n1 = n2 = 1;
        i  = crntClust[1];

        loc1 = crntClust[0]*this->ndim;
        loc2 = crntClust[1]*this->ndim;
        for(j = 0; j < this->ndim; j++)
        {
            arrayD[clabel*this->ndim+j] = this->data[loc1+j];
            arrayD[nwlbl*this->ndim+j]  = this->data[loc2+j];
        }
        this->labels[i] = nwlbl;
        crntClust.clear();
        return true;
    }

    int *mylabels   = new int[numb];
    memset(mylabels,   0, sizeof(int)*numb);

    unsigned sed1 = 0, sed2 = 0, iter = 0;
    double tmpEg = 0, optEg = 0, tmpE1, tmpE2;
    int v = 0, tmpn1 = 0, tmpn2 = 0;
    float sim1 = 0,  sim2 = 0;

    bool UPDATE = true, _IMPRVED_ = false;
    n1 = n2 = 0;

    /**** Initial assigment ****************/
    this->randTwin(sed1, sed2, crntClust);
    ///this->KppTwin(sed1, sed2, crntClust);
    loc1 = sed1*this->ndim;
    memcpy(tmpC1, data+loc1, cpysize);
    loc2 = sed2*this->ndim;
    memcpy(tmpC2, data+loc2, cpysize);
    memset(D1, 0, cpysize);
    memset(D2, 0, cpysize);
    tmpn1 = tmpn2 = 0;
    optEg = 0;

    for(i = 0; i < numb; i++)
    {
        v    = crntClust[i];
        sim1 = RBKMeans::cos(tmpC1, 0, data, v, this->ndim);
        sim2 = RBKMeans::cos(tmpC2, 0, data, v, this->ndim);
        loc  = v*this->ndim;

        if(sim1 > sim2)
        {
            for(j = 0; j < this->ndim; j++)
            {
                this->D1[j] += data[loc+j];
            }
            tmpn1++;
            mylabels[i] = -1;
        }
        else
        {
            for(j = 0; j < this->ndim; j++)
            {
                D2[j] += data[loc+j];
            }
            tmpn2++;
            mylabels[i] = v;
        }
    }

    optEg = I1(D1, D2, this->ndim, tmpn1, tmpn2, tmpE1, tmpE2);
    tmpEg = 0;

    /******** Incremental optimization *****/
    iter  = 0;
    int *rl = new int[numb];

    do
    {
        _IMPRVED_ = false;
        for(i = 0; i < numb; i++)
        {
            rl[i] = i;
        }
        random_shuffle(rl, rl+numb);
        for(iter = 0; iter < numb; iter++)
        {
            memcpy(tmpD1, D1, sizeof(float)*this->ndim);
            memcpy(tmpD2, D2, sizeof(float)*this->ndim);
            r   = rl[iter];
            v   = crntClust[r];
            loc = v*this->ndim;
            if(mylabels[r] == -1)
            {
                if(tmpn1 > 1)
                {
                    for(j = 0; j < this->ndim; j++)
                    {
                        tmpD1[j] -= data[loc+j];
                        tmpD2[j] += data[loc+j];
                    }
                    tmpn1--;
                    tmpn2++;
                    mylabels[r] = v;
                    UPDATE = true;
                }
                else
                {
                    UPDATE = false;
                }
            }
            else
            {
                if(tmpn2 > 1)
                {
                    for(j = 0; j < this->ndim; j++)
                    {
                        tmpD1[j] += data[loc+j];
                        tmpD2[j] -= data[loc+j];
                    }
                    tmpn1++;
                    tmpn2--;
                    mylabels[r] = -1;
                    UPDATE = true;
                }
                else
                {
                    UPDATE = false;
                }
            }

            tmpEg = I1(tmpD1, tmpD2,  this->ndim,  tmpn1, tmpn2, tmpE1, tmpE2);

            if(tmpEg > optEg)
            {
                memcpy(D1, tmpD1, cpysize);
                memcpy(D2, tmpD2, cpysize);
                E1 = tmpE1;
                E2 = tmpE2;
                n1 = tmpn1;
                n2 = tmpn2;
                optEg = tmpEg + 0.001;
                _IMPRVED_ = true;
            }
            else
            {
                if(mylabels[r] == -1)
                {
                    ///roll-back
                    if(UPDATE)
                    {
                        tmpn1--;
                        tmpn2++;
                        mylabels[r] = v;
                    }
                }
                else
                {
                    ///roll-back
                    if(UPDATE)
                    {
                        tmpn1++;
                        tmpn2--;
                        mylabels[r] = -1;
                    }
                }
            }
        }//for(iter)
    }
    while(_IMPRVED_);

    delete [] rl;
    rl = NULL;
    crntClust.clear();

    if(tmpn1 == 0 || tmpn2 == 0)
    {
        delete [] mylabels;
        mylabels = NULL;
        return false;
    }

    optEg = this->I1(D1, D2, this->ndim,  tmpn1, tmpn2, E1, E2);

    n1 = tmpn1;
    n2 = tmpn2;
    E1 = (E1 - n1)/2;
    E2 = (E2 - n2)/2;

    loc1 = clabel*this->ndim;
    loc2 = nwlbl*this->ndim;
    cpyCenters(arrayD + loc1, D1, this->ndim);
    cpyCenters(arrayD + loc2, D2, this->ndim);

    for(i = 0; i < numb; i++)
    {
        v = mylabels[i];
        if(v != -1)
        {
            labels[v] = nwlbl;
        }
    }

    delete [] mylabels;
    mylabels = NULL;

    return true;
}

double RBKMeans::l2_norm(const float *vect1, const unsigned int dim)
{
    double  sz = 0;
    for(unsigned int i = 0; i < dim; i++)
    {
        sz += vect1[i]*vect1[i];
    }
    return sz;
}

bool RBKMeans::randTwin(unsigned int &r1, unsigned int &r2, vector<int> &vect)
{
    const unsigned int bound   = vect.size();
    unsigned int t = 0;

    r1 = (unsigned int)floor((rand()/(RAND_MAX+1.0f))*bound);
    r1 = (r1 >= bound)?(bound-1):r1;
    r2 = (unsigned int)floor((rand()/(RAND_MAX+1.0f))*bound);
    r2 = (r2 >= bound)?(bound-1):r2;

    while(r1 == r2 && t < 500)
    {
        r2 = (unsigned int)floor((rand()/(RAND_MAX+1.0f))*bound);
        r2 = (r2 >= bound)?(bound-1):r2;
        t++;
    }
    if(r1 == r2)
    {
        r1 = 0;
        r2 = 2;
    }
    r1 = vect[r1];
    r2 = vect[r2];
    return true;
}

bool RBKMeans::randTwo(const int bound, unsigned int &r1, unsigned int &r2)
{
    unsigned int sed0 = (unsigned int)time(NULL);
    r1 = (unsigned int)floor((rand_r(&sed0)/(RAND_MAX+1.0f))*bound);
    r2 = (unsigned int)floor((rand_r(&sed0)/(RAND_MAX+1.0f))*bound);

    while(r1 == r2)
    {
        r2 = (unsigned int)floor((rand_r(&sed0)/(RAND_MAX+1.0f))*bound);
    }
    assert(r1 != r2);
    return true;
}

unsigned int RBKMeans::randOne(const int bound)
{
    unsigned int sed0 = time(NULL);
    unsigned int r1   = (unsigned int)floor((rand_r(&sed0)/(RAND_MAX+1.0f))*bound);
    return r1;
}

/** random k seeds with the way of k-means++ **/
bool RBKMeans::KppTwin(unsigned int &r1, unsigned int &r2, vector<int> &vect)
{
    long i = 0, j = 0;
    double sum = 0.0f;
    const int bound = vect.size();
    long c1 = 0, c2 = 0;

    if(this->_INIT_ == false)
    {
        cout<<"Error ........................ initialize K-Means++ first!\n";
        exit(0);
    }
    int rseeds[2], sel;

    for(i = 0; i < 2; i++)
    {
        rseeds[i] = i;
    }

    float *disbest = new float[bound];
    for(i = 0; i < bound; i++)
    {
        disbest[i] = RAND_MAX;
    }
    float * distmp    = new float[bound];
    unsigned int sed0 = (unsigned int)time(NULL);
    float tmp = 0;

    rseeds[0] = rand_r(&sed0) % bound;
    for (i = 1 ; i < 2; i++)
    {
        ///cout<<"step "<<i<<endl;
        sel = rseeds[i - 1];
        c1  = vect[sel];
        for (j = 0 ; j < bound; j++)
        {
            c2  = vect[j];
            tmp = this->l2(this->data, c2, this->data, c1, this->ndim);
            if(tmp < disbest[j])
            {
                disbest[j] = tmp;
            }
        }

        /** convert the best distances to probabilities **/
        memcpy (distmp, disbest, bound * sizeof (*distmp));

        sum = fvec_norm(distmp, bound, 1);
        fvec_scale(distmp, bound, 1.0/sum);
        double rd = rand_r(&sed0)/(RAND_MAX+1.0f);

        for (j = 0 ; j < bound-1; j++)
        {
            rd -= distmp[j];
            if (rd < 0)
                break;
        }

        rseeds[i] = j;
    }
    r1 = rseeds[0];
    r2 = rseeds[1];

    r1 = vect[r1];
    r2 = vect[r2];

    if(r1 == r2)
    {
        for (j = 0 ; j < bound; j++)
        {
            if(vect[j] != r1)
            {
                break;
            }
        }
        r2 = vect[j];
    }

    delete [] disbest;
    delete [] distmp;
    disbest = distmp = NULL;
    return true;
}

bool RBKMeans::cpyCenters(double *dstArray, float *srcArray, const unsigned int ndim)
{
    unsigned int i = 0;
    for(i = 0; i < ndim; i++)
    {
        dstArray[i] = srcArray[i];
    }
    return true;
}

void RBKMeans::saveCenters(const char *dstfn, bool append)
{
    unsigned int i = 0, j = 0, loc = 0, cloc = 0;
    unsigned int clabel = 0, rCNum = 0;
    float val = 0;

    if(!this->_INIT_||this->clnumb == 0)
    {
        return ;
    }

    for(i = 0; i < this->clnumb; i++)
    {
        if(this->infoMap[clabel].n > 0)
        {
            rCNum++;
        }
    }

    ofstream *outStrm = NULL;
    if(append)
    {
        outStrm = new ofstream(dstfn, ios::app);
    }
    else
    {
        outStrm = new ofstream(dstfn, ios::out);
        (*outStrm)<<rCNum<<" "<<this->ndim<<endl;;
    }

    float *centers = new float[this->clnumb*this->ndim];
    memset(centers, 0, this->clnumb*this->ndim*sizeof(float));
    for(i = 0; i < this->count; i++)
    {
        clabel = this->labels[i];
        if(clabel >= 0 && this->infoMap[clabel].n > 0)
        {
            loc  = i*this->ndim;
            cloc = clabel*this->ndim;
            val  = this->lens[i];
            for(j = 0; j < this->ndim; j++)
            {
                centers[cloc+j] += (this->data[loc+j]*val)/this->infoMap[clabel].n;
            }
        }
    }

    for(clabel = 0; clabel < this->clnumb; clabel++)
    {
        cloc   = clabel*this->ndim;
        if(this->infoMap[clabel].n > 0)
        {
            for(j = 0; j < this->ndim; j++)
            {
                (*outStrm)<<centers[cloc+j]<<" ";
            }
            (*outStrm)<<endl;
        }
        else
        {
            cout<<clabel<<"\t cluster is empty!\n";
        }
    }
    outStrm->close();
    delete [] centers;
    centers = NULL;
    return ;
}

int RBKMeans::fetchCenters(float *centers)
{
    unsigned int i, j, loc = 0, cloc = 0;
    unsigned int clabel = 0;
    float val  = 0;
    assert(centers);
    memset(centers, 0, this->clnumb*this->ndim*sizeof(float));

    if(!this->_INIT_||this->clnumb == 0)
        return 0;

    for(clabel = 0; clabel < this->clnumb; clabel++)
    {
        if(this->infoMap[clabel].n == 0)
        {
            cout<<clabel<<"\t cluster is empty!\n";
        }
    }

    for(i = 0; i < this->count; i++)
    {
        clabel = this->labels[i];
        if(clabel >= 0 && this->infoMap[clabel].n > 0)
        {
            loc  = i*this->ndim;
            cloc = clabel*this->ndim;
            val  = this->lens[i];
            for(j = 0; j < this->ndim; j++)
            {
                centers[cloc+j] += (float)(val*this->data[loc+j])/this->infoMap[clabel].n;
            }
        }
    }

    return this->clnumb;
}

float RBKMeans::cos(const float v1[], const unsigned int s1, const float v2[],
                    const unsigned int s2, const unsigned int d0)
{
    float w1  = 0.0f, w2 = 0.0f;
    float val = 0.0f;
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

float RBKMeans::l2(const float v1[], const unsigned int s1, const float v2[], const unsigned int s2, const unsigned int d0)
{
    float w1 = 0.0f, w2 = 0.0f;
    unsigned int loc1 = d0 * s1;
    unsigned int loc2 = d0 * s2;

    for(unsigned int i = 0; i < d0; i++)
    {
        w1  = v1[loc1+i] - v2[loc2+i];
        w2 += w1*w1;
    }

    w2 = sqrt(w2);

    return w2;
}

double RBKMeans::fvec_norm(const float *v, const long n, const double norm)
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

bool RBKMeans::fvec_scale(float *v, const long n, const double sc)
{
    assert(sc!=0.0f);

    for(long i = 0; i < n; i++)
    {
        v[i] = v[i]*sc;
    }
    return true;
}

RBKMeans::~RBKMeans()
{
    if(this->D1 != NULL)
    {
        delete [] this->D1;
        this->D1 = NULL;
    }

    if(this->D2 != NULL)
    {
        delete [] this->D2;
        this->D2 = NULL;
    }
    if(this->tmpD1 != NULL)
    {
        delete [] this->tmpD1;
        this->tmpD1 = NULL;
    }
    if(this->tmpD2 != NULL)
    {
        delete [] this->tmpD2;
        this->tmpD2 = NULL;
    }
    if(this->tmpC1 != NULL)
    {
        delete [] this->tmpC1;
        this->tmpC1 = NULL;
    }
    if(this->tmpC2 != NULL)
    {
        delete [] this->tmpC2;
        this->tmpC2 = NULL;
    }

}

void RBKMeans::test()
{
    /**/
    const char *srcMatFn1 = "/home/wlzhao/datasets/bignn/sift1m/sift_base.txt";
    const char *srcMatFn2 = "/home/wlzhao/datasets/bignn/sift1m/sift_learn.txt";
    const char *srcMatFn3 = "itm_vlad_hesaff_flickr10m.";
    const char *dstfn  = "/home/wlzhao/datasets/bignn/vocab/rslt.txt";
    const char *dstfn1 = "/home/wlzhao/datasets/bignn/rbk_clustxx.txt";

    RBKMeans *mykmean = new RBKMeans();
    mykmean->buildcluster(srcMatFn2, dstfn1, "non", "large", "i2", 1024, false);
    delete mykmean;
}
