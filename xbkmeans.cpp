#include "xbkmeans.h"

#include "print2scrn.h"
#include "randseed.h"
#include "vstring.h"
#include "cleaner.h"
#include "ioagent.h"
#include "pqmath.h"
#include "timer.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>
#include <cmath>

const int   XBKMeans::NTRAILS = 1;
const float XBKMeans::Err0    = 0.01f;
const int   XBKMeans::NIter0  = 100;
const int   XBKMeans::NRfn    = 10;
const unsigned interval = 100;
using namespace std;

XBKMeans::XBKMeans()
{
    D1 = D2 = NULL;
    C1 = C2 = NULL;
    tmpD1   = tmpD2 = NULL;
    tmpC1   = tmpC2 = NULL;
    innDcts = NULL;
    kmMtd   = _xbkmn_;
    strcpy(mthStr, "_xbk_");
    this->_REFER_ = false;
    cout<<"Method ........................... XB K-Means\n";
}

bool  XBKMeans::init(const char *srcfn)
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
          this->dataType = 0;
    }
    else if(VString::endWith(srcfn, ".fvecs"))
    {
         this->data = IOAgent::load_fvecs(srcfn, this->ndim, this->count);
         this->dataType = 0;
    }
    else if(VString::endWith(srcfn, ".mat"))
    {
         this->data = IOAgent::loadDat(srcfn, this->count, this->ndim);
         /** this->sdata    = IOAgent::loadSparse(srcfn, this->count, this->ndim);
         this->dataType = 1;/**/
    }
    else
    {
         this->data = IOAgent::loadItms(srcfn, "fsvtab", this->count, this->ndim);
         //cout<<"Unrecognizable input file format!!!\n";
         //this->data = NULL;
         //exit(0);
    }

    cout<<this->count<<"x"<<this->ndim<<endl;
    if(dataType)
    {
        if(this->sdata.data == NULL)
        {
            cout<<"Exceptions ....................... Loading matrix failed!\n";
            exit(0);
        }
    }
    else if(this->data == NULL)
    {
        cout<<"Exceptions ....................... Loading matrix failed!\n";
        exit(0);
    }

    this->D1     = new double[this->ndim];
    this->D2     = new double[this->ndim];
    this->C1     = new double[this->ndim];
    this->C2     = new double[this->ndim];
    this->tmpD1  = new double[this->ndim];
    this->tmpD2  = new double[this->ndim];
    this->tmpC1  = new double[this->ndim];
    this->tmpC2  = new double[this->ndim];

    this->_REFER_= false;

    return true;
}

bool  XBKMeans::init(float *mat, const int row, const int dim)
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

    this->D1     = new double[this->ndim];
    this->D2     = new double[this->ndim];
    this->C1     = new double[this->ndim];
    this->C2     = new double[this->ndim];
    this->tmpD1  = new double[this->ndim];
    this->tmpD2  = new double[this->ndim];
    this->tmpC1  = new double[this->ndim];
    this->tmpC2  = new double[this->ndim];

    this->_INIT_  = true;
    this->_REFER_ = true;
    return true;
}

bool XBKMeans::refresh()
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

    if(this->C1 != NULL)
    {
        delete [] this->C1;
        this->C1 = NULL;
    }

    if(this->C2 != NULL)
    {
        delete [] this->C2;
        this->C2 = NULL;
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

XBKMeans::~XBKMeans()
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

    if(this->C1 != NULL)
    {
        delete [] this->C1;
        this->C1 = NULL;
    }

    if(this->C2 != NULL)
    {
        delete [] this->C2;
        this->C2 = NULL;
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

bool  XBKMeans::config(const char *_seed_, const char *lg_first, const char *crtrn, int verbose)
{
    ///if(verbose)
    cout<<"Bisecting stradegy ............... ";
    if(!strcmp(lg_first, "best"))
    {
        this->LG_FIRST = false;
        cout<<"best\n";
    }
    else if(!strcmp(lg_first, "large"))
    {
        this->LG_FIRST = true;
        cout<<"large\n";
    }
    else
    {
        this->LG_FIRST = true;
        cout<<"large\n";
    }

    if(verbose)
        cout<<"Seeding approach ................. ";

    if(!strcmp(_seed_, "rnd"))
    {
        if(verbose)
            cout<<"rand\n";
        sedFunc = &XBKMeans::rndTwin;
        seed    = _rnd_;
    }
    else if(!strcmp(_seed_, "kpp"))
    {
        if(verbose)
            cout<<"kpp\n";
        sedFunc = &XBKMeans::kppTwin;
        seed    = _kpp_;
    }
    else if(!strcmp(_seed_, "non"))
    {
        if(verbose)
            cout<<"non\n";
        sedFunc = NULL;
        seed    = _non_;
    }
    else
    {
        if(verbose)
            cout<<"kpp\n";
        sedFunc = &XBKMeans::kppTwin;
        seed    = _kpp_;
    }

    if(verbose)
        cout<<"Distance function ................ l2\n";


    if(verbose)
    cout<<"Optimization function ............ ";
    if(!strcmp(crtrn, "i1"))
    {
        optFunc = &XBKMeans::I1Func;
        myoptz     = _I1_;
        if(verbose)
        cout<<"I1\n";
    }
    else if(!strcmp(crtrn, "i2"))
    {
        optFunc = &XBKMeans::I2Func;
        myoptz  = _I2_;
        if(verbose)
        cout<<"I2\n";
    }
    else if(!strcmp(crtrn, "i3"))
    {
        optFunc = &XBKMeans::I3Func;
        myoptz  = _I3_;
        if(verbose)
        cout<<"I3\n";
    }
    else if(!strcmp(crtrn, "i4") )
    {
        ///optFunc = &XBKMeans::E1Func;
        myoptz  = _I4_;

        if(verbose)
        cout<<"e1\n";
    }
    else if(!strcmp(crtrn, "e1") )
    {
        optFunc = &XBKMeans::E1Func;
        myoptz  = _E1_;

        if(verbose)
        cout<<"e1\n";
    }
    else if(!strcmp(crtrn, "e2") )
    {
        optFunc = &XBKMeans::E2Func;
        myoptz  = _E2_;

        if(verbose)
        cout<<"e2\n";
    }
    else if(!strcmp(crtrn, "t1"))
    {
        ///optFunc = &XBKMeans::tkmOptz;
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

int XBKMeans::clust(const unsigned int clust_num, const char *dstfn, const int verbose)
{
    int clabel = 0, n1 = 0, n2 = 0;
    unsigned int i = 0, j = 0;
    double avgD1 = 0, avgD2 = 0, optEg = 0, maxEg = -RAND_MAX;
    double *tmpDs  = new double[clust_num*this->ndim];
    innDcts= new double[2];
    int *tmpLabels = new int[this->count];
    int *tmpNs     = new int[clust_num];
    vector<NNItem*>::iterator vit;
    NNItem *crntItm = NULL;
    bool OPTM = false;

    for(j = 0; j < clust_num; j++)
    {
        this->infoMap[j].E = 0.0f;
        this->infoMap[j].n = 0;
    }

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
        ///Timer *mytm = new Timer();
        ///mytm->start();
        while(i < clust_num)
        {
            if(dataType)
            {
                switch(myoptz)
                {
                case _I1_:
                    OPTM = incrOptzI1s(clabel, avgD1, avgD2, n1, n2, i);
                    break;
                case _I2_:
                {
                    OPTM = incrOptzI2s(clabel, avgD1, avgD2, n1, n2, i);
                    break;
                }
                case _E1_:
                    OPTM = incrOptzE1s(clabel, avgD1, avgD2, n1, n2, i);
                    break;
                case _E2_:
                    OPTM = incrOptzE2s(clabel, avgD1, avgD2, n1, n2, i);
                    break;
                case _I4_:
                    OPTM = incrOptzI4s(clabel, avgD1, avgD2, n1, n2, i);
                    break;
                }

            }
            else
                switch(myoptz)
                {
                case _I1_:
                    OPTM = incrOptzI1(clabel, avgD1, avgD2, n1, n2, i);
                    break;
                case _I2_:
                {
                    OPTM = incrOptzI2(clabel, avgD1, avgD2, n1, n2, i);
                    break;
                }
                case _I4_:
                {
                    OPTM = incrOptzI4(clabel, avgD1, avgD2, n1, n2, i);
                    break;
                }
                case _E1_:
                    OPTM = incrOptzE1(clabel, avgD1, avgD2, n1, n2, i);
                    break;
                case _E2_:
                    OPTM = incrOptzE2(clabel, avgD1, avgD2, n1, n2, i);
                    break;
                case _T1_:
                {
                    OPTM = tkmOptz(clabel, avgD1, avgD2, n1, n2, i);
                    break;
                }
                }
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
            // cout<<i<<"\t"<<clabel<<"\t"<<crntItm->val<<"\t"<<crntItm->size<<endl;
        }///while(i)
        ///mytm->end(true);
        switch(myoptz)
        {
        case _I1_:
            optEg = getI1(this->arrayD, clust_num, this->ndim, this->Ns);
            break;
        case _I2_:
            optEg = getI2(this->arrayD, clust_num, this->ndim, this->Ns);
            break;
        case _I3_:
            break;
        case _E1_:
            optEg = getE1(this->arrayD, clust_num, this->ndim, this->Ns);
            break;
        case _E2_:
        {
            optEg = getE2(this->arrayD, clust_num, this->ndim, this->Ns);
            break;
        }
        case _T1_:
        {
            optEg = getI2(this->arrayD, clust_num, this->ndim, this->Ns);
            break;
        }
        case _I4_:
        {
            optEg = getI2(this->arrayD, clust_num, this->ndim, this->Ns);
            break;
        }
        }

        if(false)
        {
            if(dataType)
                optEg = refineSparse(clust_num, XBKMeans::NRfn);
            else
                optEg = refine(clust_num, XBKMeans::NRfn);
        }
        if(optEg > maxEg)
        {
            memcpy(tmpDs,     this->arrayD, sizeof(double)*clust_num*this->ndim);
            memcpy(tmpLabels, this->labels, sizeof(int)*this->count);
            memcpy(tmpNs,     this->Ns,     sizeof(int)*clust_num);
            maxEg = optEg;
        }
        Cleaner::clearVector(sorted_stack);
        ///cout<<"Moves: "<<mvs[0]<<"\t"<<mvs[1]<<endl;
    }

    memcpy(this->arrayD, tmpDs,     sizeof(double)*clust_num*this->ndim);
    memcpy(this->labels, tmpLabels, sizeof(int)*this->count);
    memcpy(this->Ns,     tmpNs,     sizeof(int)*clust_num);

    for(i = 0 ; i < clust_num; i ++)
    {
        infoMap[i].n = Ns[i];
    }

    ///cout<<"Moves: "<<mvs[0]<<"\t"<<mvs[1]<<endl;
    double sumEg = 0;
    sumEg = this->calAVGDist(this->arrayD, clust_num, infoMap);
    this->r2refine(clust_num, clust_num/10);
    sumEg = this->calAVGDist(this->arrayD, clust_num, infoMap);

    if(strlen(dstfn) > 0)
    {
        char infofn[512];
        VString::parsePath(infofn, dstfn);
        strcat(infofn, "_clust_info.txt");
        //Print2Scrn::printvect(sorted_stack, infofn);
        XBKMeans::save_clust(dstfn);
    }
    delete [] tmpDs;
    delete [] tmpNs;
    delete [] tmpLabels;
    delete [] innDcts;
    tmpDs = NULL;
    tmpNs = NULL;
    tmpLabels = NULL;
    innDcts = NULL;
    return clust_num;
}

bool XBKMeans::incrOptzI2(const int clabel, double &E1, double &E2, int &n1, int &n2, const int nwlbl)
{
    unsigned int _loc1 = 0, _loc2 = 0, loc = 0;
    unsigned int i = 0, j = 0, r = 0;
    unsigned int cpysize = sizeof(double)*this->ndim;
    double tmp1 = 0, tmp2 = 0, len = 0;

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
        E2 = E1 = 0.0f;
        n1 = n2 = 1;
        i  = crntClust[1];
        this->labels[i] = nwlbl;

        _loc1 = crntClust[0]*this->ndim;
        _loc2 = crntClust[1]*this->ndim;
        for(j = 0; j < this->ndim; j++)
        {
            arrayD[clabel*this->ndim+j] = this->data[_loc1+j];
            arrayD[nwlbl*this->ndim+j]  = this->data[_loc2+j];
        }
        crntClust.clear();

        return true;
    }

    int *nwlabels = new int[numb];
    memset(nwlabels, 0, sizeof(int)*numb);
    double *tmpData = new double[ndim];

    int v = 0, tmpn1 = 0, tmpn2 = 0;
    unsigned int sed1 = 0, sed2 = 0, iter = 0;
    double tmpEg = 0, optEg = 0, dst1 = 0,  dst2 = 0;
    bool _UPDATE_ = true, F2S = false;
    n1 = n2 = 0;
    /**** Initial assigment ****************/
    if(this->seed != _non_)
    {
        (this->*sedFunc)(sed1, sed2, crntClust);
        _loc1 = sed1*this->ndim;
        _loc2 = sed2*this->ndim;
        for(j = 0; j < this->ndim; j++)
        {
            tmpC1[j] = this->data[_loc1+j];
            tmpC2[j] = this->data[_loc2+j];
        }
    }

    memset(D1, 0, cpysize);
    memset(D2, 0, cpysize);
    innDcts[0] = innDcts[1] = 0;

    for(i = 0; i < numb; i++)
    {
        v = crntClust[i];
        if(this->seed != _non_)
        {
            mvs[0]++;
            dst1 = PQMath::l2d(tmpC1, 0, this->data, v, this->ndim);
            dst2 = PQMath::l2d(tmpC2, 0, this->data, v, this->ndim);
        }
        else
        {
            dst1 = rand();
            dst2 = rand();
        }
        loc  = v*this->ndim;

        if(dst1 < dst2)
        {
            for(j = 0; j < this->ndim; j++)
            {
                this->D1[j] += this->data[loc+j];
            }
            tmpn1++;
            nwlabels[i] = -1;
        }
        else
        {
            for(j = 0; j < this->ndim; j++)
            {
                D2[j] += this->data[loc+j];
            }
            tmpn2++;
            nwlabels[i] = v;
        }
    }

    if(tmpn1 == 0)
    {
        optEg = 0;
        n1 = tmpn2;
        n2 = tmpn1;
        E1 = 0;
        E2 = 0;
        _loc1 = nwlbl*this->ndim;
        memset(this->arrayD+_loc1, 0, cpysize);
        delete [] nwlabels;
        nwlabels = NULL;

        crntClust.clear();
        return false;
    }

    if(tmpn2 == 0)
    {
        optEg = 0;
        n1 = tmpn1;
        n2 = tmpn2;
        E2 = 0;
        E1 = 0;
        _loc1 = nwlbl*this->ndim;
        memset(this->arrayD+_loc1, 0, cpysize);
        delete [] nwlabels;
        nwlabels = NULL;

        crntClust.clear();
        return false;
    }
    innDcts[0] = PQMath::dvec_norm(D1, this->ndim, 2);
    innDcts[1] = PQMath::dvec_norm(D2, this->ndim, 2);
    innDcts[0] = innDcts[0]*innDcts[0];
    innDcts[1] = innDcts[1]*innDcts[1];

    optEg = innDcts[0]/tmpn1 +  innDcts[1]/tmpn2;

    int *rl = new int[numb];

    /******** Incremental optimization *****/
    do
    {
        for(i = 0; i < numb; i++)
        {
            rl[i] = i;
        }
        random_shuffle(rl, rl + numb);
        _UPDATE_ = false;
        for(iter = 0; iter < numb; iter++)
        {
            mvs[0]++;
            memcpy(tmpD1, D1, cpysize);
            memcpy(tmpD2, D2, cpysize);
            r   = rl[iter];
            v   = crntClust[r];
            loc = v*this->ndim;
            for(j = 0; j < this->ndim; j++)
            {
                tmpData[j] = data[loc+j];
            }
            F2S = false;

            if(nwlabels[r] == -1)
            {
                if(tmpn1 > 1)
                {
                    tmpn1--;
                    tmpn2++;
                    nwlabels[r] = v;
                    F2S = true;
                }
                else
                {
                    continue;
                }
            }
            else
            {
                if(tmpn2 > 1)
                {

                    tmpn1++;
                    tmpn2--;
                    nwlabels[r] = -1;
                    F2S = false;
                }
                else
                {
                    continue;
                }
            }

            tmp1 = tmp2 = 0;
            for(j = 0; j < this->ndim; j++)
            {
                tmp1 += tmpD1[j]*tmpData[j];
                tmp2 += tmpD2[j]*tmpData[j];
            }

            len = this->lens[v];
            if(F2S)
            {
                tmp1 = innDcts[0] - 2*tmp1 + len;
                tmp2 = innDcts[1] + 2*tmp2 + len;

            }
            else
            {
                tmp1 = innDcts[0] + 2*tmp1 + len;
                tmp2 = innDcts[1] - 2*tmp2 + len;
            }

            tmpEg = tmp1/tmpn1 + tmp2/tmpn2;

            if(tmpEg > (optEg+0.00001))
            {
                if(F2S)
                {
                    for(j = 0; j < this->ndim; j++)
                    {
                        tmpD1[j] -= tmpData[j];
                        tmpD2[j] += tmpData[j];
                    }
                }
                else
                {
                    for(j = 0; j < this->ndim; j++)
                    {
                        tmpD1[j] += tmpData[j];
                        tmpD2[j] -= tmpData[j];
                    }
                }
                memcpy(D1, tmpD1, cpysize);
                memcpy(D2, tmpD2, cpysize);
                optEg = tmpEg;
                innDcts[0] = tmp1;
                innDcts[1] = tmp2;
                _UPDATE_ = true;
                mvs[1]++;
            }
            else
            {
                ///roll-back if being updated
                if(nwlabels[r] == -1)
                {
                    tmpn1--;
                    tmpn2++;
                    nwlabels[r] = v;
                }
                else
                {
                    tmpn1++;
                    tmpn2--;
                    nwlabels[r] = -1;
                }
            }
        }
    }
    while(_UPDATE_);

    _loc1 = clabel*this->ndim;
    memcpy(this->arrayD+_loc1, D1, cpysize);
    _loc1 = nwlbl*this->ndim;
    memcpy(this->arrayD+_loc1, D2, cpysize);

    if(optEg > 0)
    {
        n1 = tmpn1;
        n2 = tmpn2;
    }

    for(i = 0; i < numb; i++)
    {
        v = nwlabels[i];
        if(v != -1)
        {
            labels[v] = nwlbl;
        }
    }

    double len1, len2, tmpE1 = 0, tmpE2 = 0;
    len1 = len2 = 0;
    int k = 0;
    for(i = 0; i < numb; i++)
    {
        v = nwlabels[i];
        k = crntClust[i];
        if(v != -1)
        {
            len2 += lens[k];
        }
        else
        {
            len1 += lens[k];
        }
    }

    /**estimate average intra-cluster distance*/
    for(j = 0; j < this->ndim; j++)
    {
        C1[j] = D1[j]/n1;
        C2[j] = D2[j]/n2;
    }
    tmpE1 = PQMath::dvec_norm(C1, ndim, 2);
    tmpE2 = PQMath::dvec_norm(C2, ndim, 2);
    if(n1 > 1)
    {
        E1 = len1 - n1*tmpE1*tmpE1;
        E1 = 2*(E1/(n1-1));
    }
    else
    {
        E1 = 0;
    }

    if(n2 > 1)
    {
        E2 = len2 - n2*tmpE2*tmpE2;
        E2 = 2*(E2/(n2-1));
    }
    else
    {
        E2 = 0;
    }

    crntClust.clear();
    delete [] nwlabels;
    delete [] tmpData;
    delete [] rl;
    nwlabels = NULL;
    tmpData = NULL;
    return true;
}

bool XBKMeans::incrOptzI4(const int clabel, double &E1, double &E2, int &n1, int &n2, const int nwlbl)
{
    unsigned int _loc1 = 0, _loc2 = 0, loc = 0;
    unsigned int i = 0, j = 0, r = 0;
    unsigned int cpysize = sizeof(double)*this->ndim;
    double tmp1 = 0, tmp2 = 0, len = 0;

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
        E2 = E1 = 0.0f;
        n1 = n2 = 1;
        i  = crntClust[1];
        this->labels[i] = nwlbl;

        _loc1 = crntClust[0]*this->ndim;
        _loc2 = crntClust[1]*this->ndim;
        for(j = 0; j < this->ndim; j++)
        {
            arrayD[clabel*this->ndim+j] = this->data[_loc1+j];
            arrayD[nwlbl*this->ndim+j]  = this->data[_loc2+j];
        }
        crntClust.clear();

        return true;
    }

    int *nwlabels = new int[numb];
    memset(nwlabels, 0, sizeof(int)*numb);

    int v = 0, tmpn1 = 0, tmpn2 = 0;
    unsigned int sed1 = 0, sed2 = 0, iter = 0;
    double tmpEg = 0, optEg = 0, dst1 = 0,  dst2 = 0;
    bool _UPDATE_ = true, F2S = false;
    n1 = n2 = 0;
    /**** Initial assigment ****************/
    if(this->seed != _non_)
    {
        (this->*sedFunc)(sed1, sed2, crntClust);
        _loc1 = sed1*this->ndim;
        _loc2 = sed2*this->ndim;
        for(j = 0; j < this->ndim; j++)
        {
            tmpC1[j] = this->data[_loc1+j];
            tmpC2[j] = this->data[_loc2+j];
        }
    }

    memset(D1, 0, cpysize);
    memset(D2, 0, cpysize);
    innDcts[0] = innDcts[1] = 0;

    for(i = 0; i < numb; i++)
    {
        v = crntClust[i];
        if(this->seed != _non_)
        {
            dst1 = PQMath::l2d(tmpC1, 0, this->data, v, this->ndim);
            dst2 = PQMath::l2d(tmpC2, 0, this->data, v, this->ndim);
        }
        else
        {
            dst1 = rand();
            dst2 = rand();
        }
        loc  = v*this->ndim;

        if(dst1 < dst2)
        {
            for(j = 0; j < this->ndim; j++)
            {
                this->D1[j] += this->data[loc+j];
            }
            tmpn1++;
            nwlabels[i] = -1;
        }
        else
        {
            for(j = 0; j < this->ndim; j++)
            {
                D2[j] += this->data[loc+j];
            }
            tmpn2++;
            nwlabels[i] = v;
        }
    }

    if(tmpn1 == 0)
    {
        optEg = 0;
        n1 = tmpn2;
        n2 = tmpn1;
        E1 = 0;
        E2 = 0;
        _loc1 = nwlbl*this->ndim;
        memset(this->arrayD+_loc1, 0, cpysize);
        delete [] nwlabels;
        nwlabels = NULL;
        crntClust.clear();
        return false;
    }

    if(tmpn2 == 0)
    {
        optEg = 0;
        n1 = tmpn1;
        n2 = tmpn2;
        E2 = 0;
        E1 = 0;
        _loc1 = nwlbl*this->ndim;
        memset(this->arrayD+_loc1, 0, cpysize);
        delete [] nwlabels;
        nwlabels = NULL;

        crntClust.clear();
        return false;
    }
    innDcts[0] = PQMath::dvec_norm(D1, this->ndim, 2);
    innDcts[1] = PQMath::dvec_norm(D2, this->ndim, 2);
    innDcts[0] = innDcts[0]*innDcts[0];
    innDcts[1] = innDcts[1]*innDcts[1];

    optEg = sqrt(innDcts[0]) +  sqrt(innDcts[1]);

    int *rl = new int[numb];
    int nRuns = 0;

    /******** Incremental optimization *****/
    do
    {
        for(i = 0; i < numb; i++)
        {
            rl[i] = i;
        }
        random_shuffle(rl, rl + numb);
        _UPDATE_ = false;
        for(iter = 0; iter < numb; iter++)
        {
            mvs[0]++;
            memcpy(tmpD1, D1, cpysize);
            memcpy(tmpD2, D2, cpysize);
            r   = rl[iter];
            v   = crntClust[r];
            loc = v*this->ndim;
            F2S = false;

            if(nwlabels[r] == -1)
            {
                if(tmpn1 > 1)
                {
                    tmpn1--;
                    tmpn2++;
                    nwlabels[r] = v;
                    F2S = true;
                }
                else
                {
                    continue;
                }
            }
            else
            {
                if(tmpn2 > 1)
                {

                    tmpn1++;
                    tmpn2--;
                    nwlabels[r] = -1;
                    F2S = false;
                }
                else
                {
                    continue;
                }
            }

            tmp1 = tmp2 = 0;
            for(j = 0; j < this->ndim; j++)
            {
                tmp1 += tmpD1[j]*data[loc+j];
                tmp2 += tmpD2[j]*data[loc+j];
            }

            len = this->lens[v];
            if(F2S)
            {
                tmp1 = innDcts[0] - 2*tmp1 + len;
                tmp2 = innDcts[1] + 2*tmp2 + len;

            }
            else
            {
                tmp1 = innDcts[0] + 2*tmp1 + len;
                tmp2 = innDcts[1] - 2*tmp2 + len;
            }

            tmpEg = sqrt(tmp1) + sqrt(tmp2);

            if(tmpEg > (optEg+0.00001))
            {
                if(F2S)
                {
                    for(j = 0; j < this->ndim; j++)
                    {
                        tmpD1[j] -= this->data[loc+j];
                        tmpD2[j] += this->data[loc+j];
                    }
                }
                else
                {
                    for(j = 0; j < this->ndim; j++)
                    {
                        tmpD1[j] += this->data[loc+j];
                        tmpD2[j] -= this->data[loc+j];
                    }
                }
                memcpy(D1, tmpD1, cpysize);
                memcpy(D2, tmpD2, cpysize);
                optEg = tmpEg;
                innDcts[0] = tmp1;
                innDcts[1] = tmp2;
                _UPDATE_ = true;
                mvs[1]++;
            }
            else
            {
                ///roll-back if being updated
                if(nwlabels[r] == -1)
                {
                    tmpn1--;
                    tmpn2++;
                    nwlabels[r] = v;
                }
                else
                {
                    tmpn1++;
                    tmpn2--;
                    nwlabels[r] = -1;
                }
            }
        }
    }
    while(_UPDATE_);

    _loc1 = clabel*this->ndim;
    memcpy(this->arrayD+_loc1, D1, cpysize);
    _loc1 = nwlbl*this->ndim;
    memcpy(this->arrayD+_loc1, D2, cpysize);

    if(optEg > 0)
    {
        n1 = tmpn1;
        n2 = tmpn2;
    }

    for(j = 0; j < numb; j++)
    {
        labels[crntClust[j]] = clabel;
    }
    for(i = 0; i < numb; i++)
    {
        v = nwlabels[i];
        if(v != -1)
        {
            labels[v] = nwlbl;
        }
    }

    double len1, len2, tmpE1 = 0, tmpE2 = 0;
    len1 = len2 = 0;
    int k = 0;
    for(i = 0; i < numb; i++)
    {
        v = nwlabels[i];
        k = crntClust[i];
        if(v != -1)
        {
            len2 += lens[k];
        }
        else
        {
            len1 += lens[k];
        }
    }

    /**estimate average intra-cluster distance*/
    for(j = 0; j < this->ndim; j++)
    {
        C1[j] = D1[j]/n1;
        C2[j] = D2[j]/n2;
    }
    tmpE1 = PQMath::dvec_norm(C1, ndim, 2);
    tmpE2 = PQMath::dvec_norm(C2, ndim, 2);
    if(n1 > 1)
    {
        E1 = len1 - n1*tmpE1*tmpE1;
        E1 = 2*(E1/(n1-1));
    }
    else
    {
        E1 = 0;
    }

    if(n2 > 1)
    {
        E2 = len2 - n2*tmpE2*tmpE2;
        E2 = 2*(E2/(n2-1));
    }
    else
    {
        E2 = 0;
    }

    crntClust.clear();
    delete [] nwlabels;
    delete [] rl;
    nwlabels = NULL;
    return true;
}

bool XBKMeans::incrOptzI1(const int clabel, double &E1, double &E2, int &n1, int &n2, const int nwlbl)
{
    unsigned int _loc1 = 0, _loc2 = 0, loc = 0;
    unsigned int i = 0, j = 0, r = 0;
    unsigned int cpysize = sizeof(double)*this->ndim;
    double tmp1 = 0, tmp2 = 0, len = 0, x = 0, y = 0, tmpx = 0, tmpy = 0, tmpx1 = 0, tmpy1 = 0;
    double tmp11 = 0, tmp21 = 0;

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
        E2 = E1 = 0.0f;
        n1 = n2 = 1;
        i  = crntClust[1];
        this->labels[i] = nwlbl;

        _loc1 = crntClust[0]*this->ndim;
        _loc2 = crntClust[1]*this->ndim;
        for(j = 0; j < this->ndim; j++)
        {
            arrayD[clabel*this->ndim+j] = this->data[_loc1+j];
            arrayD[nwlbl*this->ndim+j]  = this->data[_loc2+j];
        }

        crntClust.clear();

        return true;
    }

    int *nwlabels = new int[numb];
    memset(nwlabels, 0, sizeof(int)*numb);

    int v = 0, tmpn1 = 0, tmpn2 = 0;
    unsigned int sed1 = 0, sed2 = 0, iter = 0;
    double tmpEg = 0, optEg = 0, dst1 = 0,  dst2 = 0;
    bool _UPDATE_ = true, F2S = false;
    n1 = n2 = 0;
    /**** Initial assigment ****************/
    if(this->seed != _non_)
    {
        (this->*sedFunc)(sed1, sed2, crntClust);
        _loc1 = sed1*this->ndim;
        _loc2 = sed2*this->ndim;
        for(j = 0; j < this->ndim; j++)
        {
            tmpC1[j] = this->data[_loc1+j];
            tmpC2[j] = this->data[_loc2+j];
        }
    }

    memset(D1, 0, cpysize);
    memset(D2, 0, cpysize);
    innDcts[0] = innDcts[1] = 0;

    for(i = 0; i < numb; i++)
    {
        v    = crntClust[i];
        if(this->seed != _non_)
        {
            dst1 = PQMath::l2d(tmpC1, 0, this->data, v, this->ndim);
            dst2 = PQMath::l2d(tmpC2, 0, this->data, v, this->ndim);
        }
        else
        {
            dst1 = rand();
            dst2 = rand();
        }
        loc  = v*this->ndim;

        if(dst1 < dst2)
        {
            for(j = 0; j < this->ndim; j++)
            {
                this->D1[j] += this->data[loc+j];
            }
            tmpn1++;
            nwlabels[i] = -1;
        }
        else
        {
            for(j = 0; j < this->ndim; j++)
            {
                D2[j] += this->data[loc+j];
            }
            tmpn2++;
            nwlabels[i] = v;
        }
    }

    if(tmpn1 == 0)
    {
        optEg = 0;
        n1 = tmpn2;
        n2 = tmpn1;
        E1 = 0;
        E2 = 0;
        _loc1 = nwlbl*this->ndim;
        memset(this->arrayD+_loc1, 0, cpysize);
        delete [] nwlabels;
        nwlabels = NULL;
        crntClust.clear();
        return false;
    }

    if(tmpn2 == 0)
    {
        optEg = 0;
        n1 = tmpn1;
        n2 = tmpn2;
        E2 = 0;
        E1 = 0;
        _loc1 = nwlbl*this->ndim;
        memset(this->arrayD+_loc1, 0, cpysize);
        delete [] nwlabels;
        nwlabels = NULL;
        crntClust.clear();
        return false;
    }
    innDcts[0] = PQMath::dvec_norm(D1, this->ndim, 2);
    innDcts[1] = PQMath::dvec_norm(D2, this->ndim, 2);
    innDcts[0] = innDcts[0]*innDcts[0];
    innDcts[1] = innDcts[1]*innDcts[1];
    x = 0;
    y = 0;
    for(i = 0; i < numb; i++)
    {
        if(nwlabels[i] == -1)
            x += this->lens[crntClust[i]];
        else
            y += this->lens[crntClust[i]];
    }

    if(tmpn1 != 1)
    {
        tmpx1 = x*tmpn1/(tmpn1 - 1);
        tmp11 = innDcts[0]/(tmpn1-1);
    }
    else
    {
        tmpx1 = len;
        tmp11 = 0;
    }
    if(tmpn2 != 1)
    {
        tmpy1 = y*tmpn2/(tmpn2 - 1);
        tmp21 = innDcts[1]/(tmpn2-1);
    }
    else
    {
        tmpy1 = len;
        tmp21 = 0;
    }

    optEg = tmp11 + tmp21 - tmpx1 - tmpy1;

    /******** Incremental optimization *****/
    int *rl = new int[numb];

    do
    {
        for(i = 0; i < numb; i++)
        {
            rl[i] = i;
        }
        random_shuffle(rl, rl + numb);
        _UPDATE_ = false;
        for(iter = 0; iter < numb; iter++)
        {
            mvs[0]++;
            memcpy(tmpD1, D1, cpysize);
            memcpy(tmpD2, D2, cpysize);
            ///r   = (unsigned int)floor((rand()/(RAND_MAX+1.0f))*numb);
            ///r   = r>=numb?(numb-1):r;
            r   = rl[iter];
            v   = crntClust[r];
            loc = v*this->ndim;
            F2S = false;

            if(nwlabels[r] == -1)
            {
                if(tmpn1 > 1)
                {
                    tmpn1--;
                    tmpn2++;
                    nwlabels[r] = v;
                    F2S = true;
                }
                else
                {
                    continue;
                }
            }
            else
            {
                if(tmpn2 > 1)
                {

                    tmpn1++;
                    tmpn2--;
                    nwlabels[r] = -1;
                    F2S = false;
                }
                else
                {
                    continue;
                }
            }

            tmp1 = tmp2 = 0;
            for(j = 0; j < this->ndim; j++)
            {
                tmp1 += tmpD1[j]*data[loc+j];
                tmp2 += tmpD2[j]*data[loc+j];
            }

            ///tmpEg = (this->*optFunc)(tmpD1, tmpD2, this->ndim, tmpn1, tmpn2);
            len = this->lens[v];
            if(F2S)
            {
                tmp1 = innDcts[0] - 2*tmp1 + len;
                tmp2 = innDcts[1] + 2*tmp2 + len;
                tmpx = x - len;
                tmpy = y + len;
                if(tmpn1 != 1)
                {
                    tmpx1 = tmpx*tmpn1/(tmpn1 - 1);
                    tmp11 = tmp1/(tmpn1-1);
                }
                else
                {
                    tmpx1 = len;
                    tmp11 = 0;
                }
                if(tmpn2 != 1)
                {
                    tmpy1 = tmpy*tmpn2/(tmpn2 - 1);
                    tmp21 = tmp2/(tmpn2 - 1);
                }
                else
                {
                    tmpy1 = len;
                    tmp21 = 0;
                }
            }
            else
            {
                tmp1 = innDcts[0] + 2*tmp1 + len;
                tmp2 = innDcts[1] - 2*tmp2 + len;
                tmpx = x + len;
                tmpy = y - len;
                if(tmpn1 != 1)
                {
                    tmpx1 = tmpx*tmpn1/(tmpn1 - 1);
                    tmp11 = tmp1/(tmpn1-1);
                }
                else
                {
                    tmpx1 = len;
                    tmp11 = 0;
                }
                if(tmpn2 != 1)
                {
                    tmpy1 = tmpy*tmpn2/(tmpn2 - 1);
                    tmp21 = tmp2/(tmpn2 - 1);
                }
                else
                {
                    tmpy1 = len;
                    tmp21 = 0;
                }
            }

            tmpEg = tmp11 + tmp21 - tmpx1 - tmpy1;

            if(tmpEg > (optEg+0.00001))
            {
                if(F2S)
                {
                    for(j = 0; j < this->ndim; j++)
                    {
                        tmpD1[j] -= this->data[loc+j];
                        tmpD2[j] += this->data[loc+j];
                    }
                }
                else
                {
                    for(j = 0; j < this->ndim; j++)
                    {
                        tmpD1[j] += this->data[loc+j];
                        tmpD2[j] -= this->data[loc+j];
                    }
                }
                memcpy(D1, tmpD1, cpysize);
                memcpy(D2, tmpD2, cpysize);
                optEg = tmpEg;
                x = tmpx;
                y = tmpy;
                innDcts[0] = tmp1;
                innDcts[1] = tmp2;
                _UPDATE_ = true;
                mvs[1]++;
            }
            else
            {
                ///roll-back if being updated
                if(nwlabels[r] == -1)
                {
                    tmpn1--;
                    tmpn2++;
                    nwlabels[r] = v;
                }
                else
                {
                    tmpn1++;
                    tmpn2--;
                    nwlabels[r] = -1;
                }
            }
        }
    }
    while(_UPDATE_);
    _loc1 = clabel*this->ndim;
    memcpy(this->arrayD+_loc1, D1, cpysize);
    _loc1 = nwlbl*this->ndim;
    memcpy(this->arrayD+_loc1, D2, cpysize);

    // if(optEg > 0)
    {
        n1 = tmpn1;
        n2 = tmpn2;
    }

    for(i = 0; i < numb; i++)
    {
        v = nwlabels[i];
        if(v != -1)
        {
            labels[v] = nwlbl;
        }
    }

    double len1, len2, tmpE1 = 0, tmpE2 = 0;
    len1 = len2 = 0;
    int k = 0;
    for(i = 0; i < numb; i++)
    {
        v = nwlabels[i];
        k = crntClust[i];
        if(v != -1)
        {
            len2 += lens[k];
        }
        else
        {
            len1 += lens[k];
        }
    }

    /**estimate average intra-cluster distance*/
    for(j = 0; j < this->ndim; j++)
    {
        C1[j] = D1[j]/n1;
        C2[j] = D2[j]/n2;
    }
    tmpE1 = PQMath::dvec_norm(C1, ndim, 2);
    tmpE2 = PQMath::dvec_norm(C2, ndim, 2);
    if(n1 > 1)
    {
        E1 = len1 - n1*tmpE1*tmpE1;
        E1 = 2*(E1/(n1-1));
    }
    else
    {
        E1 = 0;
    }

    if(n2 > 1)
    {
        E2 = len2 - n2*tmpE2*tmpE2;
        E2 = 2*(E2/(n2-1));
    }
    else
    {
        E2 = 0;
    }

    crntClust.clear();
    delete [] nwlabels;
    delete [] rl;
    nwlabels = NULL;
    return true;
}


bool XBKMeans::incrOptzE1(const int clabel, double &E1, double &E2, int &n1, int &n2, const int nwlbl)
{
    unsigned int _loc1 = 0, _loc2 = 0, loc = 0;
    unsigned int i = 0, j = 0, r = 0;
    unsigned int cpysize = sizeof(double)*this->ndim;
    double tmp1 = 0, tmp2 = 0, len = 0, tmpij = 0, x = 0, y = 0, tmpx = 0, tmpy = 0;
    double innerij = 0;

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
        E2 = E1 = 0.0f;
        n1 = n2 = 1;
        i  = crntClust[1];
        this->labels[i] = nwlbl;
        _loc1 = crntClust[0]*this->ndim;
        _loc2 = crntClust[1]*this->ndim;
        for(j = 0; j < this->ndim; j++)
        {
            arrayD[clabel*this->ndim+j] = this->data[_loc1+j];
            arrayD[nwlbl*this->ndim+j]  = this->data[_loc2+j];
        }
        crntClust.clear();
        return true;
    }

    int *nwlabels = new int[numb];
    memset(nwlabels, 0, sizeof(int)*numb);

    int v = 0, tmpn1 = 0, tmpn2 = 0;
    unsigned int sed1 = 0, sed2 = 0, iter = 0;
    double tmpEg = 0, optEg = 0, dst1 = 0,  dst2 = 0;
    bool _UPDATE_ = true, F2S = false;
    n1 = n2 = 0;
    /**** Initial assigment ****************/

    if(this->seed != _non_)
    {
        (this->*sedFunc)(sed1, sed2, crntClust);
        _loc1 = sed1*this->ndim;
        _loc2 = sed2*this->ndim;
        for(j = 0; j < this->ndim; j++)
        {
            tmpC1[j] = this->data[_loc1+j];
            tmpC2[j] = this->data[_loc2+j];
        }
    }
    memset(D1, 0, cpysize);
    memset(D2, 0, cpysize);
    innDcts[0] = innDcts[1] = 0;

    for(i = 0; i < numb; i++)
    {
        v    = crntClust[i];
        if(this->seed == _non_)
        {
            dst1 = PQMath::l2d(tmpC1, 0, this->data, v, this->ndim);
            dst2 = PQMath::l2d(tmpC2, 0, this->data, v, this->ndim);
        }
        else
        {
            dst1 = rand();
            dst2 = rand();
        }
        loc  = v*this->ndim;

        if(dst1 < dst2)
        {
            for(j = 0; j < this->ndim; j++)
            {
                this->D1[j] += this->data[loc+j];
            }
            tmpn1++;
            nwlabels[i] = -1;
        }
        else
        {
            for(j = 0; j < this->ndim; j++)
            {
                D2[j] += this->data[loc+j];
            }
            tmpn2++;
            nwlabels[i] = v;
        }
    }

    if(tmpn1 == 0)
    {
        optEg = 0;
        n1 = tmpn2;
        n2 = tmpn1;
        E1 = 0;
        E2 = 0;
        _loc1 = nwlbl*this->ndim;
        memset(this->arrayD+_loc1, 0, cpysize);
        delete [] nwlabels;
        nwlabels = NULL;
        crntClust.clear();
        return false;
    }


    if(tmpn2 == 0)
    {
        optEg = 0;
        n1 = tmpn1;
        n2 = tmpn2;
        E2 = 0;
        E1 = 0;
        _loc1 = nwlbl*this->ndim;
        memset(this->arrayD+_loc1, 0, cpysize);
        delete [] nwlabels;
        nwlabels = NULL;
        crntClust.clear();
        return false;
    }
    innDcts[0] = PQMath::dvec_norm(D1, this->ndim, 2);
    innDcts[1] = PQMath::dvec_norm(D2, this->ndim, 2);
    innDcts[0] = innDcts[0]*innDcts[0];
    innDcts[1] = innDcts[1]*innDcts[1];
    innerij = PQMath::innerProduct(D1, D2, this->ndim);
    x = 0;
    y = 0;
    for(i = 0; i < numb; i++)
    {
        if(nwlabels[i] == -1)
            x += this->lens[crntClust[i]];
        else
            y += this->lens[crntClust[i]];
    }

    optEg = (tmpn1*x + tmpn2*y - 2*innerij)*tmpn1*tmpn2;

    /******** Incremental optimization *****/

    int *rl = new int[numb];
    do
    {
        for(i = 0; i < numb; i++)
        {
            rl[i] = i;
        }
        random_shuffle(rl, rl + numb);
        _UPDATE_ = false;
        for(iter = 0; iter < numb; iter++)
        {
            mvs[0]++;
            memcpy(tmpD1, D1, cpysize);
            memcpy(tmpD2, D2, cpysize);
            r   = (unsigned int)floor((rand()/(RAND_MAX+1.0f))*numb);
            r   = r>=numb?(numb-1):r;
            //r = rl[iter];
            v   = crntClust[r];
            loc = v*this->ndim;
            F2S = false;

            if(nwlabels[r] == -1)
            {
                if(tmpn1 > 1)
                {
                    tmpn1--;
                    tmpn2++;
                    nwlabels[r] = v;
                    F2S = true;
                }
                else
                {
                    continue;
                }
            }
            else
            {
                if(tmpn2 > 1)
                {

                    tmpn1++;
                    tmpn2--;
                    nwlabels[r] = -1;
                    F2S = false;
                }
                else
                {
                    continue;
                }
            }

            tmp1 = tmp2 = 0;
            for(j = 0; j < this->ndim; j++)
            {
                tmp1 += tmpD1[j]*data[loc+j];
                tmp2 += tmpD2[j]*data[loc+j];
            }

            ///tmpEg = (this->*optFunc)(tmpD1, tmpD2, this->ndim, tmpn1, tmpn2);
            len = this->lens[v];
            if(F2S)
            {
                tmpij = innerij + tmp1 - tmp2 - len;
                tmpx = x - len;
                tmpy = y + len;
            }
            else
            {
                tmpij = innerij - tmp1 + tmp2 - len;
                tmpx = x + len;
                tmpy = y - len;
            }

            tmpEg = (tmpn1*tmpx + tmpn2*tmpy - 2*tmpij)*tmpn1*tmpn2;

            if(tmpEg > (optEg+0.00001))
            {
                if(F2S)
                {
                    for(j = 0; j < this->ndim; j++)
                    {
                        tmpD1[j] -= this->data[loc+j];
                        tmpD2[j] += this->data[loc+j];
                    }
                }
                else
                {
                    for(j = 0; j < this->ndim; j++)
                    {
                        tmpD1[j] += this->data[loc+j];
                        tmpD2[j] -= this->data[loc+j];
                    }
                }
                memcpy(D1, tmpD1, cpysize);
                memcpy(D2, tmpD2, cpysize);
                innerij  = tmpij;
                optEg = tmpEg;
                x = tmpx;
                y = tmpy;
                _UPDATE_ = true;
                mvs[1]++;
            }
            else
            {
                ///roll-back if being updated
                if(nwlabels[r] == -1)
                {
                    tmpn1--;
                    tmpn2++;
                    nwlabels[r] = v;
                }
                else
                {
                    tmpn1++;
                    tmpn2--;
                    nwlabels[r] = -1;
                }
            }
        }
    }
    while(_UPDATE_);
    _loc1 = clabel*this->ndim;
    memcpy(this->arrayD+_loc1, D1, cpysize);
    _loc1 = nwlbl*this->ndim;
    memcpy(this->arrayD+_loc1, D2, cpysize);

    // if(optEg > 0)
    {
        n1 = tmpn1;
        n2 = tmpn2;
    }

    for(i = 0; i < numb; i++)
    {
        v = nwlabels[i];
        if(v != -1)
        {
            labels[v] = nwlbl;
        }
    }

    double len1, len2, tmpE1 = 0, tmpE2 = 0;
    len1 = len2 = 0;
    int k = 0;
    for(i = 0; i < numb; i++)
    {
        v = nwlabels[i];
        k = crntClust[i];
        if(v != -1)
        {
            len2 += lens[k];
        }
        else
        {
            len1 += lens[k];
        }
    }

    /**estimate average intra-cluster distance*/
    for(j = 0; j < this->ndim; j++)
    {
        C1[j] = D1[j]/n1;
        C2[j] = D2[j]/n2;
    }
    tmpE1 = PQMath::dvec_norm(C1, ndim, 2);
    tmpE2 = PQMath::dvec_norm(C2, ndim, 2);
    if(n1 > 1)
    {
        E1 = len1 - n1*tmpE1*tmpE1;
        E1 = 2*(E1/(n1-1));
    }
    else
    {
        E1 = 0;
    }

    if(n2 > 1)
    {
        E2 = len2 - n2*tmpE2*tmpE2;
        E2 = 2*(E2/(n2-1));
    }
    else
    {
        E2 = 0;
    }

    crntClust.clear();
    delete [] nwlabels;
    delete [] rl;
    nwlabels = NULL;
    return true;
}

bool XBKMeans::incrOptzE2(const int clabel, double &E1, double &E2, int &n1, int &n2, const int nwlbl)
{
    unsigned int _loc1 = 0, _loc2 = 0, loc = 0;
    unsigned int i = 0, j = 0, r = 0;
    unsigned int cpysize = sizeof(double)*this->ndim;
    double tmp1 = 0, tmp2 = 0, len = 0, tmpij = 0;
    double innerij ;

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
        E2 = E1 = 0.0f;
        n1 = n2 = 1;
        i  = crntClust[1];
        this->labels[i] = nwlbl;

        _loc1 = crntClust[0]*this->ndim;
        _loc2 = crntClust[1]*this->ndim;
        for(j = 0; j < this->ndim; j++)
        {
            arrayD[clabel*this->ndim+j] = this->data[_loc1+j];
            arrayD[nwlbl*this->ndim+j]  = this->data[_loc2+j];
        }
        crntClust.clear();
        return true;
    }

    int *nwlabels = new int[numb];
    memset(nwlabels, 0, sizeof(int)*numb);

    int v = 0, tmpn1 = 0, tmpn2 = 0;
    unsigned int sed1 = 0, sed2 = 0, iter = 0;
    double tmpEg = 0, optEg = 0, dst1 = 0,  dst2 = 0;
    bool _UPDATE_ = true, F2S = false;
    n1 = n2 = 0;
    /**** Initial assigment ****************/
    if(this->seed != _non_)
    {
        (this->*sedFunc)(sed1, sed2, crntClust);
        _loc1 = sed1*this->ndim;
        _loc2 = sed2*this->ndim;
        for(j = 0; j < this->ndim; j++)
        {
            tmpC1[j] = this->data[_loc1+j];
            tmpC2[j] = this->data[_loc2+j];
        }
    }

    memset(D1, 0, cpysize);
    memset(D2, 0, cpysize);
    innDcts[0] = innDcts[1] = 0;

    for(i = 0; i < numb; i++)
    {
        v    = crntClust[i];
        if(this->seed == _non_)
        {
            dst1 = PQMath::l2d(tmpC1, 0, this->data, v, this->ndim);
            dst2 = PQMath::l2d(tmpC2, 0, this->data, v, this->ndim);
        }
        else
        {
            dst1 = rand();
            dst2 = rand();
        }
        loc  = v*this->ndim;

        if(dst1 < dst2)
        {
            for(j = 0; j < this->ndim; j++)
            {
                this->D1[j] += this->data[loc+j];
            }
            tmpn1++;
            nwlabels[i] = -1;
        }
        else
        {
            for(j = 0; j < this->ndim; j++)
            {
                D2[j] += this->data[loc+j];
            }
            tmpn2++;
            nwlabels[i] = v;
        }
    }

    if(tmpn1 == 0)
    {
        optEg = 0;
        n1 = tmpn2;
        n2 = tmpn1;
        E1 = 0;
        E2 = 0;
        _loc1 = nwlbl*this->ndim;
        memset(this->arrayD+_loc1, 0, cpysize);
        delete [] nwlabels;
        nwlabels = NULL;
        crntClust.clear();
        return false;
    }


    if(tmpn2 == 0)
    {
        optEg = 0;
        n1 = tmpn1;
        n2 = tmpn2;
        E2 = 0;
        E1 = 0;
        _loc1 = nwlbl*this->ndim;
        memset(this->arrayD+_loc1, 0, cpysize);
        delete [] nwlabels;
        nwlabels = NULL;
        crntClust.clear();
        return false;
    }
    innDcts[0] = PQMath::dvec_norm(D1, this->ndim, 2);
    innDcts[1] = PQMath::dvec_norm(D2, this->ndim, 2);
    innDcts[0] = innDcts[0]*innDcts[0];
    innDcts[1] = innDcts[1]*innDcts[1];
    innerij = PQMath::innerProduct(D1, D2, this->ndim);

    optEg = innDcts[0]/tmpn1*tmpn2 + innDcts[1]/tmpn2*tmpn1 - 2*innerij;

    /******** Incremental optimization *****/

    int *rl = new int[numb];
    do
    {
        for(i = 0; i < numb; i++)
        {
            rl[i] = i;
        }
        random_shuffle(rl, rl + numb);
        _UPDATE_ = false;
        for(iter = 0; iter < numb; iter++)
        {
            mvs[0]++;
            memcpy(tmpD1, D1, cpysize);
            memcpy(tmpD2, D2, cpysize);
            r   = (unsigned int)floor((rand()/(RAND_MAX+1.0f))*numb);
            r   = r>=numb?(numb-1):r;
            //r = rl[iter];
            v   = crntClust[r];
            loc = v*this->ndim;
            F2S = false;

            if(nwlabels[r] == -1)
            {
                if(tmpn1 > 1)
                {
                    tmpn1--;
                    tmpn2++;
                    nwlabels[r] = v;
                    F2S = true;
                }
                else
                {
                    continue;
                }
            }
            else
            {
                if(tmpn2 > 1)
                {

                    tmpn1++;
                    tmpn2--;
                    nwlabels[r] = -1;
                    F2S = false;
                }
                else
                {
                    continue;
                }
            }

            tmp1 = tmp2 = 0;
            for(j = 0; j < this->ndim; j++)
            {
                tmp1 += tmpD1[j]*data[loc+j];
                tmp2 += tmpD2[j]*data[loc+j];
            }

            ///tmpEg = (this->*optFunc)(tmpD1, tmpD2, this->ndim, tmpn1, tmpn2);
            len = this->lens[v];
            if(F2S)
            {
                tmpij = innerij + tmp1 - tmp2 - len;
                tmp1 = innDcts[0] - 2*tmp1 + len;
                tmp2 = innDcts[1] + 2*tmp2 + len;
            }
            else
            {
                tmpij = innerij - tmp1 + tmp2 - len;
                tmp1 = innDcts[0] + 2*tmp1 + len;
                tmp2 = innDcts[1] - 2*tmp2 + len;

            }

            tmpEg = tmp1/tmpn1*tmpn2 + tmp2/tmpn2*tmpn1 - 2*tmpij;

            if(tmpEg > (optEg+0.00001))
            {
                if(F2S)
                {
                    for(j = 0; j < this->ndim; j++)
                    {
                        tmpD1[j] -= this->data[loc+j];
                        tmpD2[j] += this->data[loc+j];
                    }
                }
                else
                {
                    for(j = 0; j < this->ndim; j++)
                    {
                        tmpD1[j] += this->data[loc+j];
                        tmpD2[j] -= this->data[loc+j];
                    }
                }
                memcpy(D1, tmpD1, cpysize);
                memcpy(D2, tmpD2, cpysize);
                innerij  = tmpij;
                optEg = tmpEg;
                innDcts[0] = tmp1;
                innDcts[1] = tmp2;
                _UPDATE_ = true;
                mvs[1]++;
            }
            else
            {
                ///roll-back if being updated
                if(nwlabels[r] == -1)
                {
                    tmpn1--;
                    tmpn2++;
                    nwlabels[r] = v;
                }
                else
                {
                    tmpn1++;
                    tmpn2--;
                    nwlabels[r] = -1;
                }
            }
        }
    }
    while(_UPDATE_);
    _loc1 = clabel*this->ndim;
    memcpy(this->arrayD+_loc1, D1, cpysize);
    _loc1 = nwlbl*this->ndim;
    memcpy(this->arrayD+_loc1, D2, cpysize);

    if(optEg > 0)
    {
        n1 = tmpn1;
        n2 = tmpn2;
    }

    for(i = 0; i < numb; i++)
    {
        v = nwlabels[i];
        if(v != -1)
        {
            labels[v] = nwlbl;
        }
    }

    double len1, len2, tmpE1 = 0, tmpE2 = 0;
    len1 = len2 = 0;
    int k = 0;
    for(i = 0; i < numb; i++)
    {
        v = nwlabels[i];
        k = crntClust[i];
        if(v != -1)
        {
            len2 += lens[k];
        }
        else
        {
            len1 += lens[k];
        }
    }

    /**estimate average intra-cluster distance*/
    for(j = 0; j < this->ndim; j++)
    {
        C1[j] = D1[j]/n1;
        C2[j] = D2[j]/n2;
    }
    tmpE1 = PQMath::dvec_norm(C1, ndim, 2);
    tmpE2 = PQMath::dvec_norm(C2, ndim, 2);
    if(n1 > 1)
    {
        E1 = len1 - n1*tmpE1*tmpE1;
        E1 = 2*(E1/(n1-1));
    }
    else
    {
        E1 = 0;
    }

    if(n2 > 1)
    {
        E2 = len2 - n2*tmpE2*tmpE2;
        E2 = 2*(E2/(n2-1));
    }
    else
    {
        E2 = 0;
    }

    crntClust.clear();
    delete [] nwlabels;
    delete [] rl;
    nwlabels = NULL;
    return true;
}

bool XBKMeans::incrOptzI2s(const int clabel, double &E1, double &E2, int &n1, int &n2, const int nwlbl)
{
    unsigned int _loc1 = 0, i = 0, j = 0, r = 0;
    unsigned int cpysize = sizeof(double)*this->ndim;
    double tmp1 = 0, tmp2 = 0, len = 0;
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
        E2 = E1 = 0.0f;
        n1 = n2 = 1;
        i  = crntClust[1];
        this->labels[i] = nwlbl;
        crntClust.clear();

        return true;
    }
    int *nwlabels = new int[numb];
    memset(nwlabels, 0, sizeof(int)*numb);
    int v = 0, tmpn1 = 0, tmpn2 = 0;
    unsigned int sed1 = 0, iter = 0, sed2 = 0;
    double tmpEg = 0, optEg = 0, dst1 = 0,  dst2 = 0;
    bool _UPDATE_ = true, F2S = false;
    n1 = n2 = 0;
    /**** Initial assigment ****************/
    if(this->seed != _non_)
    {
        (this->*sedFunc)(sed1, sed2, crntClust);
        memset(tmpC1, 0, sizeof(double)*this->ndim);
        memset(tmpC2, 0, sizeof(double)*this->ndim);
        for(j = sdata.col[sed1]; j < sdata.col[sed1 + 1]; j++)
        {
            tmpC1[sdata.index[j]] = sdata.data[j];
        }
        for(j = sdata.col[sed2]; j < sdata.col[sed2 + 1]; j++)
        {
            tmpC2[sdata.index[j]] = sdata.data[j];
        }
    }

    memset(D1, 0, cpysize);
    memset(D2, 0, cpysize);
    innDcts[0] = innDcts[1] = 0;
    for(i = 0; i < numb; i++)
    {
        v    = crntClust[i];
        if(this->seed != _non_)
        {
            dst1 = PQMath::l2d(tmpC1, 0, sdata, v, this->ndim);
            dst2 = PQMath::l2d(tmpC2, 0, sdata, v, this->ndim);
        }
        else
        {
            dst1 = rand();
            dst2 = rand();
        }

        if(dst1 < dst2)
        {
            for(j = sdata.col[v]; j < sdata.col[v+1]; j++)
            {
                this->D1[sdata.index[j]] += sdata.data[j];
            }
            tmpn1++;
            nwlabels[i] = -1;
        }
        else
        {
            for(j = sdata.col[v]; j < sdata.col[v+1]; j++)
            {
                this->D2[sdata.index[j]] += sdata.data[j];
            }
            tmpn2++;
            nwlabels[i] = v;
        }
    }
    if(tmpn1 == 0)
    {
        optEg = 0;
        n1 = tmpn2;
        n2 = tmpn1;
        E1 = 0;
        E2 = 0;
        _loc1 = nwlbl*this->ndim;
        memset(this->arrayD+_loc1, 0, cpysize);
        delete [] nwlabels;
        nwlabels = NULL;
        crntClust.clear();
        return false;
    }


    if(tmpn2 == 0)
    {
        optEg = 0;
        n1 = tmpn1;
        n2 = tmpn2;
        E2 = 0;
        E1 = 0;
        _loc1 = nwlbl*this->ndim;
        memset(this->arrayD+_loc1, 0, cpysize);
        delete [] nwlabels;
        nwlabels = NULL;
        crntClust.clear();
        return false;
    }
    innDcts[0] = PQMath::dvec_norm(D1, this->ndim, 2);
    innDcts[1] = PQMath::dvec_norm(D2, this->ndim, 2);
    innDcts[0] = innDcts[0]*innDcts[0];
    innDcts[1] = innDcts[1]*innDcts[1];

    optEg = innDcts[0]/tmpn1 +  innDcts[1]/tmpn2;


    int *rl = new int[numb];
    /******** Incremental optimization *****/
    do
    {
        for(i = 0; i < numb; i++)
        {
            rl[i] = i;
        }
        random_shuffle(rl, rl + numb);
        _UPDATE_ = false;
        for(iter = 0; iter < numb; iter++)
        {
            mvs[0]++;
            memcpy(tmpD1, D1, cpysize);
            memcpy(tmpD2, D2, cpysize);
            r   = (unsigned int)floor((rand()/(RAND_MAX+1.0f))*numb);
            r   = r>=numb?(numb-1):r;
            //r = rl[iter];
            v   = crntClust[r];
            F2S = false;

            if(nwlabels[r] == -1)
            {
                if(tmpn1 > 1)
                {
                    tmpn1--;
                    tmpn2++;
                    nwlabels[r] = v;
                    F2S = true;
                }
                else
                {
                    continue;
                }
            }
            else
            {
                if(tmpn2 > 1)
                {

                    tmpn1++;
                    tmpn2--;
                    nwlabels[r] = -1;
                    F2S = false;
                }
                else
                {
                    continue;
                }
            }

            tmp1 = tmp2 = 0;
            for(j = sdata.col[v]; j < sdata.col[v+1]; j++)
            {
                tmp1 += tmpD1[sdata.index[j]]*sdata.data[j];
                tmp2 += tmpD2[sdata.index[j]]*sdata.data[j];
            }

            ///tmpEg = (this->*optFunc)(tmpD1, tmpD2, this->ndim, tmpn1, tmpn2);
            len = this->lens[v];
            if(F2S)
            {
                tmp1 = innDcts[0] - 2*tmp1 + len;
                tmp2 = innDcts[1] + 2*tmp2 + len;

            }
            else
            {
                tmp1 = innDcts[0] + 2*tmp1 + len;
                tmp2 = innDcts[1] - 2*tmp2 + len;
            }

            tmpEg = tmp1/tmpn1 + tmp2/tmpn2;

            if(tmpEg > (optEg+0.00001))
            {
                if(F2S)
                {
                    for(j = sdata.col[v]; j < sdata.col[v+1]; j++)
                    {
                        tmpD1[sdata.index[j]] -= sdata.data[j];
                        tmpD2[sdata.index[j]] += sdata.data[j];
                    }
                }
                else
                {
                    for(j = sdata.col[v]; j < sdata.col[v+1]; j++)
                    {
                        tmpD1[sdata.index[j]] += sdata.data[j];
                        tmpD2[sdata.index[j]] -= sdata.data[j];
                    }
                }
                memcpy(D1, tmpD1, cpysize);
                memcpy(D2, tmpD2, cpysize);
                optEg = tmpEg;
                innDcts[0] = tmp1;
                innDcts[1] = tmp2;
                _UPDATE_ = true;
                mvs[1]++;
            }
            else
            {
                ///roll-back if being updated
                if(nwlabels[r] == -1)
                {
                    tmpn1--;
                    tmpn2++;
                    nwlabels[r] = v;
                }
                else
                {
                    tmpn1++;
                    tmpn2--;
                    nwlabels[r] = -1;
                }
            }
        }
    }
    while(_UPDATE_);

    _loc1 = clabel*this->ndim;
    memcpy(this->arrayD+_loc1, D1, cpysize);
    _loc1 = nwlbl*this->ndim;
    memcpy(this->arrayD+_loc1, D2, cpysize);

    if(optEg > 0)
    {
        n1 = tmpn1;
        n2 = tmpn2;
    }

    for(i = 0; i < numb; i++)
    {
        v = nwlabels[i];
        if(v != -1)
        {
            labels[v] = nwlbl;
        }
    }

    double len1, len2, tmpE1 = 0, tmpE2 = 0;
    len1 = len2 = 0;
    int k = 0;
    for(i = 0; i < numb; i++)
    {
        v = nwlabels[i];
        k = crntClust[i];
        if(v != -1)
        {
            len2 += lens[k];
        }
        else
        {
            len1 += lens[k];
        }
    }

    /**estimate average intra-cluster distance*/
    for(j = 0; j < this->ndim; j++)
    {
        C1[j] = D1[j]/n1;
        C2[j] = D2[j]/n2;
    }
    tmpE1 = PQMath::dvec_norm(C1, ndim, 2);
    tmpE2 = PQMath::dvec_norm(C2, ndim, 2);
    if(n1 > 1)
    {
        E1 = len1 - n1*tmpE1*tmpE1;
        E1 = 2*(E1/(n1-1));
    }
    else
    {
        E1 = 0;
    }

    if(n2 > 1)
    {
        E2 = len2 - n2*tmpE2*tmpE2;
        E2 = 2*(E2/(n2-1));
    }
    else
    {
        E2 = 0;
    }

    crntClust.clear();
    delete [] nwlabels;
    delete [] rl;
    nwlabels = NULL;
    return true;
}

bool XBKMeans::incrOptzI4s(const int clabel, double &E1, double &E2, int &n1, int &n2, const int nwlbl)
{
    unsigned int _loc1 = 0, i = 0, j = 0, r = 0;
    unsigned int cpysize = sizeof(double)*this->ndim;
    double tmp1 = 0, tmp2 = 0, len = 0;
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
        E2 = E1 = 0.0f;
        n1 = n2 = 1;
        i  = crntClust[1];
        this->labels[i] = nwlbl;
        crntClust.clear();

        return true;
    }
    int *nwlabels = new int[numb];
    memset(nwlabels, 0, sizeof(int)*numb);
    int v = 0, tmpn1 = 0, tmpn2 = 0;
    unsigned int sed1 = 0, iter = 0, sed2 = 0;
    double tmpEg = 0, optEg = 0, dst1 = 0,  dst2 = 0;
    bool _UPDATE_ = true, F2S = false;
    n1 = n2 = 0;
    /**** Initial assigment ****************/
    if(this->seed != _non_)
    {
        (this->*sedFunc)(sed1, sed2, crntClust);
        memset(tmpC1, 0, sizeof(double)*this->ndim);
        memset(tmpC2, 0, sizeof(double)*this->ndim);
        for(j = sdata.col[sed1]; j < sdata.col[sed1 + 1]; j++)
        {
            tmpC1[sdata.index[j]] = sdata.data[j];
        }
        for(j = sdata.col[sed2]; j < sdata.col[sed2 + 1]; j++)
        {
            tmpC2[sdata.index[j]] = sdata.data[j];
        }
    }

    memset(D1, 0, cpysize);
    memset(D2, 0, cpysize);
    innDcts[0] = innDcts[1] = 0;
    for(i = 0; i < numb; i++)
    {
        v    = crntClust[i];
        if(this->seed != _non_)
        {
            dst1 = PQMath::l2d(tmpC1, 0, sdata, v, this->ndim);
            dst2 = PQMath::l2d(tmpC2, 0, sdata, v, this->ndim);
        }
        else
        {
            dst1 = rand();
            dst2 = rand();
        }

        if(dst1 < dst2)
        {
            for(j = sdata.col[v]; j < sdata.col[v+1]; j++)
            {
                this->D1[sdata.index[j]] += sdata.data[j];
            }
            tmpn1++;
            nwlabels[i] = -1;
        }
        else
        {
            for(j = sdata.col[v]; j < sdata.col[v+1]; j++)
            {
                this->D2[sdata.index[j]] += sdata.data[j];
            }
            tmpn2++;
            nwlabels[i] = v;
        }
    }
    if(tmpn1 == 0)
    {
        optEg = 0;
        n1 = tmpn2;
        n2 = tmpn1;
        E1 = 0;
        E2 = 0;
        _loc1 = nwlbl*this->ndim;
        memset(this->arrayD+_loc1, 0, cpysize);
        delete [] nwlabels;
        nwlabels = NULL;
        crntClust.clear();
        return false;
    }


    if(tmpn2 == 0)
    {
        optEg = 0;
        n1 = tmpn1;
        n2 = tmpn2;
        E2 = 0;
        E1 = 0;
        _loc1 = nwlbl*this->ndim;
        memset(this->arrayD+_loc1, 0, cpysize);
        delete [] nwlabels;
        nwlabels = NULL;
        crntClust.clear();
        return false;
    }
    innDcts[0] = PQMath::dvec_norm(D1, this->ndim, 2);
    innDcts[1] = PQMath::dvec_norm(D2, this->ndim, 2);
    innDcts[0] = innDcts[0]*innDcts[0];
    innDcts[1] = innDcts[1]*innDcts[1];

    optEg = sqrt(innDcts[0]) +  sqrt(innDcts[1]);
    int *rl = new int[numb];
    /******** Incremental optimization *****/
    do
    {
        for(i = 0; i < numb; i++)
        {
            rl[i] = i;
        }
        random_shuffle(rl, rl + numb);
        _UPDATE_ = false;
        for(iter = 0; iter < numb; iter++)
        {
            mvs[0]++;
            memcpy(tmpD1, D1, cpysize);
            memcpy(tmpD2, D2, cpysize);
            r   = (unsigned int)floor((rand()/(RAND_MAX+1.0f))*numb);
            r   = r>=numb?(numb-1):r;
            //r = rl[iter];
            v   = crntClust[r];
            F2S = false;

            if(nwlabels[r] == -1)
            {
                if(tmpn1 > 1)
                {
                    tmpn1--;
                    tmpn2++;
                    nwlabels[r] = v;
                    F2S = true;
                }
                else
                {
                    continue;
                }
            }
            else
            {
                if(tmpn2 > 1)
                {

                    tmpn1++;
                    tmpn2--;
                    nwlabels[r] = -1;
                    F2S = false;
                }
                else
                {
                    continue;
                }
            }

            tmp1 = tmp2 = 0;
            for(j = sdata.col[v]; j < sdata.col[v+1]; j++)
            {
                tmp1 += tmpD1[sdata.index[j]]*sdata.data[j];
                tmp2 += tmpD2[sdata.index[j]]*sdata.data[j];
            }

            ///tmpEg = (this->*optFunc)(tmpD1, tmpD2, this->ndim, tmpn1, tmpn2);
            len = this->lens[v];
            if(F2S)
            {
                tmp1 = innDcts[0] - 2*tmp1 + len;
                tmp2 = innDcts[1] + 2*tmp2 + len;

            }
            else
            {
                tmp1 = innDcts[0] + 2*tmp1 + len;
                tmp2 = innDcts[1] - 2*tmp2 + len;
            }

            tmpEg = sqrt(tmp1) + sqrt(tmp2);

            if(tmpEg > (optEg+0.00001))
            {
                if(F2S)
                {
                    for(j = sdata.col[v]; j < sdata.col[v+1]; j++)
                    {
                        tmpD1[sdata.index[j]] -= sdata.data[j];
                        tmpD2[sdata.index[j]] += sdata.data[j];
                    }
                }
                else
                {
                    for(j = sdata.col[v]; j < sdata.col[v+1]; j++)
                    {
                        tmpD1[sdata.index[j]] += sdata.data[j];
                        tmpD2[sdata.index[j]] -= sdata.data[j];
                    }
                }
                memcpy(D1, tmpD1, cpysize);
                memcpy(D2, tmpD2, cpysize);
                optEg = tmpEg;
                innDcts[0] = tmp1;
                innDcts[1] = tmp2;
                _UPDATE_ = true;
                mvs[1]++;
            }
            else
            {
                ///roll-back if being updated
                if(nwlabels[r] == -1)
                {
                    tmpn1--;
                    tmpn2++;
                    nwlabels[r] = v;
                }
                else
                {
                    tmpn1++;
                    tmpn2--;
                    nwlabels[r] = -1;
                }
            }
        }
    }
    while(_UPDATE_);
    _loc1 = clabel*this->ndim;
    memcpy(this->arrayD+_loc1, D1, cpysize);
    _loc1 = nwlbl*this->ndim;
    memcpy(this->arrayD+_loc1, D2, cpysize);

    if(optEg > 0)
    {
        n1 = tmpn1;
        n2 = tmpn2;
    }

    for(i = 0; i < numb; i++)
    {
        v = nwlabels[i];
        if(v != -1)
        {
            labels[v] = nwlbl;
        }
    }

    double len1, len2, tmpE1 = 0, tmpE2 = 0;
    len1 = len2 = 0;
    int k = 0;
    for(i = 0; i < numb; i++)
    {
        v = nwlabels[i];
        k = crntClust[i];
        if(v != -1)
        {
            len2 += lens[k];
        }
        else
        {
            len1 += lens[k];
        }
    }

    /**estimate average intra-cluster distance*/
    for(j = 0; j < this->ndim; j++)
    {
        C1[j] = D1[j]/n1;
        C2[j] = D2[j]/n2;
    }
    tmpE1 = PQMath::dvec_norm(C1, ndim, 2);
    tmpE2 = PQMath::dvec_norm(C2, ndim, 2);
    if(n1 > 1)
    {
        E1 = len1 - n1*tmpE1*tmpE1;
        E1 = 2*(E1/(n1-1));
    }
    else
    {
        E1 = 0;
    }

    if(n2 > 1)
    {
        E2 = len2 - n2*tmpE2*tmpE2;
        E2 = 2*(E2/(n2-1));
    }
    else
    {
        E2 = 0;
    }

    crntClust.clear();
    delete [] nwlabels;
    delete [] rl;
    nwlabels = NULL;
    return true;
}


bool XBKMeans::incrOptzI1s(const int clabel, double &E1, double &E2, int &n1, int &n2, const int nwlbl)
{
    unsigned int _loc1 = 0, i = 0, j = 0, r = 0;
    unsigned int cpysize = sizeof(double)*this->ndim;
    double tmp1 = 0, tmp2 = 0, len = 0, x = 0, y = 0, tmpx = 0, tmpy = 0, tmpx1 = 0, tmpy1 = 0;
    double tmp11 = 0, tmp21 = 0;

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
        E2 = E1 = 0.0f;
        n1 = n2 = 1;
        i  = crntClust[1];
        this->labels[i] = nwlbl;
        crntClust.clear();
        return true;
    }

    int *nwlabels = new int[numb];
    memset(nwlabels, 0, sizeof(int)*numb);

    int v = 0, tmpn1 = 0, tmpn2 = 0;
    unsigned int sed1 = 0, sed2 = 0, iter = 0;
    double tmpEg = 0, optEg = 0, dst1 = 0,  dst2 = 0;
    bool _UPDATE_ = true, F2S = false;
    n1 = n2 = 0;
    /**** Initial assigment ****************/
    if(this->seed != _non_)
    {
        (this->*sedFunc)(sed1, sed2, crntClust);
        memset(tmpC1, 0, sizeof(double)*this->ndim);
        memset(tmpC2, 0, sizeof(double)*this->ndim);
        for(j = sdata.col[sed1]; j < sdata.col[sed1 + 1]; j++)
        {
            tmpC1[sdata.index[j]] = sdata.data[j];
        }
        for(j = sdata.col[sed2]; j < sdata.col[sed2 + 1]; j++)
        {
            tmpC2[sdata.index[j]] = sdata.data[j];
        }
    }

    memset(D1, 0, cpysize);
    memset(D2, 0, cpysize);
    innDcts[0] = innDcts[1] = 0;

    for(i = 0; i < numb; i++)
    {
        v    = crntClust[i];
        if(this->seed != _non_)
        {
            dst1 = PQMath::l2d(tmpC1, 0, sdata, v, this->ndim);
            dst2 = PQMath::l2d(tmpC2, 0, sdata, v, this->ndim);
        }
        else
        {
            dst1 = rand();
            dst2 = rand();
        }

        if(dst1 < dst2)
        {
            for(j = sdata.col[v]; j < sdata.col[v+1]; j++)
            {
                this->D1[sdata.index[j]] += sdata.data[j];
            }
            tmpn1++;
            nwlabels[i] = -1;
        }
        else
        {
            for(j = sdata.col[v]; j < sdata.col[v+1]; j++)
            {
                this->D2[sdata.index[j]] += sdata.data[j];
            }
            tmpn2++;
            nwlabels[i] = v;
        }
    }

    if(tmpn1 == 0)
    {
        optEg = 0;
        n1 = tmpn2;
        n2 = tmpn1;
        E1 = 0;
        E2 = 0;
        _loc1 = nwlbl*this->ndim;
        memset(this->arrayD+_loc1, 0, cpysize);
        delete [] nwlabels;
        nwlabels = NULL;
        crntClust.clear();
        return false;
    }


    if(tmpn2 == 0)
    {
        optEg = 0;
        n1 = tmpn1;
        n2 = tmpn2;
        E2 = 0;
        E1 = 0;
        _loc1 = nwlbl*this->ndim;
        memset(this->arrayD+_loc1, 0, cpysize);
        delete [] nwlabels;
        nwlabels = NULL;
        crntClust.clear();
        return false;
    }
    innDcts[0] = PQMath::dvec_norm(D1, this->ndim, 2);
    innDcts[1] = PQMath::dvec_norm(D2, this->ndim, 2);
    innDcts[0] = innDcts[0]*innDcts[0];
    innDcts[1] = innDcts[1]*innDcts[1];
    x = 0;
    y = 0;
    for(i = 0; i < numb; i++)
    {
        if(nwlabels[i] == -1)
            x += this->lens[crntClust[i]];
        else
            y += this->lens[crntClust[i]];
    }

    if(tmpn1 != 1)
    {
        tmpx1 = x*tmpn1/(tmpn1 - 1);
        tmp11 = innDcts[0]/(tmpn1-1);
    }
    else
    {
        tmpx1 = len;
        tmp11 = 0;
    }
    if(tmpn2 != 1)
    {
        tmpy1 = y*tmpn2/(tmpn2 - 1);
        tmp21 = innDcts[1]/(tmpn2-1);
    }
    else
    {
        tmpy1 = len;
        tmp21 = 0;
    }

    optEg = tmp11 + tmp21 - tmpx1 - tmpy1;

    /******** Incremental optimization *****/
    int *rl = new int[numb];

    do
    {

        for(i = 0; i < numb; i++)
        {
            rl[i] = i;
        }
        random_shuffle(rl, rl + numb);
        _UPDATE_ = false;
        for(iter = 0; iter < numb; iter++)
        {
            mvs[0]++;
            memcpy(tmpD1, D1, cpysize);
            memcpy(tmpD2, D2, cpysize);
            r   = rl[iter];
            v   = crntClust[r];
            F2S = false;

            if(nwlabels[r] == -1)
            {
                if(tmpn1 > 1)
                {
                    tmpn1--;
                    tmpn2++;
                    nwlabels[r] = v;
                    F2S = true;
                }
                else
                {
                    continue;
                }
            }
            else
            {
                if(tmpn2 > 1)
                {

                    tmpn1++;
                    tmpn2--;
                    nwlabels[r] = -1;
                    F2S = false;
                }
                else
                {
                    continue;
                }
            }

            tmp1 = tmp2 = 0;
            for(j = sdata.col[v]; j < sdata.col[v+1]; j++)
            {
                tmp1 += tmpD1[sdata.index[j]]*sdata.data[j];
                tmp2 += tmpD2[sdata.index[j]]*sdata.data[j];
            }
            ///tmpEg = (this->*optFunc)(tmpD1, tmpD2, this->ndim, tmpn1, tmpn2);
            len = this->lens[v];
            if(F2S)
            {
                tmp1 = innDcts[0] - 2*tmp1 + len;
                tmp2 = innDcts[1] + 2*tmp2 + len;
                tmpx = x - len;
                tmpy = y + len;
                if(tmpn1 != 1)
                {
                    tmpx1 = tmpx*tmpn1/(tmpn1 - 1);
                    tmp11 = tmp1/(tmpn1-1);
                }
                else
                {
                    tmpx1 = len;
                    tmp11 = 0;
                }
                if(tmpn2 != 1)
                {
                    tmpy1 = tmpy*tmpn2/(tmpn2 - 1);
                    tmp21 = tmp2/(tmpn2 - 1);
                }
                else
                {
                    tmpy1 = len;
                    tmp21 = 0;
                }
            }
            else
            {
                tmp1 = innDcts[0] + 2*tmp1 + len;
                tmp2 = innDcts[1] - 2*tmp2 + len;
                tmpx = x + len;
                tmpy = y - len;
                if(tmpn1 != 1)
                {
                    tmpx1 = tmpx*tmpn1/(tmpn1 - 1);
                    tmp11 = tmp1/(tmpn1-1);
                }
                else
                {
                    tmpx1 = len;
                    tmp11 = 0;
                }
                if(tmpn2 != 1)
                {
                    tmpy1 = tmpy*tmpn2/(tmpn2 - 1);
                    tmp21 = tmp2/(tmpn2 - 1);
                }
                else
                {
                    tmpy1 = len;
                    tmp21 = 0;
                }
            }

            tmpEg = tmp11 + tmp21 - tmpx1 - tmpy1;

            if(tmpEg > (optEg+0.00001))
            {
                if(F2S)
                {
                    for(j = sdata.col[v]; j < sdata.col[v+1]; j++)
                    {
                        tmpD1[sdata.index[j]] -= sdata.data[j];
                        tmpD2[sdata.index[j]] += sdata.data[j];
                    }
                }
                else
                {
                    for(j = sdata.col[v]; j < sdata.col[v+1]; j++)
                    {
                        tmpD1[sdata.index[j]] += sdata.data[j];
                        tmpD2[sdata.index[j]] -= sdata.data[j];
                    }
                }
                memcpy(D1, tmpD1, cpysize);
                memcpy(D2, tmpD2, cpysize);
                optEg = tmpEg;
                x = tmpx;
                y = tmpy;
                innDcts[0] = tmp1;
                innDcts[1] = tmp2;
                _UPDATE_ = true;
                mvs[1]++;
            }
            else
            {
                ///roll-back if being updated
                if(nwlabels[r] == -1)
                {
                    tmpn1--;
                    tmpn2++;
                    nwlabels[r] = v;
                }
                else
                {
                    tmpn1++;
                    tmpn2--;
                    nwlabels[r] = -1;
                }
            }
        }
    }
    while(_UPDATE_);
    _loc1 = clabel*this->ndim;
    memcpy(this->arrayD+_loc1, D1, cpysize);
    _loc1 = nwlbl*this->ndim;
    memcpy(this->arrayD+_loc1, D2, cpysize);

    // if(optEg > 0)
    {
        n1 = tmpn1;
        n2 = tmpn2;
    }

    for(i = 0; i < numb; i++)
    {
        v = nwlabels[i];
        if(v != -1)
        {
            labels[v] = nwlbl;
        }
    }

    double len1, len2, tmpE1 = 0, tmpE2 = 0;
    len1 = len2 = 0;
    int k = 0;
    for(i = 0; i < numb; i++)
    {
        v = nwlabels[i];
        k = crntClust[i];
        if(v != -1)
        {
            len2 += lens[k];
        }
        else
        {
            len1 += lens[k];
        }
    }

    /**estimate average intra-cluster distance*/
    for(j = 0; j < this->ndim; j++)
    {
        C1[j] = D1[j]/n1;
        C2[j] = D2[j]/n2;
    }
    tmpE1 = PQMath::dvec_norm(C1, ndim, 2);
    tmpE2 = PQMath::dvec_norm(C2, ndim, 2);
    if(n1 > 1)
    {
        E1 = len1 - n1*tmpE1*tmpE1;
        E1 = 2*(E1/(n1-1));
    }
    else
    {
        E1 = 0;
    }

    if(n2 > 1)
    {
        E2 = len2 - n2*tmpE2*tmpE2;
        E2 = 2*(E2/(n2-1));
    }
    else
    {
        E2 = 0;
    }

    crntClust.clear();
    delete [] nwlabels;
    delete [] rl;
    nwlabels = NULL;
    return true;
}

bool XBKMeans::incrOptzE1s(const int clabel, double &E1, double &E2, int &n1, int &n2, const int nwlbl)
{
    unsigned int _loc1 = 0, i = 0, j = 0, r = 0;
    unsigned int cpysize = sizeof(double)*this->ndim;
    double tmp1 = 0, tmp2 = 0, len = 0, tmpij = 0, x = 0, y = 0, tmpx = 0, tmpy = 0;
    double innerij = 0;

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
        E2 = E1 = 0.0f;
        n1 = n2 = 1;
        i  = crntClust[1];
        this->labels[i] = nwlbl;
        crntClust.clear();
        return true;
    }

    int *nwlabels = new int[numb];
    memset(nwlabels, 0, sizeof(int)*numb);

    int v = 0, tmpn1 = 0, tmpn2 = 0;
    unsigned int iter = 0, sed1 = 0, sed2 = 0;
    double tmpEg = 0, optEg = 0, dst1 = 0,  dst2 = 0;
    bool _UPDATE_ = true, F2S = false;
    n1 = n2 = 0;
    /**** Initial assigment ****************/
    if(this->seed != _non_)
    {
        (this->*sedFunc)(sed1, sed2, crntClust);
        memset(tmpC1, 0, sizeof(double)*this->ndim);
        memset(tmpC2, 0, sizeof(double)*this->ndim);
        for(j = sdata.col[sed1]; j < sdata.col[sed1 + 1]; j++)
        {
            tmpC1[sdata.index[j]] = sdata.data[j];
        }
        for(j = sdata.col[sed2]; j < sdata.col[sed2 + 1]; j++)
        {
            tmpC2[sdata.index[j]] = sdata.data[j];
        }
    }

    memset(D1, 0, cpysize);
    memset(D2, 0, cpysize);
    innDcts[0] = innDcts[1] = 0;

    for(i = 0; i < numb; i++)
    {
        v    = crntClust[i];
        if(this->seed != _non_)
        {
            dst1 = PQMath::l2d(tmpC1, 0, sdata, v, this->ndim);
            dst2 = PQMath::l2d(tmpC2, 0, sdata, v, this->ndim);
        }
        else
        {
            dst1 = rand();
            dst2 = rand();
        }

        if(dst1 < dst2)
        {
            for(j = sdata.col[v]; j < sdata.col[v+1]; j++)
            {
                this->D1[sdata.index[j]] += sdata.data[j];
            }
            tmpn1++;
            nwlabels[i] = -1;
        }
        else
        {
            for(j = sdata.col[v]; j < sdata.col[v+1]; j++)
            {
                this->D2[sdata.index[j]] += sdata.data[j];
            }
            tmpn2++;
            nwlabels[i] = v;
        }
    }

    if(tmpn1 == 0)
    {
        optEg = 0;
        n1 = tmpn2;
        n2 = tmpn1;
        E1 = 0;
        E2 = 0;
        _loc1 = nwlbl*this->ndim;
        memset(this->arrayD+_loc1, 0, cpysize);
        delete [] nwlabels;
        nwlabels = NULL;
        crntClust.clear();
        return false;
    }


    if(tmpn2 == 0)
    {
        optEg = 0;
        n1 = tmpn1;
        n2 = tmpn2;
        E2 = 0;
        E1 = 0;
        _loc1 = nwlbl*this->ndim;
        memset(this->arrayD+_loc1, 0, cpysize);
        delete [] nwlabels;
        nwlabels = NULL;
        crntClust.clear();
        return false;
    }
    innDcts[0] = PQMath::dvec_norm(D1, this->ndim, 2);
    innDcts[1] = PQMath::dvec_norm(D2, this->ndim, 2);
    innDcts[0] = innDcts[0]*innDcts[0];
    innDcts[1] = innDcts[1]*innDcts[1];
    innerij = PQMath::innerProduct(D1, D2, this->ndim);
    x = 0;
    y = 0;
    for(i = 0; i < numb; i++)
    {
        if(nwlabels[i] == -1)
            x += this->lens[crntClust[i]];
        else
            y += this->lens[crntClust[i]];
    }

    optEg = (tmpn1*x + tmpn2*y - 2*innerij)*tmpn1*tmpn2;

    /******** Incremental optimization *****/

    int *rl = new int[numb];
    do
    {
        for(i = 0; i < numb; i++)
        {
            rl[i] = i;
        }
        random_shuffle(rl, rl + numb);
        _UPDATE_ = false;
        for(iter = 0; iter < numb; iter++)
        {
            mvs[0]++;
            memcpy(tmpD1, D1, cpysize);
            memcpy(tmpD2, D2, cpysize);
            r   = rl[iter];
            v   = crntClust[r];
            F2S = false;

            if(nwlabels[r] == -1)
            {
                if(tmpn1 > 1)
                {
                    tmpn1--;
                    tmpn2++;
                    nwlabels[r] = v;
                    F2S = true;
                }
                else
                {
                    continue;
                }
            }
            else
            {
                if(tmpn2 > 1)
                {

                    tmpn1++;
                    tmpn2--;
                    nwlabels[r] = -1;
                    F2S = false;
                }
                else
                {
                    continue;
                }
            }

            tmp1 = tmp2 = 0;
            for(j = sdata.col[v]; j < sdata.col[v+1]; j++)
            {
                tmp1 += tmpD1[sdata.index[j]]*sdata.data[j];
                tmp2 += tmpD2[sdata.index[j]]*sdata.data[j];
            }
            ///tmpEg = (this->*optFunc)(tmpD1, tmpD2, this->ndim, tmpn1, tmpn2);
            len = this->lens[v];
            if(F2S)
            {
                tmpij = innerij + tmp1 - tmp2 - len;
                tmpx = x - len;
                tmpy = y + len;
            }
            else
            {
                tmpij = innerij - tmp1 + tmp2 - len;
                tmpx = x + len;
                tmpy = y - len;
            }

            tmpEg = (tmpn1*tmpx + tmpn2*tmpy - 2*tmpij)*tmpn1*tmpn2;
            // if(tmpEg < 0)
            // cout<<tmpEg<<"\n";

            if(tmpEg > (optEg+0.00001))
            {
                if(F2S)
                {
                    for(j = sdata.col[v]; j < sdata.col[v+1]; j++)
                    {
                        tmpD1[sdata.index[j]] -= sdata.data[j];
                        tmpD2[sdata.index[j]] += sdata.data[j];
                    }
                }
                else
                {
                    for(j = sdata.col[v]; j < sdata.col[v+1]; j++)
                    {
                        tmpD1[sdata.index[j]] += sdata.data[j];
                        tmpD2[sdata.index[j]] -= sdata.data[j];
                    }
                }
                memcpy(D1, tmpD1, cpysize);
                memcpy(D2, tmpD2, cpysize);
                innerij  = tmpij;
                optEg = tmpEg;
                x = tmpx;
                y = tmpy;
                _UPDATE_ = true;
                mvs[1]++;
            }
            else
            {
                ///roll-back if being updated
                if(nwlabels[r] == -1)
                {
                    tmpn1--;
                    tmpn2++;
                    nwlabels[r] = v;
                }
                else
                {
                    tmpn1++;
                    tmpn2--;
                    nwlabels[r] = -1;
                }
            }
        }
    }
    while(_UPDATE_);
    _loc1 = clabel*this->ndim;
    memcpy(this->arrayD+_loc1, D1, cpysize);
    _loc1 = nwlbl*this->ndim;
    memcpy(this->arrayD+_loc1, D2, cpysize);

    // if(optEg > 0)
    {
        n1 = tmpn1;
        n2 = tmpn2;
    }

    for(i = 0; i < numb; i++)
    {
        v = nwlabels[i];
        if(v != -1)
        {
            labels[v] = nwlbl;
        }
    }

    double len1, len2, tmpE1 = 0, tmpE2 = 0;
    len1 = len2 = 0;
    int k = 0;
    for(i = 0; i < numb; i++)
    {
        v = nwlabels[i];
        k = crntClust[i];
        if(v != -1)
        {
            len2 += lens[k];
        }
        else
        {
            len1 += lens[k];
        }
    }

    /**estimate average intra-cluster distance*/
    for(j = 0; j < this->ndim; j++)
    {
        C1[j] = D1[j]/n1;
        C2[j] = D2[j]/n2;
    }
    tmpE1 = PQMath::dvec_norm(C1, ndim, 2);
    tmpE2 = PQMath::dvec_norm(C2, ndim, 2);
    if(n1 > 1)
    {
        E1 = len1 - n1*tmpE1*tmpE1;
        E1 = 2*(E1/(n1-1));
    }
    else
    {
        E1 = 0;
    }

    if(n2 > 1)
    {
        E2 = len2 - n2*tmpE2*tmpE2;
        E2 = 2*(E2/(n2-1));
    }
    else
    {
        E2 = 0;
    }

    crntClust.clear();
    delete [] nwlabels;
    delete [] rl;
    nwlabels = NULL;
    rl = NULL;
    return true;
}


bool XBKMeans::incrOptzE2s(const int clabel, double &E1, double &E2, int &n1, int &n2, const int nwlbl)
{
    unsigned int _loc1 = 0, i = 0, j = 0, r = 0;
    unsigned int cpysize = sizeof(double)*this->ndim;
    double tmp1 = 0, tmp2 = 0, len = 0, tmpij = 0;
    double innerij ;

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
        E2 = E1 = 0.0f;
        n1 = n2 = 1;
        i  = crntClust[1];
        this->labels[i] = nwlbl;
        crntClust.clear();
        return true;
    }

    int *nwlabels = new int[numb];
    memset(nwlabels, 0, sizeof(int)*numb);

    int v = 0, tmpn1 = 0, tmpn2 = 0;
    unsigned int sed1 = 0, sed2 = 0, iter = 0;
    double tmpEg = 0, optEg = 0, dst1 = 0,  dst2 = 0;
    bool _UPDATE_ = true, F2S = false;
    n1 = n2 = 0;
    /**** Initial assigment ****************/
    if(this->seed != _non_)
    {
        (this->*sedFunc)(sed1, sed2, crntClust);
        memset(tmpC1, 0, sizeof(double)*this->ndim);
        memset(tmpC2, 0, sizeof(double)*this->ndim);
        for(j = sdata.col[sed1]; j < sdata.col[sed1 + 1]; j++)
        {
            tmpC1[sdata.index[j]] = sdata.data[j];
        }
        for(j = sdata.col[sed2]; j < sdata.col[sed2 + 1]; j++)
        {
            tmpC2[sdata.index[j]] = sdata.data[j];
        }
    }

    memset(D1, 0, cpysize);
    memset(D2, 0, cpysize);
    innDcts[0] = innDcts[1] = 0;

    for(i = 0; i < numb; i++)
    {
        v    = crntClust[i];
        if(this->seed != _non_)
        {
            dst1 = PQMath::l2d(tmpC1, 0, sdata, v, this->ndim);
            dst2 = PQMath::l2d(tmpC2, 0, sdata, v, this->ndim);
        }
        else
        {
            dst1 = rand();
            dst2 = rand();
        }

        if(dst1 < dst2)
        {
            for(j = sdata.col[v]; j < sdata.col[v+1]; j++)
            {
                this->D1[sdata.index[j]] += sdata.data[j];
            }
            tmpn1++;
            nwlabels[i] = -1;
        }
        else
        {
            for(j = sdata.col[v]; j < sdata.col[v+1]; j++)
            {
                this->D2[sdata.index[j]] += sdata.data[j];
            }
            tmpn2++;
            nwlabels[i] = v;
        }
    }

    if(tmpn1 == 0)
    {
        optEg = 0;
        n1 = tmpn2;
        n2 = tmpn1;
        E1 = 0;
        E2 = 0;
        _loc1 = nwlbl*this->ndim;
        memset(this->arrayD+_loc1, 0, cpysize);
        delete [] nwlabels;
        nwlabels = NULL;
        crntClust.clear();
        return false;
    }


    if(tmpn2 == 0)
    {
        optEg = 0;
        n1 = tmpn1;
        n2 = tmpn2;
        E2 = 0;
        E1 = 0;
        _loc1 = nwlbl*this->ndim;
        memset(this->arrayD+_loc1, 0, cpysize);
        delete [] nwlabels;
        nwlabels = NULL;
        crntClust.clear();
        return false;
    }
    innDcts[0] = PQMath::dvec_norm(D1, this->ndim, 2);
    innDcts[1] = PQMath::dvec_norm(D2, this->ndim, 2);
    innDcts[0] = innDcts[0]*innDcts[0];
    innDcts[1] = innDcts[1]*innDcts[1];
    innerij = PQMath::innerProduct(D1, D2, this->ndim);

    optEg = innDcts[0]/tmpn1*tmpn2 + innDcts[1]/tmpn2*tmpn1 - 2*innerij;


    /******** Incremental optimization *****/

    int *rl = new int[numb];
    do
    {
        for(i = 0; i < numb; i++)
        {
            rl[i] = i;
        }
        random_shuffle(rl, rl + numb);
        _UPDATE_ = false;
        for(iter = 0; iter < numb; iter++)
        {
            mvs[0]++;
            memcpy(tmpD1, D1, cpysize);
            memcpy(tmpD2, D2, cpysize);
            r   = (unsigned int)floor((rand()/(RAND_MAX+1.0f))*numb);
            r   = r>=numb?(numb-1):r;
            //r = rl[iter];
            v   = crntClust[r];
            F2S = false;

            if(nwlabels[r] == -1)
            {
                if(tmpn1 > 1)
                {
                    tmpn1--;
                    tmpn2++;
                    nwlabels[r] = v;
                    F2S = true;
                }
                else
                {
                    continue;
                }
            }
            else
            {
                if(tmpn2 > 1)
                {

                    tmpn1++;
                    tmpn2--;
                    nwlabels[r] = -1;
                    F2S = false;
                }
                else
                {
                    continue;
                }
            }

            tmp1 = tmp2 = 0;
            for(j = sdata.col[v]; j < sdata.col[v+1]; j++)
            {
                tmp1 += tmpD1[sdata.index[j]]*sdata.data[j];
                tmp2 += tmpD2[sdata.index[j]]*sdata.data[j];
            }

            ///tmpEg = (this->*optFunc)(tmpD1, tmpD2, this->ndim, tmpn1, tmpn2);
            len = this->lens[v];
            if(F2S)
            {
                tmpij = innerij + tmp1 - tmp2 - len;
                tmp1 = innDcts[0] - 2*tmp1 + len;
                tmp2 = innDcts[1] + 2*tmp2 + len;
            }
            else
            {
                tmpij = innerij - tmp1 + tmp2 - len;
                tmp1 = innDcts[0] + 2*tmp1 + len;
                tmp2 = innDcts[1] - 2*tmp2 + len;

            }

            tmpEg = tmp1/tmpn1*tmpn2 + tmp2/tmpn2*tmpn1 - 2*tmpij;

            if(tmpEg > (optEg+0.00001))
            {
                if(F2S)
                {
                    for(j = sdata.col[v]; j < sdata.col[v+1]; j++)
                    {
                        tmpD1[sdata.index[j]] -= sdata.data[j];
                        tmpD2[sdata.index[j]] += sdata.data[j];
                    }
                }
                else
                {
                    for(j = sdata.col[v]; j < sdata.col[v+1]; j++)
                    {
                        tmpD1[sdata.index[j]] += sdata.data[j];
                        tmpD2[sdata.index[j]] -= sdata.data[j];
                    }
                }
                memcpy(D1, tmpD1, cpysize);
                memcpy(D2, tmpD2, cpysize);
                innerij  = tmpij;
                optEg = tmpEg;
                innDcts[0] = tmp1;
                innDcts[1] = tmp2;
                _UPDATE_ = true;
                mvs[1]++;

            }
            else
            {
                ///roll-back if being updated
                if(nwlabels[r] == -1)
                {
                    tmpn1--;
                    tmpn2++;
                    nwlabels[r] = v;
                }
                else
                {
                    tmpn1++;
                    tmpn2--;
                    nwlabels[r] = -1;
                }
            }
        }
    }
    while(_UPDATE_);
    _loc1 = clabel*this->ndim;
    memcpy(this->arrayD+_loc1, D1, cpysize);
    _loc1 = nwlbl*this->ndim;
    memcpy(this->arrayD+_loc1, D2, cpysize);

    if(optEg > 0)
    {
        n1 = tmpn1;
        n2 = tmpn2;
    }

    for(i = 0; i < numb; i++)
    {
        v = nwlabels[i];
        if(v != -1)
        {
            labels[v] = nwlbl;
        }
    }

    double len1, len2, tmpE1 = 0, tmpE2 = 0;
    len1 = len2 = 0;
    int k = 0;
    for(i = 0; i < numb; i++)
    {
        v = nwlabels[i];
        k = crntClust[i];
        if(v != -1)
        {
            len2 += lens[k];
        }
        else
        {
            len1 += lens[k];
        }
    }

    /**estimate average intra-cluster distance*/
    for(j = 0; j < this->ndim; j++)
    {
        C1[j] = D1[j]/n1;
        C2[j] = D2[j]/n2;
    }
    tmpE1 = PQMath::dvec_norm(C1, ndim, 2);
    tmpE2 = PQMath::dvec_norm(C2, ndim, 2);
    if(n1 > 1)
    {
        E1 = len1 - n1*tmpE1*tmpE1;
        E1 = 2*(E1/(n1-1));
    }
    else
    {
        E1 = 0;
    }

    if(n2 > 1)
    {
        E2 = len2 - n2*tmpE2*tmpE2;
        E2 = 2*(E2/(n2-1));
    }
    else
    {
        E2 = 0;
    }

    crntClust.clear();
    delete [] nwlabels;
    delete [] rl;
    nwlabels = NULL;
    return true;
}


bool XBKMeans::tkmOptz(const int clabel, double &E1, double &E2, int &n1, int &n2, const int nwlbl)
{
    /**static unsigned int step = 0;
    char tmpfile[1024];/**/

    unsigned int i = 0, j = 0, k = 0, loc = 0;
    unsigned int _loc1 = 0, _loc2 = 0;
    for(i = 0; i < count; i++)
    {
        if(clabel == labels[i])
        {
            crntClust.push_back(i);
        }
    }
    unsigned int numb = crntClust.size();
    if(numb == 2)
    {
        E2 = E1 = 0.0f;
        n1 = n2 = 1;
        i  = crntClust[1];
        this->labels[i] = nwlbl;
        crntClust.clear();
        return true;
    }
    nwlabels   = new int[numb];
    memset(nwlabels, 0, sizeof(int)*numb);

    double _sumDst = 0, dst1 = 0, dst2 = 0, sumDst = 0, _ERR = 0;
    double tmpE1 = 0, tmpE2 = 0, len1 = 0, len2 = 0;

    int v = 0, iter = 0, tmpn1 = 0, tmpn2 = 0;
    unsigned int sed1 = 0, sed2 = 0;
    n1 = n2 = 0;
    unsigned int cpysize = sizeof(double)*this->ndim;

    (this->*sedFunc)(sed1, sed2, crntClust);
    tmpn1 = tmpn2 = 0;
    tmpE1 = tmpE2 = 0;
    _loc1 = sed1*this->ndim;
    _loc2 = sed2*this->ndim;

    for(j = 0; j < this->ndim; j++)
    {
        tmpC1[j] = data[_loc1+j];
        tmpC2[j] = data[_loc2+j];
    }

    memset(tmpD1, 0, cpysize);
    memset(tmpD2, 0, cpysize);


    for(i = 0; i < numb; i++)
    {
        v = crntClust[i];
        dst1 = PQMath::l2d(tmpC1, 0, data, v, this->ndim);
        dst2 = PQMath::l2d(tmpC2, 0, data, v, this->ndim);
        loc = v*this->ndim;
        if(dst1 < dst2)
        {
            for(j = 0; j < this->ndim; j++)
            {
                tmpD1[j] += data[loc+j];
            }
            tmpn1++;
            nwlabels[i] = -1;
            tmpE1 += dst1;
        }
        else
        {
            for(j = 0; j < this->ndim; j++)
            {
                tmpD2[j] += data[loc+j];
            }
            tmpn2++;
            nwlabels[i] = v;
            tmpE2 += dst2;
        }
    }

    sumDst = tmpE1 + tmpE2;
    for(j = 0; j < this->ndim; j++)
    {
        tmpC1[j] = tmpD1[j]/tmpn1;
        tmpC2[j] = tmpD2[j]/tmpn2;
    }

    iter = 0;
    do
    {
        memset(tmpD1, 0, cpysize);
        memset(tmpD2, 0, cpysize);

        tmpn1 = tmpn2 = 0;
        tmpE1 = tmpE2 = 0;

        for(i = 0; i < numb; i++)
        {
            v = crntClust[i];
            dst1 = PQMath::l2d(tmpC1, 0, data, v, this->ndim);
            dst2 = PQMath::l2d(tmpC2, 0, data, v, this->ndim);
            loc  = v*this->ndim;
            mvs[0]++;
            if(dst1 > dst2)
            {
                for(j = 0; j < this->ndim; j++)
                {
                    tmpD1[j] += data[loc+j];
                }
                tmpn1++;
                nwlabels[i] = -1;
                tmpE1 += dst1;
            }
            else
            {
                for(j = 0; j < this->ndim; j++)
                {
                    tmpD2[j] += data[loc+j];
                }
                tmpn2++;
                nwlabels[i]  = v;
                tmpE2 += dst2;
            }

        }

        for(j = 0; j < this->ndim; j++)
        {
            tmpC1[j] = tmpD1[j]/tmpn1;
            tmpC2[j] = tmpD2[j]/tmpn2;
        }
        iter++;
        sumDst  = tmpE1 + tmpE2;
        _ERR    = fabs(sumDst - _sumDst);
        _sumDst = sumDst;
    }
    while(_ERR > XBKMeans::Err0);

    _loc1 = clabel*this->ndim;
    _loc2 = nwlbl*this->ndim;
    E1 = tmpE1;
    E2 = tmpE2;
    n1 = tmpn1;
    n2 = tmpn2;

    if(n1 == 0 || n2 == 0)
    {
        crntClust.clear();
        delete [] nwlabels;
        nwlabels = NULL;
        E1 = 0;
        E2 = 0;
        return false;
    }

    memcpy(arrayD + _loc1, tmpD1, cpysize);
    memcpy(arrayD + _loc2, tmpD2, cpysize);

    for(i = 0; i < numb; i++)
    {
        labels[crntClust[i]] = clabel;
    }
    for(i = 0; i < numb; i++)
    {
        v = nwlabels[i];
        if(v != -1)
        {
            labels[v] = nwlbl;
        }
    }

    len1 = len2 = 0;
    for(i = 0; i < numb; i++)
    {
        v = nwlabels[i];
        k = crntClust[i];
        if(v != -1)
        {
            len2 += lens[k];
        }
        else
        {
            len1 += lens[k];
        }
    }

    delete [] nwlabels;
    nwlabels = NULL;

    memcpy(C1, tmpD1, cpysize);
    memcpy(C2, tmpD2, cpysize);

    /**estimate average intra-cluster distance*/
    for(j = 0; j < this->ndim; j++)
    {
        C1[j] = C1[j]/(n1+0.0f);
        C2[j] = C2[j]/(n2+0.0f);
    }
    tmpE1 = PQMath::dvec_norm(C1, ndim, 2);
    tmpE2 = PQMath::dvec_norm(C2, ndim, 2);
    if(n1 > 1)
    {
        E1 = len1 - n1*tmpE1*tmpE1;
        E1 = 2*(E1/(n1-1));
    }
    else
    {
        E1 = 0;
    }
    if(n2 > 1)
    {
        E2 = len2 - n2*tmpE2*tmpE2;
        E2 = 2*(E2/(n2-1));
    }
    else
    {
        E2 = 0;
    }
    crntClust.clear();

    return true;
}

double XBKMeans::I1Func(const double *D1, const double *D2,
                        const unsigned int dim,   const unsigned int tmpn1,
                        const unsigned int tmpn2)
{
    double sumDst = 0, dst1 = 0, dst2 = 0;
    unsigned int i = 0;

    if(tmpn1 == 1)
    {
        dst1 = 0;
    }
    else
    {
        for(i = 0; i < dim; i++)
        {
            dst1 += D1[i]*D1[i];
        }
        dst1 = dst1/(tmpn1-1);
    }

    if(tmpn2 == 1)
    {
        dst2 = 0;
    }
    else
    {
        for(i = 0; i < dim; i++)
        {
            dst2 += D2[i]*D2[i];
        }
        dst2 = dst2/(tmpn2-1);
    }

    sumDst = dst1 + dst2;
    return sumDst;
}


double XBKMeans::I2Func(const double *D1, const double *D2,
                        const unsigned int dim, const unsigned int tmpn1,
                        const unsigned int tmpn2)
{
    double sumDst = 0, dst1 = 0, dst2 = 0;
    unsigned int i = 0;

    for(i = 0; i < dim; i++)
    {
        dst1 += D1[i]*D1[i];
        dst2 += D2[i]*D2[i];
    }
    sumDst = dst1/tmpn1 + dst2/tmpn2;
    return sumDst;
}

double XBKMeans::I3Func(const double *D1, const double *D2, const unsigned int dim,
                        const unsigned int tmpn1, const unsigned int tmpn2)
{
    double sumDst = 0, dst1 = 0, dst2 = 0;
    unsigned int i = 0, numb = 0, k;
    int v;
    double eg1 = 0, eg2 = 0;

    numb = tmpn1 + tmpn2;
    for(i = 0; i < numb; i++)
    {
        v = nwlabels[i];
        k = crntClust[i];
        if(v != -1)
        {
            eg2 += lens[k];
        }
        else
        {
            eg1 += lens[k];
        }
    }


    if(tmpn1 == 1)
    {
        dst1 = 0;
        eg1  = 0;
    }
    else
    {
        for(i = 0; i < dim; i++)
        {
            dst1 += D1[i]*D1[i];
        }
        dst1 = dst1/(tmpn1-1);
        eg1  = eg1/(tmpn1-1);
    }

    if(tmpn2 == 1)
    {
        dst2 = 0;
        eg2  = 0;
    }
    else
    {
        for(i = 0; i < dim; i++)
        {
            dst2 += D2[i]*D2[i];
        }
        dst2 = dst2/(tmpn2-1);
        eg2  = eg2/(tmpn2-1);
    }

    sumDst = dst1 + dst2 - eg1 - eg2;
    return sumDst;
}

double XBKMeans::E1Func(const double *D1, const double *D2,
                        const unsigned int dim, const unsigned int tmpn1,
                        const unsigned int tmpn2)
{
    double sumDst = 0, dst1 = 0, dst2 = 0;
    unsigned int i = 0;
    for(i = 0; i < dim; i++)
    {
        sumDst += D1[i]*D2[i];
    }

    sumDst = dst1/tmpn1 + dst2/tmpn2 - 2*sumDst;
    return sumDst;
}

double XBKMeans::E2Func(const double *D1, const double *D2,
                        const unsigned int dim,   const unsigned int tmpn1,
                        const unsigned int tmpn2)
{
    double sumDst = 0, dx = 0;
    unsigned int i = 0;
    for(i = 0; i < dim; i++)
    {
        dx = D1[i]/tmpn1 - D2[i]/tmpn2;
        sumDst += dx*dx;
    }
    //sumDst = sqrt(sumDst);
    return sumDst;
}

bool XBKMeans::kppTwin(unsigned int &r1, unsigned int &r2, vector<int> &vect)
{
    if(dataType)
        assert(sdata.data);
    else
        assert(data);

    long i = 0, j = 0;
    double sum = 0.0f;
    const int bound = vect.size();
    long c1 = 0, c2 = 0;
    int rseeds[2], sel;

    for(i = 0; i < 2; i++)
    {
        rseeds[i] = i;
    }

    double *disbest = new double[bound];
    double  *distmp = new double[bound];
    for(i = 0; i < bound; i++)
    {
        disbest[i] = RAND_MAX+0.0f;
    }
    double tmp = 0;

    rseeds[0] = (int)floor((rand()/(RAND_MAX+1.0))*bound);
    rseeds[0] = rseeds[0]>=bound?(bound-1):rseeds[0];

    for (i = 1 ; i < 2; i++)
    {
        sel = rseeds[i - 1];
        c1  = vect[sel];
        for (j = 0 ; j < bound; j++)
        {
            c2  = vect[j];
            if(!dataType)
                tmp = PQMath::l2f(this->data, c2, this->data, c1, ndim);
            else
                tmp = PQMath::l2f(this->sdata, c2, this->sdata, c1, ndim);
            if(tmp < disbest[j])
            {
                disbest[j] = tmp;
            }
        }
        memcpy (distmp, disbest, bound * sizeof(double));

        sum = PQMath::dvec_norm(distmp, bound, 1);
        PQMath::dvec_scale(distmp, bound, 1.0/sum);
        double rd = rand()/(RAND_MAX+1.0f);

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
            if(vect[j] != (int)r1)
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


bool XBKMeans::rndTwin(unsigned int &r1, unsigned int &r2, vector<int> &vect)
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


void XBKMeans::saveCenters(const char *dstfn, bool append)
{
    unsigned int clabel = 0, j, loc, rCNum = 0;

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

    ofstream *outStrm  = NULL;
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
        if(this->infoMap[clabel].n > 0)
        {
            for(j = 0; j < this->ndim; j++)
            {
                (*outStrm)<<this->arrayD[loc+j]/this->infoMap[clabel].n<<" ";
            }
            (*outStrm)<<endl;
        }
    }

    outStrm->close();
    cout<<"done\n";
    return ;
}


int XBKMeans::fetchCenters(float *centers)
{
    unsigned int clabel = 0, j = 0, loc = 0, idxi = 0, rCNum = 0;
    assert(centers);
    memset(centers, 0, this->clnumb*this->ndim*sizeof(float));

    for(clabel = 0; clabel < this->clnumb; clabel++)
    {
        if(this->infoMap[clabel].n > 0)
        {
            rCNum++;
        }
    }

    if(!this->_INIT_||rCNum == 0)
    {
        return 0;
    }

    for(clabel = 0; clabel < this->clnumb; clabel++)
    {
        loc   = clabel*this->ndim;

        if(this->infoMap[clabel].n > 0)
        {
            for(j = 0; j < this->ndim; j++)
            {
                centers[idxi + j] = (float)(this->arrayD[loc+j]/this->infoMap[clabel].n);
            }
            rCNum++;
        }
        else
        {
            cout<<"Empty cluster:\t"<<clabel<<"\t"<<this->infoMap[clabel].n<<endl;
        }

        idxi += this->ndim;
    }
    return rCNum;
}

void XBKMeans::hyeTest()
{
   const char *srcMatFn1 = "/home/hye/datasets/clust/sift100k/sift_learn.txt";
   const char *srcMatFn2 = "/home/wlzhao/datasets/bignn/mat/sift1m/sift_base.txt";
   const char *srcMatFn3 = "itm_vlad_hesaff_flickr10m.itm";
   const char *dstFn = "/home/hye/result/xbk1_rslt1_100k.txt";
   int i = 1024;

    XBKMeans *mykm = new XBKMeans();
    mykm->buildcluster(srcMatFn1, dstFn, "non", "large", "i2", 1000, true);
    delete mykm;
}
void XBKMeans::test()
{
   const char *srcMatFn1 = "/home/wlzhao/datasets/bignn/sift1m/sift_learn.txt";
   const char *srcMatFn2 = "/home/wlzhao/datasets/bignn/sift1m/sift_base.txt";
   const char *srcMatFn3 = "itm_vlad_hesaff_flickr10m.itm";
   const char *dstFn = "/home/wlzhao/datasets/clust/xbk1_rslt1_xx.txt";
   int i = 1024;
   //do
   //{
     XBKMeans *mykm = new XBKMeans();
     mykm->buildcluster(srcMatFn1, dstFn, "non", "large", "i2", 1024, true);
     delete mykm;
     //i = i*2;
   //}while(i <= 8192);

}


