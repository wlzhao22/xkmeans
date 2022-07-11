#include "xtkmeans.h"

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


const int XTKMeans::NTRAILS = 128;
const int XTKMeans::NRef    = 1;
const int XTKMeans::Round   = 40;
const int XTKMeans::Pubnum  = 100;


XTKMeans::XTKMeans()
{
    Ds = Cs = tmpCs   = NULL;
    kmMtd   = _xtkmn_;
    strcpy(mthStr, "_xtk_");
    this->_REFER_ = false;
    cout<<"Method ........................... XT K-Means\n";
}

XTKMeans::~XTKMeans()
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
}

bool  XTKMeans::init(const char *srcfn)
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
    ///this->count = 100000;

    cout<<this->count<<"x"<<this->ndim<<endl;
#ifdef _WATCH_
    for(int i = 0; i < 100000; i++)
    {
        for(int j = 0; j < 50; j++)
        {
            this->hits[i][j] = 0;
        }
    }
#endif // _WATCH_

    if(this->data == NULL && this->sdata.data == NULL)
    {
        cout<<"Exceptions ....................... Loading matrix failed!\n";
        exit(0);
    }

    this->_REFER_= false;
    return true;
}

bool XTKMeans::allocMem(const unsigned int clustNum)
{
    this->Ds     = new double[clustNum*this->ndim];
    this->Cs     = new double[clustNum*this->ndim];
    this->tmpCs  = new double[clustNum*this->ndim];

#ifdef _WATCH_
    this->hists  = new int[Round];
    memset(this->hists, 0, sizeof(int)*Round);
#endif

    return true;
}

bool XTKMeans::deallocMem()
{
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

#ifdef _WATCH_
    if(this->hists != NULL)
    {
        delete [] this->hists;
        this->hists = NULL;
    }
#endif

    return true;
}

bool  XTKMeans::init(float *mat, const int row, const int dim)
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

bool XTKMeans::refresh()
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

    deallocMem();

    return true;
}

bool  XTKMeans::config(const char *_seed_, const char *lg_first, const char *crtrn, int verbose)
{

    if(verbose)
        cout<<"Bisecting stradegy ............... ";
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
    {
        if(this->LG_FIRST)
            cout<<"large\n";
        else
            cout<<"best\n";
    }

    if(verbose)
        cout<<"Seeding approach ................. ";

    if(!strcmp(_seed_, "rnd"))
    {
        if(verbose)
            cout<<"rand\n";
        //sedFunc = &XBKMeans::rndTwin;
        seed    = _rnd_;
    }
    else if(!strcmp(_seed_, "kpp"))
    {
        if(verbose)
            cout<<"kpp\n";
        //sedFunc = &XBKMeans::kppTwin;
        seed    = _kpp_;
    }
    else if(!strcmp(_seed_, "non"))
    {
        if(verbose)
            cout<<"non\n";
        //sedFunc = NULL;
        seed    = _non_;
    }
    else
    {
        if(verbose)
            cout<<"kpp\n";
        //sedFunc = &XBKMeans::kppTwin;
        seed    = _kpp_;
    }

    if(verbose)
        cout<<"Distance function ................ l2\n";


    if(verbose)
        cout<<"Optimization function ............ ";
    if(!strcmp(crtrn, "i1"))
    {
        //optFunc = &XTKMeans::I1Func;
        myoptz     = _I1_;
        if(verbose)
            cout<<"I1\n";
    }
    else if(!strcmp(crtrn, "i2"))
    {
        optFunc = &XTKMeans::I2Func;
        myoptz  = _I2_;
        if(verbose)
            cout<<"I2\n";
    }
    else if(!strcmp(crtrn, "i4"))
    {
        //optFunc = &XTKMeans::I4Func;
        myoptz  = _I4_;
        if(verbose)
            cout<<"I4\n";
    }
    else if(!strcmp(crtrn, "e1") )
    {
        //optFunc = &XTKMeans::E1Func;
        myoptz  = _E1_;

        if(verbose)
            cout<<"e1\n";
    }
    else if(!strcmp(crtrn, "e2") )
    {
        //optFunc = &XTKMeans::E2Func;
        myoptz  = _E2_;

        if(verbose)
            cout<<"e2\n";
    }
    else
    {
        cout<<"Unkown optimize option '"<<crtrn<<"'!\n";
        exit(0);
    }

    return true;
}

int XTKMeans::getL2norms(const unsigned int n, const unsigned int d0)
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

double XTKMeans::pairwDst(const int nclust, int* tmpNs, int*mylabels)
{
    unsigned int i = 0, label = 0, cloc = 0, c0 = 0, nl = 0;
    double l2nsum = 0, pd = 0, *D;
    for(label = 0; label < nclust; label++)
    {
        cloc = label * ndim;
        D    = this->Ds+cloc;
        for(i = 0; i < this->ndim; i++)
        {
            pd += D[i]*D[i];
        }
    }
    for(i = 0; i < this->count; i++)
    {
        c0 = mylabels[i];
        nl = tmpNs[c0];
        l2nsum += nl * this->lens[i];
    }
    return (l2nsum-pd)/this->count;
}

int XTKMeans::clust(const unsigned int clust_num, const char *dstfn, const int verbose)
{
    unsigned int j = 0;
    allocMem(clust_num);
    this->getL2norms(this->count, this->ndim);

    for(j = 0; j < clust_num; j++)
    {
        this->infoMap[j].E = 0.0f;
        this->infoMap[j].n = 0;
        Es[j]              = 0;
    }

    if(verbose)
        cout<<"Clustering ....................... on progress\n";

    ///cout<<this->myoptz<<"\t"<<dataType<<endl;
    if(dataType)
    {
        if(this->myoptz == _I4_)
            incrOptzI4s(Es, clust_num);
        else if(this->myoptz == _I2_)
            incrOptzI2s(Es, clust_num);
    }
    else
    {
        if(this->myoptz == _I2_)
        {
            incrOptzI2(Es, clust_num);
            cout<<"i2\n";
        }
        else if(this->myoptz == _I4_)
        {
           incrOptzI4(Es, clust_num);
        }
    }
    this->calAVGDist(this->arrayD, clust_num, this->infoMap);
    ///cout<<"mvs0\t"<<mvs[0]<<"\tmvs1\t"<<mvs[1]<<endl;

    if(strlen(dstfn) > 0)
    {
        stable_sort(sorted_stack.begin(), sorted_stack.end(), NNItem::LLIDXcomparer);
        char infofn[512];
        VString::parsePath(infofn, dstfn);
        strcat(infofn, "_clust_info.txt");
        Print2Scrn::printvect(sorted_stack, infofn);
        XTKMeans::save_clust(dstfn);
    }
    Cleaner::clearVector(sorted_stack);
    return clust_num;
}

bool XTKMeans::incrOptzI2(double *Es, const unsigned int clustNum)
{
    unsigned int loc = 0, i = 0, j = 0, k = 0, t = 0, k0 = 0, loci = 0;
    unsigned int cpysize = sizeof(double)*this->ndim*clustNum;
    double rate = 1.0;
    cout<<"i am in optz i2\n";
    vector<int> crntClust;
    bool _UPDATED_ = false;

    for(i = 0; i < this->count; i++)
    {
        if(this->labels[i] == 0)
        {
            crntClust.push_back(i);
        }
        /// else
        ///cout<<"erro\t"<<labels[i]<<endl;
    }
    cout<<"Cluster num. ..................... "<<clustNum<<endl;

    unsigned int numb   = crntClust.size();
    int *tmpNs = new int[clustNum];
    double *tmpEs       = new double[clustNum];
    int *mylabels       = new int[numb];
    int *seeds          = new int[clustNum];
    int *rl             = new int[numb];
    double *tmpData     = new double[ndim];
    //ofstream outStrm;
    //outStrm.open("/home/wlzhao/xtk_brk5_distort.txt", ios::out);

    double maxEg = 0, tmpEg1 = 0, optEg = 0, tmpEg2 = 0, dct = 0, prvEg = 0, avgDst = 0;
    double dst = 0, minDst = 0, delta, delta1, mxEs, err = 0, delta0 = 1.0;
    unsigned int r = 0, tk = 0, iters = 0;
    clock_t t1;
    double t2 = 0;

    if(this->record == NULL)
    {
        this->record = new double[Round*6];
    }
    memset(this->record, 0, sizeof(double)*Round*6);
    cout<<"Initialize record ................ ";
    for(unsigned i = 0; i < Round; i++)
    {
        this->record[6*i]   = numeric_limits<double>::max();
        this->record[6*i+1] = -1;
        this->record[6*i+3] = numeric_limits<double>::max();
        this->record[6*i+4] = -1;
    }
    cout<<"done"<<endl;

    for(t = 0; t < XTKMeans::NTRAILS; t++)
    {
        memset(Ds,    0, cpysize);
        memset(tmpNs, 0, sizeof(int)*clustNum);
        memset(tmpEs, 0, sizeof(double)*clustNum);
        memset(mylabels, 0, sizeof(int)*numb);
        optEg = prvEg = 0;


        if(this->seed != _non_)
        {
            ///this->kppSeeds(clustNum, seeds, numb);
            this->rndSeeds(clustNum, seeds, numb);

            for(i = 0; i < clustNum; i++)
            {
                loc  = seeds[i]*this->ndim;
                loci = i*this->ndim;
                for(j = 0; j < this->ndim; j++)
                {
                    tmpCs[loci+j] = data[loc+j];
                }
            }

            for(i = 0; i < numb; i++)
            {
                r = crntClust[i];
                minDst = RAND_MAX;
                for(k = 0; k < clustNum; k++)
                {
                    dst = PQMath::l2d(tmpCs, k, data, r, this->ndim);
                    if(dst < minDst)
                    {
                        minDst = dst;
                        k0  = k;
                    }
                }

                loc  = r*this->ndim;
                loci = k0*ndim;
                for(j = 0; j < this->ndim; j++)
                {
                    this->Ds[loci+j] += data[loc+j];
                }

                tmpNs[k0]++;
                mylabels[i] = k0;
            }
        }
        else
        {
            for(i = 0; i < numb; i++)
            {
                r  = crntClust[i];
                k0 = rand()%clustNum;
                loc  = r*this->ndim;
                loci = k0*ndim;
                for(j = 0; j < this->ndim; j++)
                {
                    this->Ds[loci+j] += data[loc+j];
                }

                tmpNs[k0]++;
                mylabels[i] = k0;
            }
        }

        optEg = (this->*optFunc)(Ds, clustNum, this->ndim, tmpNs, tmpEs);
        ///cout<<"energy: "<<allEg - optEg<<endl;
        iters = 0;
        rate  = 1.0;
        do
        {
            t1 = clock();
            _UPDATED_ = false;
            for(i = 0; i < numb; i++)
            {
                rl[i] = i;
            }
            random_shuffle(rl, rl + numb);
            for(i = 0; i < numb; i++)
            {
                r     = rl[i];
                k0    = mylabels[r];
                loc   = r*this->ndim;
                for(j = 0; j < ndim; j++)
                {
                    tmpData[j] = data[loc+j];
                }
                dct   = lens[r];

                if(tmpNs[k0] <= 1)
                    continue;

                loci   = k0*ndim;
                tmpEg1 = 0;
                if(tmpNs[k0] > 1)
                {
                    for(j = 0; j < ndim; j++)
                    {
                        tmpEg1 += Ds[loci+j]*tmpData[j];
                    }
                    tmpEg1 = tmpEs[k0]*tmpNs[k0] - 2*tmpEg1 + dct;
                    tmpEg1 = tmpEg1/(tmpNs[k0]- 1);
                }

                delta1 = mxEs  = delta0 = 0;
                for(k = 0; k < clustNum; k++)
                {
                    if(k == k0)
                        continue;

                    loci   = k*ndim;
                    tmpEg2 = 0;

                    for(j = 0; j < ndim; j++)
                    {
                        tmpEg2 += Ds[loci+j]*tmpData[j];
                    }

                    tmpEg2 = tmpEs[k]*tmpNs[k] + 2*tmpEg2 + dct;
                    tmpEg2 = tmpEg2/(tmpNs[k]+1);

                    delta  = tmpEg2 - tmpEs[k] + tmpEg1 - tmpEs[k0];
                    if(delta > delta1)
                    {
                        delta0 = delta1;
                        delta1 = delta;
                        tk     = k;
                        mxEs   = tmpEg2;
                    }
                }

                if(delta1 > 0)
                {
                    tmpNs[k0]--;
                    tmpNs[tk]++;
                    mylabels[r] = tk;
                    loci = k0*ndim;
                    loc  = tk*ndim;
                    for(j = 0; j < this->ndim; j++)
                    {
                        this->Ds[loci+j] -= tmpData[j];
                        this->Ds[loc+j]  += tmpData[j];
                    }
                    optEg = optEg + delta1;
                    tmpEs[k0] = tmpEg1;
                    tmpEs[tk] = mxEs;
                    _UPDATED_ = true;
#ifdef _WATCH_
                    this->hits[r][iters] = tk; ///for observation
                    this->hists[iters]= this->hists[iters] + 1;
#endif
                }
            }///end-for(i)

            optEg = getI2(this->Ds, clustNum, ndim, tmpNs);
            //getI2(this->Ds, this->clnumb, this->ndim, Ns);
            //optEg = I2Func(this->Ds, clustNum, ndim, tmpNs, tmpEs);
            if(this->nLogs > 0 && iters < this->nLogs)
            {
                this->kmLogs[0][iters] = optEg;
                t2 = (clock() - t1 + 0.0)/CLOCKS_PER_SEC;
                if(iters == 0)
                {
                    this->kmLogs[1][iters] = t2;
                }
                else
                {
                    this->kmLogs[1][iters] = this->kmLogs[1][iters-1] + t2;
                }
            }
            double avgDistI1 = (allEg - optEg)/this->count;
            avgDst = pairwDst(clustNum, tmpNs, mylabels);
            cout<<"hello iters\t"<<iters<<"\t"<<avgDistI1<<"\t"<<this->kmLogs[1][iters]<<"\t"<<avgDst<<endl;
            if(this->record[6*iters] > avgDst)
                this->record[6*iters] = avgDst;
            if(this->record[6*iters+1] < avgDst)
                this->record[6*iters+1] = avgDst;
            this->record[6*iters+2] += avgDst;

            if(this->record[6*iters+3] > avgDistI1)
                this->record[6*iters+3] = avgDistI1;
            if(this->record[6*iters+4] < avgDistI1)
                this->record[6*iters+4] = avgDistI1;
            this->record[6*iters+5] += avgDistI1;

            err   = optEg - prvEg;
            prvEg = optEg;
            rate  = err/optEg;
            iters++;
        }
        while(_UPDATED_/**/ && iters < Round/**/);


        if(maxEg < optEg)
        {
            maxEg = optEg;
            memcpy(Ns,     tmpNs,    sizeof(unsigned int)*clustNum);
            memcpy(labels, mylabels, sizeof(int)*numb);
            memcpy(Es,     tmpEs, sizeof(double)*clustNum);
            memcpy(Cs,     Ds,    cpysize);
            memcpy(arrayD, Ds,    cpysize);
        }
    }///for(t)

    for(unsigned i = 0; i < Round; i++)
    {
        cout<<"avgI2: "<<i<<"\t"<<this->record[6*i]<<"\t"<<this->record[6*i+1]<<"\t"<<this->record[6*i+2]/NTRAILS<<endl;
        cout<<"avgI1: "<<i<<"\t"<<this->record[6*i+3]<<"\t"<<this->record[6*i+4]<<"\t"<<this->record[6*i+5]/NTRAILS<<endl;
    }

    optEg = I2Func(arrayD, clustNum, ndim, Ns, Es);
    ///this->saveHits("/home/wlzhao/look_hits.txt", Round);
    ///this->saveHist("/home/wlzhao/datasets/clust/result/sift100k_ITmoves.txt", Round);

    crntClust.clear();
    delete [] mylabels;
    delete [] tmpData;
    delete [] seeds;
    delete [] tmpNs;
    delete [] tmpEs;
    delete [] rl;
    mylabels = NULL;
    seeds    = NULL;
    tmpNs    = NULL;
    tmpEs    = NULL;
    rl       = NULL;
    tmpData  = NULL;

    return true;
}

bool XTKMeans::incrOptzI4(double *Es, const unsigned int clustNum)
{
    unsigned int loc = 0, i = 0, j = 0, k = 0, t = 0, k0 = 0, loci = 0;
    unsigned int cpysize = sizeof(double)*this->ndim*clustNum;
    double rate = 1.0, valI2 = 0;
    vector<int> crntClust;
    bool _UPDATED_ = false;

    for(i = 0; i < count; i++)
    {
        if(this->labels[i] == 0)
        {
            crntClust.push_back(i);
        }
    }
    unsigned int numb   = crntClust.size();
    int *tmpNs = new int[clustNum];
    double *tmpEs       = new double[clustNum];
    int *mylabels       = new int[numb];
    int *seeds          = new int[clustNum];
    int *rl             = new int[numb];
    float *pData        = NULL;

    double maxEg = 0, tmpEg1 = 0, optEg = 0, tmpEg2 = 0, dct = 0, prvEg = 0;
    double dst = 0, minDst = 0, delta, delta1, mxEs, err = 0;
    unsigned int r = 0, tk = 0, iters = 0;

    for(t = 0; t < XTKMeans::NTRAILS; t++)
    {
        memset(Ds,    0, cpysize);
        memset(tmpNs, 0, sizeof(int)*clustNum);
        memset(tmpEs, 0, sizeof(double)*clustNum);
        memset(mylabels, 0, sizeof(int)*numb);
        optEg = prvEg = 0;

        if(this->seed != _non_)
        {
            ///this->kppSeeds(clustNum, seeds, numb);
            this->rndSeeds(clustNum, seeds, numb);

            for(i = 0; i < clustNum; i++)
            {
                loc  = seeds[i]*this->ndim;
                loci = i*this->ndim;
                for(j = 0; j < this->ndim; j++)
                {
                    tmpCs[loci+j] = data[loc+j];
                }
            }

            for(i = 0; i < numb; i++)
            {
                r = crntClust[i];
                minDst = RAND_MAX;
                for(k = 0; k < clustNum; k++)
                {
                    dst = PQMath::l2d(tmpCs, k, data, r, this->ndim);
                    if(dst < minDst)
                    {
                        minDst = dst;
                        k0  = k;
                    }
                }

                loc  = r*this->ndim;
                loci = k0*ndim;
                for(j = 0; j < this->ndim; j++)
                {
                    this->Ds[loci+j] += data[loc+j];
                }

                tmpNs[k0]++;
                mylabels[i] = k0;
            }
        }
        else
        {
            for(i = 0; i < numb; i++)
            {
                r  = crntClust[i];
                k0 = rand()%clustNum;
                loc  = r*this->ndim;
                loci = k0*ndim;
                for(j = 0; j < this->ndim; j++)
                {
                    this->Ds[loci+j] += data[loc+j];
                }

                tmpNs[k0]++;
                mylabels[i] = k0;
            }
        }

        ///optEg = I4Func(Ds, clustNum, this->ndim, tmpNs, tmpEs);
        ///cout<<0<<"\t"<<optEg<<endl;
        iters = 0;
        rate  = 1.0;
        do
        {
            _UPDATED_ = false;
            for(i = 0; i < numb; i++)
            {
                rl[i] = i;
            }
            random_shuffle(rl, rl + numb);
            for(i = 0; i < numb; i++)
            {
                r     = rl[i];
                k0    = mylabels[r];
                loc   = r*this->ndim;
                pData = data+loc;
                dct   = lens[r];

                if(tmpNs[k0] <= 0)
                    continue;

                loci   = k0*ndim;
                tmpEg1 = 0;
                if(tmpNs[k0] > 1)
                {
                    for(j = 0; j < ndim; j++)
                    {
                        tmpEg1 += Ds[loci+j]*pData[j];
                    }
                    tmpEg1 = tmpEs[k0] - 2*tmpEg1 + dct;
                }
                ///cout<<tmpEg1<<endl;

                delta1 = mxEs  = 0;
                for(k = 0; k < clustNum; k++)
                {
                    if(k == k0)
                        continue;

                    loci   = k*ndim;
                    tmpEg2 = 0;

                    for(j = 0; j < ndim; j++)
                    {
                        tmpEg2 += Ds[loci+j]*pData[j];
                    }
                    tmpEg2 = tmpEs[k] + 2*tmpEg2 + dct;

                    delta  = sqrt(tmpEg2) + sqrt(tmpEg1) - sqrt(tmpEs[k]) - sqrt(tmpEs[k0]);
                    if(delta > delta1)
                    {
                        delta1 = delta;
                        tk     = k;
                        mxEs   = tmpEg2;
                    }
                }
                //cout<<delta1<<endl;

                if(delta1 > 0)
                {
                    tmpNs[k0]--;
                    tmpNs[tk]++;
                    mylabels[r] = tk;
                    loci = k0*ndim;
                    loc  = tk*ndim;
                    for(j = 0; j < this->ndim; j++)
                    {
                        this->Ds[loci+j] -= pData[j];
                        this->Ds[loc+j]  += pData[j];
                    }
                    optEg = optEg + delta1;
                    tmpEs[k0] = tmpEg1;
                    tmpEs[tk] = mxEs;
                    _UPDATED_ = true;
                }
            }///end-for(i)

            err   = optEg - prvEg;
            prvEg = optEg;
            valI2 = getI2(this->Ds, clustNum, ndim, tmpNs);
            //rate  = err/optEg;
            ///cout<<iters<<"\t"<<optEg<<"\t"<<(allEg-valI2)/this->count<<"\t"<<err<<endl;
            iters++;
        }
        while(_UPDATED_);

        if(maxEg < optEg)
        {
            maxEg = optEg;
            memcpy(Ns,     tmpNs,    sizeof(unsigned int)*clustNum);
            memcpy(labels, mylabels, sizeof(int)*numb);
            memcpy(Es,     tmpEs, sizeof(double)*clustNum);
            memcpy(Cs,     Ds,    cpysize);
            memcpy(arrayD, Ds,    cpysize);
        }
    }///for(t)

    optEg = getI4(arrayD, clustNum, ndim, Ns);

    crntClust.clear();
    delete [] mylabels;
    delete [] seeds;
    delete [] tmpNs;
    delete [] tmpEs;
    delete [] rl;
    mylabels = NULL;
    seeds    = NULL;
    tmpNs    = NULL;
    tmpEs    = NULL;
    rl       = NULL;

    return true;
}

bool XTKMeans::incrOptzI2s(double* Es, const unsigned int clustNum)
{
    unsigned int i = 0, j = 0, t = 0, centerid = 0, loc1 = 0, loc2 = 0, r = 0, k = 0, k0 = 0, loc = 0, loci = 0, tk = 0, iter = 0;
    double optEg = 0, maxEg = 0, dst1 = 0, dst2 = 0, tmpEg1 = 0, tmpEg2 = 0, delta = 0, delta1 = 0, maxEg1 = 0;
    int *tmplabels = NULL;
    unsigned int *tmpNs = NULL;
    bool _UPDATED_ = 0;
    int *seeds = NULL;
    int *rl = NULL;
    rl = new int[this->count];
    tmplabels = new int[this->count];
    tmpNs = new unsigned int[this->clnumb];
    seeds = new int[this->clnumb];
    int cpysize = sizeof(double)*this->ndim*this->clnumb;
    int numb = this->count;


    mvs[0] = mvs[1] = 0;

    for(i = 0; i < this->count; i++)
    {
        rl[i] = i;
    }

    maxEg1 = 0;
    for(t = 0; t < XTKMeans::NTRAILS; t++)
    {
        memset(Ds, 0, cpysize);
        memset(tmplabels, 0, sizeof(int)*numb);
        memset(tmpNs, 0, sizeof(unsigned int)*this->clnumb);

        ///**initialize **********************///


        if(this->seed != _non_)
        {
            if(this->seed == _kpp_)
                this->kppSeeds(this->clnumb, seeds, this->count);
            else if(this->seed == _rnd_)
                this->rndSeeds(this->clnumb, seeds, this->count);
            for(i = 0; i < this->clnumb; i++)
            {
                loc  = seeds[i];
                loci = i*this->ndim;
                for(j = sdata.col[loc]; j < sdata.col[loc+1]; j++)
                {
                    tmpCs[loci+sdata.index[j]] = sdata.data[j];
                }
            }
            for(i = 0; i < this->count; i++)
            {
                dst2 = RAND_MAX;
                for(j = 0; j < this->clnumb; j++)
                {
                    mvs[0]++;
                    dst1 = PQMath::l2d(tmpCs, j, sdata, i, this->ndim);
                    if(dst1 < dst2)
                    {
                        dst2 = dst1;
                        centerid = j;
                    }
                }
                if(dst2 != RAND_MAX)
                {
                    tmpNs[centerid]++;
                    tmplabels[i] = centerid;
                    loc2 = centerid*ndim;
                    loc1 = i*ndim;
                    for(j = sdata.col[i]; j < sdata.col[i+1]; j++)
                    {
                        Ds[loc2+sdata.index[j]] += sdata.data[j];
                    }
                }
            }
        }
        else
        {
            for(i = 0; i < this->count; i++)
            {
                k0 = rand()%this->clnumb;
                tmplabels[i] = k0;
                tmpNs[k0]++;
                loc1 = i*ndim;
                loc2 = k0*ndim;
                for(j = sdata.col[i]; j < sdata.col[i+1]; j++)
                {
                    Ds[loc2 + sdata.index[j]] += sdata.data[j];
                }
            }
        }
        optEg = 0;
        for(i = 0; i < this->clnumb; i++)
        {
            loci  = i*ndim;
            Es[i] = 0;
            for(j = 0; j < this->ndim; j++)
            {
                Es[i] += Ds[loci+j]*Ds[loci+j];
            }
            if(tmpNs[i] != 0)
                optEg += Es[i]/tmpNs[i];
        }

        ///*****incremental********///
        centerid = 0;

        do
        {
            _UPDATED_ = false;
            iter = 0;
            random_shuffle(rl, rl + this->count);
            for(i = 0; i < this->count; i++)
            {
                if(!SHUFFLE)
                {
                    r   = (unsigned int)floor((rand()/(RAND_MAX+1.0f))*this->clnumb);
                    r   = r>=this->clnumb?(this->clnumb-1):r;
                }
                else
                    r = rl[i];
                k0  = tmplabels[r];
                loc = r*this->ndim;
                if(tmpNs[k0] <= 1)
                    continue;

                loci   = k0*this->ndim;
                tmpEg1 = 0;
                for(j = sdata.col[r]; j < sdata.col[r+1]; j++)
                {
                    tmpEg1 += Ds[loci + sdata.index[j]]*sdata.data[j];
                }

                tmpEg1 = Es[k0] - 2*tmpEg1 + this->lens[r];
                delta1 = 0;
                for(k = 0; k < this->clnumb; k++)
                {
                    if(k == k0)
                        continue;
                    mvs[0]++;

                    loci   = k*ndim;

                    tmpEg2 = 0;
                    for(j = sdata.col[r]; j < sdata.col[r+1]; j++)
                    {
                        tmpEg2 += Ds[loci + sdata.index[j]]*sdata.data[j];
                    }

                    tmpEg2 = Es[k] + 2*tmpEg2 + this->lens[r];

                    //delta  = sqrt(tmpEg2) - sqrt(Es[k]) + sqrt(tmpEg1) - sqrt(Es[k0]);
                    delta = tmpEg2/(tmpNs[k]+1) - Es[k]/tmpNs[k] + tmpEg1/(tmpNs[k0]-1) - Es[k0]/tmpNs[k0];
                    if(delta > delta1)
                    {
                        delta1 = delta;
                        tk = k;
                        maxEg = tmpEg2;
                        if(iter > 5)
                            break;
                    }
                }

                if(delta1 > 0)
                {
                    _UPDATED_ = true;
                    mvs[1]++;
                    tmpNs[k0]--;
                    tmpNs[tk]++;
                    tmplabels[r] = tk;
                    Es[k0] = tmpEg1;
                    Es[tk] = maxEg;
                    loc1 = k0*ndim;
                    loc2 = tk*ndim;
                    for(j = sdata.col[r]; j < sdata.col[r+1]; j++)
                    {
                        Ds[loc1+sdata.index[j]] -= sdata.data[j];
                        Ds[loc2+sdata.index[j]] += sdata.data[j];
                    }
                    optEg = optEg+delta1;

                }
            }///for(i 0 - > count)
            /**
                centerid++;
             if(centerid > 30)
                break;
                /**/
        }
        while(_UPDATED_  && iter < 20);

        if(optEg > maxEg1)
        {
            maxEg1 = optEg;
            memcpy(Ns,     tmpNs,    sizeof(unsigned int)*clustNum);
            memcpy(labels, tmplabels, sizeof(int)*count);
            memcpy(arrayD, Ds,    cpysize);
        }
    }

    for(i = 0; i < this->clnumb; i++)
    {
        loci  = i*ndim;
        Es[i] = PQMath::dvec_norm(arrayD+loci, ndim, 2);
        optEg += Es[i];
        Es[i] *= Es[i];
    }


    delete []tmplabels;
    delete []tmpNs;
    delete []seeds;
    delete []rl;
    delete []tmpCs;
    tmpCs = NULL;
    tmplabels = NULL;
    tmpNs = NULL;
    seeds = NULL;
    rl = NULL;
    return true;
}

bool XTKMeans::incrOptzI4s(double* Es, const unsigned int clustNum)
{
    unsigned int i = 0, j = 0, t = 0, centerid = 0, loc1 = 0, loc2 = 0, r = 0, k = 0, k0 = 0, loc = 0, loci = 0, tk = 0;
    double optEg = 0, maxEg = 0, dst1 = 0, dst2 = 0, tmpEg1 = 0, tmpEg2 = 0, delta = 0, delta1 = 0, maxEg1 = 0;
    int *tmplabels = NULL;
    unsigned int *tmpNs = NULL;
    bool _UPDATED_ = 0;
    int *seeds = NULL;
    int *rl = NULL;
    rl = new int[this->count];
    tmplabels = new int[this->count];
    tmpNs = new unsigned int[this->clnumb];
    seeds = new int[this->clnumb];
    int cpysize = sizeof(double)*this->ndim*this->clnumb;
    int numb = this->count;


    mvs[0] = mvs[1] = 0;

    for(i = 0; i < this->count; i++)
    {
        rl[i] = i;
    }

    maxEg1 = 0;
    for(t = 0; t < XTKMeans::NTRAILS; t++)
    {
        memset(Ds, 0, cpysize);
        memset(tmplabels, 0, sizeof(int)*numb);
        memset(tmpNs, 0, sizeof(unsigned int)*this->clnumb);

        ///**initialize **********************///


        if(this->seed != _non_)
        {
            if(this->seed == _kpp_)
                this->kppSeeds(this->clnumb, seeds, this->count);
            else if(this->seed == _rnd_)
                this->rndSeeds(this->clnumb, seeds, this->count);
            for(i = 0; i < this->clnumb; i++)
            {
                loc  = seeds[i];
                loci = i*this->ndim;
                for(j = sdata.col[loc]; j < sdata.col[loc+1]; j++)
                {
                    tmpCs[loci+sdata.index[j]] = sdata.data[j];
                }
            }
            for(i = 0; i < this->count; i++)
            {
                dst2 = RAND_MAX;
                for(j = 0; j < this->clnumb; j++)
                {
                    mvs[0]++;
                    dst1 = PQMath::l2d(tmpCs, j, sdata, i, this->ndim);
                    if(dst1 < dst2)
                    {
                        dst2 = dst1;
                        centerid = j;
                    }
                }
                if(dst2 != RAND_MAX)
                {
                    tmpNs[centerid]++;
                    tmplabels[i] = centerid;
                    loc2 = centerid*ndim;
                    loc1 = i*ndim;
                    for(j = sdata.col[i]; j < sdata.col[i+1]; j++)
                    {
                        Ds[loc2+sdata.index[j]] += sdata.data[j];
                    }
                }
            }
        }
        else
        {
            for(i = 0; i < this->count; i++)
            {
                k0 = rand()%this->clnumb;
                tmplabels[i] = k0;
                tmpNs[k0]++;
                loc1 = i*ndim;
                loc2 = k0*ndim;
                for(j = sdata.col[i]; j < sdata.col[i+1]; j++)
                {
                    Ds[loc2 + sdata.index[j]] += sdata.data[j];
                }
            }
        }
        optEg = 0;
        for(i = 0; i < this->clnumb; i++)
        {
            loci  = i*ndim;
            Es[i] = 0;
            for(j = 0; j < this->ndim; j++)
            {
                Es[i] += Ds[loci+j]*Ds[loci+j];
            }
            optEg += sqrt(Es[i]);
        }

        ///*****incremental********///
        centerid = 0;

        do
        {
            _UPDATED_ = false;
            random_shuffle(rl, rl + this->count);
            for(i = 0; i < this->count; i++)
            {
                if(!SHUFFLE)
                {
                    r   = (unsigned int)floor((rand()/(RAND_MAX+1.0f))*this->clnumb);
                    r   = r>=this->clnumb?(this->clnumb-1):r;
                }
                else
                    r = rl[i];
                k0  = tmplabels[r];
                loc = r*this->ndim;
                if(tmpNs[k0] <= 1)
                    continue;

                loci   = k0*this->ndim;
                tmpEg1 = 0;
                for(j = sdata.col[r]; j < sdata.col[r+1]; j++)
                {
                    tmpEg1 += Ds[loci + sdata.index[j]]*sdata.data[j];
                }

                tmpEg1 = Es[k0] - 2*tmpEg1 + this->lens[r];
                delta1 = 0;
                for(k = 0; k < this->clnumb; k++)
                {
                    if(k == k0)
                        continue;
                    mvs[0]++;

                    loci   = k*ndim;

                    tmpEg2 = 0;
                    for(j = sdata.col[r]; j < sdata.col[r+1]; j++)
                    {
                        tmpEg2 += Ds[loci + sdata.index[j]]*sdata.data[j];
                    }

                    tmpEg2 = Es[k] + 2*tmpEg2 + this->lens[r];

                    delta  = sqrt(tmpEg2) - sqrt(Es[k]) + sqrt(tmpEg1) - sqrt(Es[k0]);
                    //delta = tmpEg2/(tmpNs[k]+1) - Es[k]/tmpNs[k] + tmpEg1/(tmpNs[k0]-1) - Es[k0]/tmpNs[k0];
                    if(delta > delta1)
                    {
                        delta1 = delta;
                        tk = k;
                        maxEg = tmpEg2;
                        break;
                    }
                }

                if(delta1 > 0)
                {
                    _UPDATED_ = true;
                    mvs[1]++;
                    tmpNs[k0]--;
                    tmpNs[tk]++;
                    tmplabels[r] = tk;
                    Es[k0] = tmpEg1;
                    Es[tk] = maxEg;
                    loc1 = k0*ndim;
                    loc2 = tk*ndim;
                    for(j = sdata.col[r]; j < sdata.col[r+1]; j++)
                    {
                        Ds[loc1+sdata.index[j]] -= sdata.data[j];
                        Ds[loc2+sdata.index[j]] += sdata.data[j];
                    }
                    optEg = optEg+delta1;
                }
            }///for(i 0 - > count)
            /**
                centerid++;
             if(centerid > 30)
                break;
                /**/
        }
        while(_UPDATED_ );

        if(optEg > maxEg1)
        {
            maxEg1 = optEg;
            memcpy(Ns,     tmpNs,    sizeof(unsigned int)*clustNum);
            memcpy(labels, tmplabels, sizeof(int)*count);
            memcpy(arrayD, Ds,    cpysize);
        }
    }

    for(i = 0; i < this->clnumb; i++)
    {
        loci  = i*ndim;
        Es[i] = PQMath::dvec_norm(arrayD+loci, ndim, 2);
        optEg += Es[i];
        Es[i] *= Es[i];
    }


    delete []tmplabels;
    delete []tmpNs;
    delete []seeds;
    delete []rl;
    delete []tmpCs;
    tmpCs = NULL;
    tmplabels = NULL;
    tmpNs = NULL;
    seeds = NULL;
    rl = NULL;
    return true;
}


bool XTKMeans::augOptz(double *Es, map<int, const char*> datafile, const unsigned int clustNum)
{

    ///1. Load unprocessed data,
    ///2. perform quantization
    ///3. perform 'refine' on these newly quantized data
    ///4. repeat step 1-3 until all data has been processed
    unsigned int id = 0, i = 0, j = 0, loc = 0, loc1 = 0, col = 0, row = 0, clid = 0;
    bool UPDATA = 0;
    double optEg;
    char filename[1024];
    map<int, const char*> dstfn;
    double *Ds = NULL;
    int *Nss;
    unsigned int *order = NULL;
    Ds = new double[clnumb*ndim];
    memset(Ds, 0, clnumb*ndim*sizeof(double));
    order = new unsigned int[datafile.size()];
    Nss = new int[clnumb];
    memset(Nss, 0, sizeof(int)*clnumb);
    for(i = 0; i < datafile.size(); i++)
    {
        order[i] = i;
    }
    /**   file procedure      **/

    for(i = 0; i < datafile.size(); i++)
    {
        sprintf(filename, "%s%s", "clustinfo_", datafile[i]);
        dstfn.insert(make_pair(i, filename));
    }
    ///load first group data and give first quantization;
    if(data != NULL)
    {
        delete [] data;
        data = NULL;
    }
    random_shuffle(order, order+datafile.size());
    id = order[0];
    data = IOAgent::load_fvecs(datafile[id], col, row);

    loc = 0;
    for(i = 0; i < col; i++)
    {
        clid = rand()%clnumb;
        loc1 = clid*ndim;
        for(j = 0; j < row; j++)
        {
            Ds[loc1 + j] += data[loc + j];
        }
        Nss[clid]++;
        loc += ndim;
    }

    /**  quantization  **/
    initialCenters(data, Ds, Nss, col, row, labels);
    optEg = extraRefine(data, Ds, Nss, col, row, labels);
    for(i = 1; i < datafile.size(); i++)
    {
        id = order[i];
        if(data != NULL)
        {
            delete [] data;
            data = NULL;
        }
        data = IOAgent::load_fvecs(datafile[id], col, row);
        initialCenters(data, arrayD, Ns, col, row, labels);
        optEg = extraRefine(data, arrayD, Ns, col, row,labels);
    }

    for(j = 1; j < Round; j++)
    {
        random_shuffle(order, order+datafile.size());
        for(i = 0; i < datafile.size(); i++)
        {
            id = order[i];
            if(data != NULL)
            {
                delete [] data;
                data = NULL;
            }
            data = IOAgent::load_fvecs(datafile[id], col, row);
            optEg = extraRefine(data, arrayD, Ns, col, row,labels);
        }
    }


    ///Optional steps:
    ///5. load a section of data in random
    ///6. perform 'refine' on this section of data
    ///7. repeat 5-6 for several rounds

    delete [] order;
    order = NULL;

    return true;
}

bool XTKMeans::augOptz(double *Es, const unsigned int clustNum)
{

    ///1. Load unprocessed data,
    ///2. perform quantization
    ///3. perform 'refine' on these newly quantized data
    ///4. repeat step 1-3 until all data has been processed
    unsigned int id = 0, i = 0, j = 0, loc = 0,
                 col = 0, row = 0, firstcol = 0;
    double optEg = 0;
    char filename[1024];
    int *tmplabels = NULL;
    float *tmpdata = NULL;
    double *Ds = NULL;
    unsigned int *order = NULL;
    int *Nss = NULL;
    firstcol = this->count;
    row = this->ndim;
    col = firstcol/Pubnum;
    order = new unsigned int[Pubnum];
    tmpdata = this->data;
    tmplabels = new int[firstcol];
    Ds = new double[clnumb*ndim];
    Nss = new int [clnumb];

    memset(arrayD, 0, sizeof(double)*clnumb*ndim);

    /**   file procedure      **/
    ///load first group data and give first quantization;
    for(i = 0; i < Pubnum; i++)
    {
        order[i] = i;
    }

    random_shuffle(order, order+Pubnum);

    /**   initialize    **/
    for(i = 0; i < Pubnum; i++)
    {
        id = order[i];
        data = tmpdata + id*col*row;
        initialCenters(data, arrayD, Ns, col, row, tmplabels + id*col);
    }
    if(firstcol%Pubnum != 0)
    {
        data = tmpdata + Pubnum*col*row;
        initialCenters(data, arrayD, Ns, firstcol%Pubnum, row, tmplabels + Pubnum*col);

    }
    ///memcpy(labels, tmplabels, firstcol*sizeof(int));

    /**  quantization  **/
    for(j = 0; j < Round; j++)
    {
        // UPDATA = 0;
        memset(Nss, 0, sizeof(int)*clnumb);
        memset(Ds, 0, sizeof(double)*clnumb*ndim);
        random_shuffle(order, order+Pubnum);
        for(i = 0; i < Pubnum; i++)
        {
            id = order[i];
            data = tmpdata + id*col*row;
            optEg = extraRefine(data, Ds, Nss, col, row, tmplabels + id*col);
            ///memcpy(tmplabels + id*col, labels, col*sizeof(int));
        }
        if(firstcol%Pubnum != 0)
        {
            data = tmpdata + Pubnum*col*row;
            optEg = extraRefine(data, Ds, Nss, firstcol%Pubnum, row, tmplabels + Pubnum*col);
            ///memcpy(tmplabels + Pubnum*col, labels, firstcol%Pubnum*sizeof(int));
        }
        for(i = 0; i < clnumb; i++)
        {
            loc = i*ndim;
            Ns[i] += Nss[i];
            for(unsigned int t = 0; t < ndim; t++)
            {
                arrayD[loc + t] += Ds[loc + t];
            }
        }
    }


    memcpy(labels, tmplabels, firstcol*sizeof(int));
    memset(Ns, 0, sizeof(int)*clnumb);
    memset(arrayD, 0, sizeof(double)*ndim*clnumb);
    for(i = 0; i < firstcol; i++)
    {
        Ns[labels[i]]++;
        data = tmpdata + i*row;
        id = labels[i]*row;
        for(j = 0; j < ndim; j++)
        {
            arrayD[id+j] += data[j];
        }
    }
    ///Optional steps:
    ///5. load a section of data in random
    ///6. perform 'refine' on this section of data
    ///7. repeat 5-6 for several rounds
    delete [] order;
    delete [] tmplabels;
    delete [] Nss;
    Nss = NULL;
    tmplabels = NULL;
    order = NULL;
    data = tmpdata;
    tmpdata = NULL;

    return true;
}

double XTKMeans::initialCenters(const float *data, double *Ds, int *Nss, const unsigned int col, const unsigned int row, int *labels)
{
    double dis, mindis;
    unsigned int i, j, loc, loc1, cpysize, tmpcenterid;
    unsigned int *centerid;
    double *centers;
    cpysize = clnumb*ndim;
    centers = new double[cpysize];
    centerid = new unsigned[col];
    /**  initial centers  **/
    //compute centers;
    loc = 0;
    if(this->seed != _non_)
    {
        for(i = 0; i < clnumb; i++)
        {
            for(j = 0; j < row; j++)
            {
                centers[loc + j] = Ds[loc + j]/Nss[i];
            }
            loc += ndim;
        }
        //assign centers to every point;
        tmpcenterid = 0;
        loc = 0;
        for(i = 0; i < col; i++)
        {
            mindis = RAND_MAX;
            for(j = 0; j < clnumb; j++)
            {
                dis = PQMath::l2d(centers, j, data, i, ndim);
                if(dis < mindis)
                {
                    mindis = dis;
                    tmpcenterid = j;
                }
            }
            centerid[i] = tmpcenterid;
            Ns[centerid[i]]++;
            loc1 = centerid[i]*ndim;
            labels[i] = centerid[i];
            for(j = 0; j < row; j++)
            {
                arrayD[loc1+j] += data[loc+j];
            }
            loc += row;
        }
        //end initial centers;
    }
    else
    {
        for(i = 0; i < col; i++)
        {
            tmpcenterid = rand()%clnumb;
            Ns[tmpcenterid]++;
            loc  = tmpcenterid*ndim;
            loc1 = i*ndim;
            labels[i]  = tmpcenterid;
            for(j = 0; j < ndim; j++)
            {
                arrayD[loc + j] += data[loc1+j];
            }
        }
    }

    delete [] centerid;
    delete [] centers;
    return 0;
}

double XTKMeans::extraRefine(const float *data, double *Ds, int *Nss, const unsigned int col, const unsigned int row, int *labels)
{
    double optEg = 0, tmpEg1 = 0, delta = 0, tmpEg2 = 0, delta1 = 0, mxEs = 0;
    unsigned int *rl = NULL,  k = 0, tk = 0, k0 = 0, loc = 0, loc1 = 0;
    unsigned int i = 0, j = 0, t = 0, iter = 0, r = 0;
    double *aD = new double[clnumb*ndim];
    rl = new unsigned int[col];
    int *sn = new int[clnumb];
    memcpy(aD, arrayD, sizeof(double)*clnumb*ndim);
    memcpy(sn, Ns, sizeof(int)*clnumb);

    /**   incremental      **/

    optEg = 0;
    for(i = 0; i < clnumb; i++)
    {
        loc1  = i*row;
        Es[i] = PQMath::dvec_norm(aD+loc1, ndim, 2);
        if(Ns[i] > 1)
        {
            Es[i] = Es[i]*(Es[i]/Ns[i]);
        }
        optEg += Es[i];
    }
    for(i = 0; i < col; i++)
    {
        rl[i] = i;
    }
    for(t = 0; t < NRef; t++)
    {
        random_shuffle(rl, rl + col);
        for(iter = 0; iter < col; iter++)
        {
            if(!SHUFFLE)
            {
                r   = (unsigned int)floor((rand()/(RAND_MAX+1.0f))*col);
                r   = r>=col?(col-1):r;
            }
            else
                r = rl[iter];
            k0  = labels[r];
            loc = r*row;

            if((Ns[k0]+Nss[k0]) <= 1)
                continue;

            loc1   = k0*row;
            tmpEg1 = I2FastM(aD+loc1, data+loc, row, sn[k0]);
            delta1 = mxEs = 0;
            for(k = 0; k < clnumb; k++)
            {
                if(k == k0)
                    continue;

                loc1   = k*row;
                tmpEg2 = I2FastP(aD+loc1, data+loc, row, sn[k]);
                delta  = tmpEg2 - Es[k] + tmpEg1 - Es[k0];

                if(delta > delta1)
                {
                    delta1 = delta;
                    tk = k;
                    mxEs = tmpEg2;
                }
            }

            if(delta1 > 0)
            {
                Nss[k0]--;
                sn[k0]--;
                sn[tk]++;
                Nss[tk]++;
                labels[r] = tk;
                Es[k0] = tmpEg1;
                Es[tk] = mxEs;
                for(j = 0; j < row; j++)
                {
                    Ds[k0*row+j] -= data[loc+j];
                    Ds[tk*row+j] += data[loc+j];
                    aD[k0*row+j] -= data[loc+j];
                    aD[tk*row+j] += data[loc+j];
                }
                optEg = optEg+delta1;
            }
        }///for(iter)
    }  ///for(t)
    /** incremental end **/
    delete [] rl;
    rl = NULL;
    delete [] aD;
    aD = NULL;
    return optEg;
}

double XTKMeans::I2Func(const double *Ds, const unsigned int k, const unsigned int dim,
                        const int *tmpns, double *dsts)
{
    double sumDst = 0;
    unsigned int i = 0, j = 0, loc = 0;
    for(j = 0; j < k; j++)
    {
        dsts[j] = 0;
        if(tmpns[j] == 0)
            continue;

        loc = dim*j;

        for(i = 0; i < dim; i++)
        {
            dsts[j] += Ds[loc+i]*Ds[loc+i];
        }
        dsts[j]= dsts[j]/tmpns[j];
        sumDst += dsts[j];
    }
    return sumDst;
}

/**
double XTKMeans::I4Func(const double *Ds, const unsigned int k, const unsigned int dim,
                        const int *tmpns, double *dsts)
{
    double sumDst = 0;
    unsigned int i = 0, j = 0, loc = 0;

    for(j = 0; j < k; j++)
    {
        dsts[j] = 0;
        if(tmpns[j] == 0)
            continue;

        loc = dim*j;

        for(i = 0; i < dim; i++)
        {
            dsts[j] += Ds[loc+i]*Ds[loc+i];
        }
        sumDst += sqrt(dsts[j]);
    }
    return sumDst;
}
/**/

bool XTKMeans::rndSeeds(const unsigned int k, int rseeds[], const unsigned int bound)
{
    srand(time(NULL));
    unsigned int NITER, it = 0, sed0 = (unsigned int)time(NULL);
    float v     = rand_r(&sed0)/(RAND_MAX+1.0f);
    bool FAILED = false;
    unsigned int i = 1, j = 0;
    int r = 0;

    rseeds[0] = (int)floor(v*bound);
    NITER = 20*k;

    while(i < k && it < NITER)
    {
        v      = rand_r(&sed0)/(RAND_MAX+1.0f);
        r      = (int)floor(v*bound);
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
bool XTKMeans::kppSeeds(const unsigned int k, int rseeds[], const unsigned int bound)
{
    unsigned long i = 0, j = 0;
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


void XTKMeans::saveCenters(const char *dstfn, bool append)
{
    unsigned int clabel = 0, j = 0, loc = 0, rCNum = 0;

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
}

void XTKMeans::saveHits(const char *dstFn, int cN)
{
#ifdef _WATCH_
    ofstream *outStrm = new ofstream(dstFn, ios::out);
    unsigned int i = 0, j = 0;
    for(i = 0; i < this->count; i++)
    {
        (*outStrm)<<i;
        for(j = 0; j < 20; j++)
        {
            (*outStrm)<<"\t"<<this->hits[i][j];
        }
        (*outStrm)<<endl;
    }

    outStrm->close();
#endif
}

void XTKMeans::saveHist(const char *dstFn, int cN)
{
#ifdef _WATCH_
    ofstream *outStrm = new ofstream(dstFn, ios::out);
    unsigned int i = 0;
    for(i = 0; i < cN; i++)
    {
        (*outStrm)<<i<<"\t"<<this->hists[i]<<endl;
    }

    outStrm->close();
#endif
}

int XTKMeans::fetchCenters(float *centers)
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

        idxi += this->ndim;
    }
    return rCNum;
}

void XTKMeans::test()
{

    const char *srcFn1 = "/home/wlzhao/datasets/bignn/sift1m/sift_base.txt";
    const char *srcFn2 = "/home/wlzhao/datasets/bignn/glove/glove1m_base.txt";
    //const char *srcFn3 = "/home/wlzhao/datasets/bignn/mat/sift_learn.txt";
    //const char *srcFn3 = "/home/wlzhao/datasets/bignn/sift1m/sift_learn.txt";
    const char *srcFn3 = "/home/rqchen/dataset/sift_learn.txt";
    const char *srcFn4 = "/home/wlzhao/datasets/clust/SUSY.csv";
    const char *srcFn5 = "/home/wlzhao/datasets/clust/uscensus1990_mat2m.txt";
    const char *srcFn6 = "/home/wlzhao/datasets/clust/msd-rh.txt";
    const char *srcFn7 = "/home/wlzhao/datasets/bignn/rand/rand100k20d.txt";

    const char *dstFn1 = "/home/wlzhao/datasets/clust/sift1m_clust.txt";
    const char *dstFn2 = "/home/wlzhao/datasets/clust/glove1m_clust_x.txt";
    //const char *dstFn3 = "/home/wlzhao/datasets/clust/sift_learn_1k_boostkm_clust.txt";
    const char *dstFn3 = "/home/rqchen/dataset/rslt/sift_learn_1k_boostkm_clust.txt";
    const char *dstFn4 = "/home/wlzhao/datasets/clust/SUSY_clust.txt";
    const char *dstFn5 = "/home/wlzhao/datasets/clust/uscensus1990_mat2m_clust.txt";
    const char *dstFn7 = "/home/wlzhao/datasets/clust/rand100k20d_1k_clust.txt";
    const char *logFn4 = "/home/wlzhao/datasets/clust/result/susy_xtkmeans#_iters.txt";
    const char *dstFn6 = "/home/wlzhao/datasets/clust/msd-rh_clust.txt";
    const char *logFn6 = "/home/wlzhao/datasets/clust/result/msd-rh_xtkmeans#_iters.txt";

    unsigned int row = 0, dim = 0, cnum = 1024, i = 0, sz = 1000;
    for(i = 0; i < 1; i++)
    {
        XTKMeans *mykm = new XTKMeans();
        mykm->setLogOn(64);
        ///mykm->realSz = sz;
        mykm->buildcluster(srcFn3, dstFn3, "non", "large", "i2", cnum, false);
        ///cout<<"\nData size: "<<mykm->count<<endl;
        ///mykm->saveLogs(logFn6, 90);
        delete mykm;
        ///sz = sz*10;
        ///cnum = cnum*2;
    }

}
