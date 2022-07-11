#include "rrclust.h"

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

const int RRClust::NTRAILS = 800;
const int RRClust::NRef = 1;
const int RRClust::Round = 20;
const int RRClust::Pubnum = 100;
const int RRClust::NIter0 = 100;

RRClust::RRClust()
{
    Ds = Cs = tmpCs   = NULL;
    kmMtd   = _rrkmn_;
    strcpy(mthStr, "_rr_");
    this->_REFER_ = false;
    cout<<"Method ........................... Round-Robin K-Means\n";
}

RRClust::~RRClust()
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

bool  RRClust::init(const char *srcfn)
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

bool RRClust::allocMem(const unsigned int clustNum)
{
    this->Ds     = new double[clustNum*this->ndim];
    this->Cs     = new double[clustNum*this->ndim];
    this->tmpCs  = new double[clustNum*this->ndim];
    return true;
}

bool RRClust::deallocMem()
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

    return true;
}

bool  RRClust::init(float *mat, const int row, const int dim)
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

bool RRClust::refresh()
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

    ///memset(this->visited, 0, nThrd*sizeof(unsigned char));
    return true;
}

bool  RRClust::config(const char *_seed_, const char *lg_first, const char *crtrn, int verbose)
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
        //optFunc = &RRClust::I1Func;
        myoptz     = _I1_;
        if(verbose)
            cout<<"I1\n";
    }
    else if(!strcmp(crtrn, "i2"))
    {
        optFunc = &RRClust::I2Func;
        myoptz  = _I2_;
        if(verbose)
            cout<<"I2\n";
    }
    else if(!strcmp(crtrn, "i4"))
    {
        //optFunc = &RRClust::I4Func;
        myoptz  = _I4_;
        if(verbose)
            cout<<"I4\n";
    }
    else if(!strcmp(crtrn, "e1") )
    {
        //optFunc = &RRClust::E1Func;
        myoptz  = _E1_;

        if(verbose)
            cout<<"e1\n";
    }
    else if(!strcmp(crtrn, "e2") )
    {
        //optFunc = &RRClust::E2Func;
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

int RRClust::clust(const unsigned int clust_num, const char *dstfn, const int verbose)
{
    unsigned int j = 0;
    allocMem(clust_num);
    for(j = 0; j < clust_num; j++)
    {
        this->infoMap[j].E = 0.0f;
        this->infoMap[j].n = 0;
        Es[j]              = 0;
    }

    if(verbose)
        cout<<"Clustering ....................... on progress\n";

    incrOptzI2(Es, clust_num);

    this->calAVGDist(this->arrayD, clust_num, this->infoMap);
    cout<<"mvs0\t"<<mvs[0]<<"\tmvs1\t"<<mvs[1]<<endl;


    return clust_num;
}

bool RRClust::incrOptzI2(double *Es, const unsigned int clustNum)
{
    int loc = 0, i = 0, j = 0, c1 = 0, c2 = 0, loc1 = 0;
    double optEg = 0;
    memset(Ns, 0, sizeof(int)*clustNum);
    memset(arrayD, 0, sizeof(double)*clustNum*ndim);
    int *cs = new int[clustNum];

    for(i = 0; i < clustNum; i++)
    {
       cs[i] = i;
    }

    for(i = 0; i < count; i++)
    {
        c1 = rand()%this->clnumb;
        loc = c1*ndim;
        for(j = 0; j < ndim; j++)
        {
            arrayD[loc + j] += this->data[loc1 + j];
        }
        Ns[c1] = Ns[c1] + 1;
        this->labels[i] = c1;
        loc1 += ndim;
    }

    /** two cluster move **/

    for(i = 0; i < NTRAILS; i++)
    {
        random_shuffle(cs, cs+clustNum);
        for(j = 0; j < this->clnumb; j++)
        {
            c1 = rand()%clnumb;
            c2 = (c1+1)%clnumb;
            c1 = cs[c1];
            c2 = cs[c2];
            BiOptz(c1, c2, arrayD+c1*ndim, arrayD+c2*ndim, Ns[c1], Ns[c2]);
        }
    }
    this->calAVGDist(arrayD, clustNum, this->infoMap);

    delete [] cs;
    cs = NULL;

    return true;
}


bool RRClust::BiOptz(int clbl1, int clbl2, double *D1, double *D2, int &n1, int &n2)
{
    unsigned int i = 0, j = 0, loc = 0, iter = 0;
    double optEg = 0, tmpEg1 = 0, tmpEg2 = 0, delta = 0;
    bool UPDATED = false;
    vector<int> crntClust;
    unsigned int num = n1 + n2;

    int *mylabels = new int[num];

    for(i = 0, j = 0; i < count; i++)
    {
        if(this->labels[i] == clbl1)
        {
           crntClust.push_back(i);
           mylabels[j] = 0;
           j++;
        }else if(this->labels[i] == clbl2)
        {
           crntClust.push_back(i);
           mylabels[j] = 1;
           j++;
        }
    }
    //cout<<num<<"\t"<<j<<endl;

    optEg = 0;
    for(i = 0; i < ndim; i++)
    {
        optEg += D1[i]*D1[i]/n1 + D2[i]*D2[i]/n2;
    }

    do
    {
        UPDATED = false;
        for(i = 0; i < num ; i++)
        {
            loc = crntClust[i]*ndim;

            if(mylabels[i] == 0)
            {
                tmpEg1 = this->I2FastM(D1, data + loc, ndim, n1);
                tmpEg2 = this->I2FastP(D2, data + loc, ndim, n2);
            }
            else if(mylabels[i] == 1)
            {
                tmpEg1 = this->I2FastP(D1, data + loc, ndim, n1);
                tmpEg2 = this->I2FastM(D2, data + loc, ndim, n2);
            }

            delta = tmpEg2 + tmpEg1 - optEg;

            if(delta > 0)
            {
                optEg += delta;
                if(mylabels[i] == 0)
                {
                    n1--;
                    n2++;
                    for(j = 0; j < ndim; j++)
                    {
                        D1[j] -= data [j + loc];
                        D2[j] += data [j + loc];
                    }
                    mylabels[i] = 1;
                }
                else if(mylabels[i] == 1)
                {
                    n1++;
                    n2--;
                    for(j = 0; j < ndim; j++)
                    {
                        D1[j] += data [j + loc];
                        D2[j] -= data [j + loc];
                    }
                    mylabels[i] = 0;
                }
                UPDATED = true;
            }//if(delta)
        }
        iter++;
        //if(iter > RRClust::NIter0)
        //    break;
    }while(UPDATED);

    for(i = 0; i < num; i++)
    {
        if(mylabels[i] == 0)
        {
           j = crntClust[i];
           this->labels[j] = clbl1;
        }else
        {
           j = crntClust[i];
           this->labels[j] = clbl2;
        }
    }

    /**/
    delete [] mylabels;
    mylabels = NULL;
    crntClust.clear();

    return true;
}

bool RRClust::BiMove(vector<int> &cluster1, vector<int > &cluster2, double *D1, double *D2)
{
    unsigned int i = 0, j = 0, loc = 0, iter = 0, n1 = 0, n2 =0;
    double optEg = 0, tmpEg1 = 0, tmpEg2 = 0, delta = 0;
    bool UPDATED = false;
    vector<int > clust;

    n1 = cluster1.size();
    n2 = cluster2.size();
    int num = cluster1.size() + cluster2.size();
    clust.reserve(num);
    int *labels = new int[num];

    for(i = 0; i < cluster1.size(); i++)
    {
        clust.push_back(cluster1[i]);
        labels[i] = 0;
    }

    j = cluster1.size();
    for(i = 0; i <cluster2.size(); i++, j++)
    {
        clust.push_back(cluster2[i]);
        labels[j] = 1;
    }


    optEg = 0;
    for(i = 0; i < ndim; i++)
    {
        optEg += D1[i]*D1[i]/n1 + D2[i]*D2[i]/n2;
    }

    do
    {
        UPDATED = false;
        for(i = 0; i < num ; i++)
        {
            loc = clust[i]*ndim;

            if(labels[i] == 0)
            {
                tmpEg1 = this->I2FastM(D1, data + loc, ndim, n1);
                tmpEg2 = this->I2FastP(D2, data + loc, ndim, n2);
            }
            else if(labels[i] == 1)
            {
                tmpEg1 = this->I2FastP(D1, data + loc, ndim, n1);
                tmpEg2 = this->I2FastM(D2, data + loc, ndim, n2);
            }

            delta = tmpEg2 + tmpEg1 - optEg;

            if(delta > 0)
            {
                optEg += delta;
                if(labels[i] == 0)
                {
                    n1--;
                    n2++;
                    for(j = 0; j < ndim; j++)
                    {
                        D1[j] -= data [j + loc];
                        D2[j] += data [j + loc];
                    }
                    labels[i] = 1;
                }
                else if(labels[i] == 1)
                {
                    n1++;
                    n2--;
                    for(j = 0; j < ndim; j++)
                    {
                        D1[j] += data [j + loc];
                        D2[j] -= data [j + loc];
                    }
                    labels[i] = 0;
                }
                UPDATED = true;
            }//if(delta)
        }
        iter++;
        //if(iter > RRClust::NIter0)
        //    break;
    }while(UPDATED);

    cluster1.clear();
    cluster2.clear();

    for(i = 0; i < num; i++)
    {
        if(labels[i] == 0)
        {
            cluster1.push_back(clust[i]);
        }
        else
            {
            cluster2.push_back(clust[i]);
        }
    }

    delete [] labels;
    labels = NULL;

    return true;
}

double RRClust::initialCenters(const float *data, double *Ds, int *Nss, const unsigned int col, const unsigned int row, int *labels)
{
    double dis, mindis;
    unsigned int i , j, loc, loc1, cpysize, tmpcenterid;
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


double RRClust::I2Func(const double *Ds, const unsigned int k, const unsigned int dim,
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


void RRClust::saveCenters(const char *dstfn, bool append)
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

int RRClust::fetchCenters(float *centers)
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

void RRClust::test()
{
    const char *srcfn = "/home/wlzhao/datasets/bignn/mat/sift_learn.txt";
    const char *dstfn = "sift_learndst.txt";
    ///const char *srcfn = "norm15dat/tr41.mat";
    ///const char *dstfn = "tr41dst.txt";

    RRClust *mykm = new RRClust();
    mykm->buildcluster(srcfn, dstfn, "non", "large", "i2", 1024, false);
    //mykm->saveCenters(ct, 0);
    delete mykm;
}
