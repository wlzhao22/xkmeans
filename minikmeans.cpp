#include "minikmeans.h"
#include "vstring.h"
#include "pqmath.h"
#include "ioagent.h"
#include "sparsematrix.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <cstring>

using namespace std;

MiniKMeans::MiniKMeans()
{
    // B0 = 100;
    rate   = 100;
    T0     = 100;
    _INIT_ = false;
    kmMtd  = _mnkmn_;

    strcpy(mthStr, "_mnk_");
    cout<<"Method ...........................  Mini-KMeans\n";
}

bool MiniKMeans::init(const char *srcfn)
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
        //this->data = IOAgent::loadSparse(srcfn, this->count, this-ndim);
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

    B0 = this->count/rate;

    this->_INIT_  = true;
    this->_REFER_ = false;

    return true;
}

bool MiniKMeans::init(float *mat, const int row, const int dim)
{
    this->refresh();

    this->data  = mat;
    this->count = row;
    this->ndim  = dim;

    B0 = this->count/rate;

    this->_INIT_  = true;
    this->_REFER_ = true;
    return true;
}

bool MiniKMeans::refresh()
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

    return true;
}

bool MiniKMeans::config(const char *_seed_, const char *crtrn, const char *lg_first, int verbose)
{
    if(verbose)
        cout<<"Distance function ................ l2\n";

    if(verbose)
        cout<<"Seeds ............................ ";
    if(!strcmp(_seed_, "rnd"))
    {
        if(verbose)
            cout<<"rand\n";
    }
    else if(!strcmp(_seed_, "kpp"))
    {

        if(verbose)
            cout<<"kpp\n";
    }
    else
    {
        if(verbose)
            cout<<"non\n";
    }

    return true;
}

int MiniKMeans::clust(const unsigned int clust_num, const char *dstfn, const int verbose)
{
    unsigned int k0 = 0, d = 0;
    double *centers = new double[clust_num*ndim];
    unsigned int i, j, loc, loc1, iter, loc2;
    unsigned int *tmpCenters = new unsigned int[clust_num];
    unsigned int *id = new unsigned int[count];
    unsigned int *vc = new unsigned int[clust_num];
    unsigned int *batch   = new unsigned[B0];
    unsigned int *b2label = new unsigned[B0];
    double *tmpDat = new double[ndim];
    double ita = 0, minDst = 0, dis = 0;

    this->clnumb = clust_num;

    for(i = 0; i < count; i++)
    {
        id[i] = i;
    }

    random_shuffle(id, id + count);
    for(i = 0; i < clust_num; i++)
    {
        tmpCenters[i] = id[i];
    }

    loc = 0;
    for(i = 0; i < clust_num; i++)
    {
        vc[i] = 0;
        loc1 = ndim*tmpCenters[i];
        for(j = 0; j < ndim; j++)
        {
            centers[loc + j] = data[loc1 + j];
        }
        loc += ndim;
    }

    //ofstream outStrm;
    //outStrm.open("/home/wlzhao/mnkm_non_fun.txt");
    cout<<"Start clustering ................. \n";

    cout<<"step-1 begin"<<endl;
    for(iter = 0; iter < T0; iter++)
    {
        ///select batch;
        random_shuffle(id, id + count);
        for(i = 0; i < B0; i++)
        {
            batch[i] = id[i];
        }

        ///find nearst centers;
        for(i = 0; i < B0; i++)
        {
            b2label[i] = 0;
            minDst = RAND_MAX;
            loc    = batch[i]*ndim;
            for(d = 0; d < ndim; d++)
            {
                tmpDat[d] = data[loc+d];
            }
            for(j = 0; j < clust_num; j++)
            {
                dis = PQMath::l2d(centers, j, tmpDat, 0, ndim);
                if(dis < minDst)
                {
                    b2label[i] = j;
                    minDst = dis;
                }
            }
        }//for(i)

        ///centers learning by gradient descent;
        for(i = 0; i < B0; i++)
        {
            loc     = b2label[i];
            loc1    = loc*ndim;
            vc[loc] = vc[loc] + 1;
            ita     = 1.0/vc[loc];
            loc2    = batch[i]*ndim;
            for(j = 0; j < ndim; j++)
            {
                centers[loc1 + j] = (1 - ita)*centers[loc1 + j] + ita*data[loc2 + j];
            }
        }
    }//for(iter)
    cout<<"step-1 end"<<endl;
    //outStrm.close();

    for(i = 0; i < clust_num; i++)
    {
        this->infoMap[i].n = 0;
    }

    ///get clust ;
    memset(this->Ns, 0, clust_num*sizeof(int));
    memset(this->arrayD, 0, clust_num*sizeof(double)*ndim);
    cout<<"step-2 begin"<<endl;
    for(i = 0; i < count; i++)
    {
        minDst = RAND_MAX;
        loc    = i*ndim;
        for(d = 0; d < ndim; d++)
        {
            tmpDat[d] = data[loc+d];
        }
        for(j = 0; j < clust_num; j++)
        {
            dis = PQMath::l2d(centers, j, tmpDat, 0, ndim);
            if(dis < minDst)
            {
                minDst = dis;
                k0     = j;
            }
        }
        labels[i] = k0;
        this->infoMap[k0].n = this->infoMap[k0].n + 1;
        this->Ns[k0] += 1;
        for(j = 0; j < ndim; j++)
        {
            this->arrayD[k0*ndim+j] += tmpDat[j];
        }
    }
    cout<<"step-2 end"<<endl;
    ///save_clust(dstfn);

    this->calAVGDist(this->arrayD, this->clnumb, infoMap);
    memcpy(this->arrayD, centers, sizeof(double)*clust_num*ndim);


    delete [] tmpCenters;
    delete [] id;
    delete [] vc;
    delete [] batch;
    delete [] b2label;
    delete [] centers;
    delete [] tmpDat;

    id =  vc = NULL;
    batch = b2label = tmpCenters = NULL;
    centers = NULL;
    return clust_num;
}

void MiniKMeans::saveCenters(const char *dstfn, bool append)
{
    unsigned int clabel = 0, j = 0, loc = 0, rCNum = 0;

    for(clabel = 0; clabel < this->clnumb; clabel++)
    {
        if(this->infoMap[clabel].n > 0)
        {
            rCNum++;
        }
    }

    cout<<clnumb<<"\t"<<rCNum<<endl;

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
                (*outStrm)<<this->arrayD[loc+j]<<" ";
            }
            (*outStrm)<<endl;
        }
    }

    outStrm->close();
    cout<<"done\n";
}

int MiniKMeans::fetchCenters(float *centers)
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
                centers[idxi + j] = (float)this->arrayD[loc+j];
            }
            rCNum++;
        }

        idxi += this->ndim;
    }
    return rCNum;
}

void MiniKMeans::setTO(const unsigned int i)
{
    T0 = i;
}



MiniKMeans::~MiniKMeans()
{


}

void MiniKMeans::pctest()
{
    const char *srcMatFn1 = "/home/pclin/nns_benchmark-master/data/gist/gist_base.fvecs";
    const char *dstFn     = "gist_test.txt";
    const char *dstEtc    = "distortion.txt";
    const unsigned int iter = 80;
    unsigned int clustnum = 1024;
    unsigned int it = 0;
    MiniKMeans *mykmean = new MiniKMeans();
    mykmean->setLogOn(iter);
    for(it = 0; it < iter; it++)
    {
        mykmean->setTO(it);
        mykmean->buildcluster(srcMatFn1, dstFn, "non", "xx", "xx", clustnum, false);
    }
    delete mykmean;
}

void MiniKMeans::test()
{
    const char *srcMatFn1 = "/home/wlzhao/datasets/bignn/sift1m/sift_base.txt";
    const char *srcMatFn2 = "/home/wlzhao/datasets/bignn/sift1m/sift_learn.txt";
    const char *srcMatFn3 = "itm_vlad_hesaff_flickr10m.itm";
    const char *srcFn6    = "/home/wlzhao/datasets/clust/msd-rh.txt";

    const char *dstFn6    = "/home/wlzhao/datasets/clust/msd-rh_clust.txt";
    const char *dstFn2    = "/home/wlzhao/datasets/clust/siftlean_minikm_clust.txt";
    const char *dstFn = "test_minik.txt";

    int i = 1024;
    ///do
    {
        MiniKMeans *mykmean = new MiniKMeans();
        mykmean->buildcluster(srcMatFn2, dstFn6, "non", "xx", "xx", i, false);
        delete mykmean;
        i = i*2;
    }///while(i <= 8192);

}
