#include "seqkmeans.h"
#include "vstring.h"
#include "ioagent.h"
#include "pqmath.h"

#include <iostream>
#include <algorithm>
#include <cstring>
#include <cassert>
#include <fstream>
#include <cmath>

using namespace std;


const int SeqKMeans::nTrails = 128;

SeqKMeans::SeqKMeans()
{
    kmMtd   = _olkmn_;
    strcpy(mthStr, "_olk_");
    this->_REFER_ = false;
    cout<<"Method ........................... OL K-Means\n";
}

SeqKMeans::~SeqKMeans()
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

bool  SeqKMeans::init(const char *srcfn)
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



bool  SeqKMeans::init(float *mat, const int row, const int dim)
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

bool SeqKMeans::refresh()
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

bool  SeqKMeans::config(const char *_seed_, const char *lg_first, const char *crtrn, int verbose)
{
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
    return true;
}

int SeqKMeans::getL2norms(const unsigned int n, const unsigned int d0)
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

double SeqKMeans::l2Dst(const unsigned int c, const unsigned int s, const unsigned int n)
{
    double w1 = 0.0f, w2 = 0.0f;
    int cloc = c * this->ndim;
    int sloc = s * this->ndim;
    for(unsigned int i = 0; i < this->ndim; i++)
    {
        w1  = this->arrayD[cloc+i]/n - this->data[sloc+i];
        w2 += w1*w1;
    }
    return w2;
}

double SeqKMeans::pairwDst(const int nclust, int seq[], int locount)
{
    unsigned int i = 0, label = 0, cloc = 0, c0 = 0, nl = 0, s = 0;
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
    for(i = 0; i < locount; i++)
    {
        s = seq[i];
        c0 = this->labels[s];
        nl = this->Ns[c0];
        l2nsum += nl * this->lens[s];
    }
    return (l2nsum-pd)/locount;
}

double SeqKMeans::l2nDst(double *D, const unsigned int x, const unsigned int n)
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

double SeqKMeans::crsDst(double Dsc, double *D, const unsigned int x, const unsigned int n)
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

double SeqKMeans::avgDist(const unsigned int nclust, int seq[], int locount, int lp)
{
    double *centers = NULL;
    double distort = 0;
    centers = new double[nclust * this->ndim];
    memset(centers, 0, sizeof(double)* nclust * this->ndim);

    unsigned int i = 0, j = 0, c = 0, dloc = 0, cloc = 0, s = 0;
    double avgDistort = 0;

    for(c = 0; c < nclust; c++)
    {
        cloc = this->ndim*c;
        for(j = 0; j < this->ndim; j++)
        {
            centers[cloc + j] = this->arrayD[cloc+j]/this->Ns[c];
        }
    }

    for(i = 0; i < locount; i++)
    {
        s = seq[i];
        c = this->labels[s];
        cloc = c*this->ndim;
        dloc = s*this->ndim;
        for(j = 0; j < this->ndim; j++)
        {
            distort = abs(centers[cloc + j] - this->data[dloc + j]);
            avgDistort += distort * distort;
        }
    }
    avgDistort /= locount;

    delete[] centers;
    centers = NULL;

    return avgDistort;
}

int SeqKMeans::clust(unsigned int clust_num, const char *dstfn, const int verbose)
{
    /*
    if(verbose)
        cout<<"Optimization function ............ I1"<<endl;
    */
    this->record_num = this->count-clust_num;
    if(this->record == NULL)
    {
        this->record = new double[this->record_num*6];
    }
    memset(this->record, 0, sizeof(double)*this->record_num*6);

    cout<<"Initialize record ................ ";
    for(unsigned i = 0; i < record_num; i++)
    {
        this->record[6*i]   = numeric_limits<double>::max();
        this->record[6*i+1] = -1;
        this->record[6*i+3] = numeric_limits<double>::max();
        this->record[6*i+4] = -1;
    }
    cout<<"done"<<endl;
    for(unsigned ntrail = 0; ntrail < nTrails; ntrail++)
    {
        cout<<"ntrail: "<<ntrail<<endl;
        opt1trail(clust_num);
    }
    for(unsigned i = clust_num; i < this->count; i++)
    {
        unsigned ri = i - clust_num;
        cout<<"avgI2: "<<i<<"\t"<<this->record[6*ri]<<"\t"<<this->record[6*ri+1]<<"\t"<<this->record[6*ri+2]/nTrails<<endl;
        cout<<"avgI1: "<<i<<"\t"<<this->record[6*ri+3]<<"\t"<<this->record[6*ri+4]<<"\t"<<this->record[6*ri+5]/nTrails<<endl;
    }
}

int SeqKMeans::opt1trail(unsigned int clust_num)
{
    unsigned int k, kt, r, q, i , j, d, loc, loc1, labelid, num, nj, s = 0, bestc, cloc, sloc, ri;
    vector<double *> center;
    vector<double > distance;
    double * tmp;
    int centerid, *seq = NULL;;
    double w, f, p, tmpdistance, mindistance, avgDistort = 0, avgDst = 0;

    if(this->Ns != NULL)
    {
        delete [] this->Ns;
        this->Ns = NULL;
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

    this->Ns      = new int[clust_num];
    this->arrayD = new double[clust_num*this->ndim];
    this->Ds      = new double[clust_num];
    memset(this->Ns, 0, sizeof(int)*clust_num);
    memset(this->arrayD, 0, sizeof(double)*clust_num*this->ndim);
    memset(this->Ds, 0, sizeof(double)*clust_num);

    this->getL2norms(this->count, this->ndim);

    cout<<"Initial centers ..................";
    seq = new int[this->count];
    for(i = 0; i < this->count; i++)
        seq[i] = i;
    random_shuffle(seq, seq+this->count);

    for(i = 0; i < clust_num; i++)
    {
        s = seq[i];
        cloc = i*this->ndim;
        sloc = s*this->ndim;
        for(j = 0; j < this->ndim; j++)
            this->arrayD[cloc + j] += this->data[sloc + j];
        this->Ns[i] += 1;
        this->labels[s] = i;
        this->Ds[i] += this->lens[s];
    }
    cout<<" done"<<endl;

    cout<<"Sequential data clustering ......."<<endl;
    for(i = clust_num; i < this->count; i++)
    {
        cout<<"data: "<<i;
        s = seq[i];
        mindistance = numeric_limits<double>::max();
        for(j = 0; j < clust_num; j++)
        {
            nj = this->Ns[j];
            cloc  = j * this->ndim;

            ///kmeans
            ///tmpdistance = l2Dst(j, s, nj);

            ///i2
            ///tmpdistance = crsDst(this->Ds[j], this->arrayD+cloc, s, nj);

            ///i1
            tmpdistance = l2nDst(arrayD+cloc, s, nj);
            tmpdistance = tmpdistance/(nj+1)/(nj+1);

            if(tmpdistance < mindistance)
            {
                bestc = j;
                mindistance = tmpdistance;
            }
        }
        this->Ns[bestc] += 1;
        this->labels[s] = bestc;
        this->Ds[bestc] += this->lens[s];
        cloc = bestc * this->ndim;
        sloc = s * this->ndim;
        for(j = 0; j < this->ndim; j++)
            this->arrayD[cloc+j] += this->data[s*ndim+j];
        avgDistort = avgDist(clust_num, seq, i+1, 2);
        avgDst     = pairwDst(clust_num, seq, i+1);
        cout<<"\t"<<avgDistort<<"\t"<<avgDst<<endl;

        ri = i - clust_num;
        if(this->record[6*ri] > avgDst)
            this->record[6*ri] = avgDst;
        if(this->record[6*ri+1] < avgDst)
            this->record[6*ri+1] = avgDst;
        this->record[6*ri+2] += avgDst;

        if(this->record[6*ri+3] > avgDistort)
            this->record[6*ri+3] = avgDistort;
        if(this->record[6*ri+4] < avgDistort)
            this->record[6*ri+4] = avgDistort;
        this->record[6*ri+5] += avgDistort;
    }

    /*
    k = ceil(double(kt - 15)/5);

    cout<<"kt\t"<<kt<<"\tk\t"<<k<<"\n";

    ///initialize center
    if(this->count < k + 10)
        cout<<"error\tdatanum=\t"<<this->count<<"\tfisrtronund clust nun=\t"<<k<<"\n";

    for(i = 0; i < k + 10; i++)
    {
        tmp = new double[this->ndim];
        loc = i*ndim;
        for(j = 0; j < ndim; j++)
        {
            tmp[j] = data[loc + j];
        }
        center.push_back(tmp);
        distance.push_back(RAND_MAX);
        labels[i] = i;
    }
    num = center.size();
    ///find w's value

    for(i = 0; i < k + 10; i++)
    {
        for(j = 0; j < k + 10; j++)
        {
            if(i == j)
                continue;
            tmpdistance = 0;
            for(d = 0; d < this->ndim; d++)
            {
                tmpdistance += pow(center[i][d] - center[j][d],2);
            }
            if(tmpdistance < distance[i])
                distance[i] = tmpdistance;
        }
    }

    sort(distance.begin(), distance.end());

    w = 0;
    for(i = 0; i < 10; i++)
    {
        w += distance[i];
    }

    q = 0;

    for(i = k + 10; i < this->count; i++)
    {
        ///probability p

        loc = i*ndim;
        for(j = 0; j < num; j++)
        {
            mindistance = RAND_MAX;
            tmpdistance = 0;
            for(d = 0; d < this->ndim; d++)
            {
                tmpdistance += pow(data[loc+d]- center[j][d],2);
            }
            if(tmpdistance < mindistance)
            {
                mindistance = tmpdistance;
                labelid = j;
            }
        }

        if(mindistance >= w || rand()*w < mindistance*RAND_MAX)
        {
            labels[i] = num;
            num++;
            tmp = new double[this->ndim];
            for(j = 0; j < ndim; j++)
            {
                tmp[j] = data[loc + j];
            }
            center.push_back(tmp);
            q++;
        }
        else
        {
            labels[i] = labelid;
        }
        if(q >= k)
        {
            q = 0;
            w *= 10;
        }
    }


    for(i = 0; i < num; i++)
    {
        this->infoMap[i].n = 0;
    }

    ///get clust ;
    float *tmpDat;
    tmpDat = new float[ndim];
    for(i = 0; i < count; i++)
    {
        mindistance = RAND_MAX;
        loc    = i*ndim;
        for(d = 0; d < ndim; d++)
        {
            tmpDat[d] = data[loc+d];
        }
        for(j = 0; j < num; j++)
        {
            tmpdistance = PQMath::l2d(center[j], 0, tmpDat, 0, ndim);
            if(tmpdistance < mindistance)
            {
                mindistance = tmpdistance;
                labelid     = j;
            }
        }
        labels[i] = labelid;
        this->infoMap[labelid].n = this->infoMap[labelid].n + 1;
        this->Ns[labelid] += 1;
        loc = labelid*ndim;
        for(j = 0; j < ndim; j++)
        {
            this->arrayD[loc+j] += tmpDat[j];
        }
    }
    ///save_clust(dstfn);

    this->calAVGDist(this->arrayD, num, infoMap);
    for(i = 0; i < num; i++)
    {
        loc = i*ndim;
        for(d = 0; d < ndim; d++)
        {
            this->arrayD[loc +d] = center[i][d];
        }

    }

    cout<<"clust num\t"<<clust_num<<"\n";

    delete [] tmpDat;

    for(i = 0; i < clust_num; i++)
    {
        delete [] center[i];
    }
    center.clear();
    distance.clear();
    */
    return 0;
}

void SeqKMeans::saveCenters(const char *dstfn, bool append)
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

int SeqKMeans::fetchCenters(float *centers)
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

void SeqKMeans::test()
{

    const char *srcMatFn1 = "/home/chenghaod/data/sift/raw/sift_base.txt";
    const char *srcMatFn2 = "/home/chenghaod/src/cpp/newxkmean/xkmean/15dat";
    ///const char *srcMatFn3 = "/home/wlzhao/sift_learn.txt";
    const char *srcMatFn3 = "/home/rqchen/dataset/sift_learn.txt";

    ///const char *dstFn = "test_olk.txt";
    const char *dstFn = "/home/rqchen/dataset/rslt/test_seqkmeans.txt";

    SeqKMeans *mykmean = new SeqKMeans();
    mykmean->buildcluster(srcMatFn3, dstFn, "non", "xx", "xx", 1024, false);
    delete mykmean;

}
