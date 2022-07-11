#include "abstractkmeans.h"
#include "randseed.h"
#include "pqmath.h"
#include "timer.h"

#include <iostream>
#include <iomanip>
#include <cassert>
#include <fstream>
#include <cstring>
#include <cmath>

#include "pqmath.h"


using namespace std;

const unsigned int AbstractKMeans::paraBound = 1024;
const float AbstractKMeans::smallVal0 = 0.0000001f;

AbstractKMeans::AbstractKMeans()
{
    srand(time(NULL));
    this->ndim     = 0;
    this->labels   = NULL;
    this->data     = NULL;
    this->infoMap  = NULL;
    this->dataMap  = NULL;
    this->Ns       = NULL;
    this->lens     = NULL;
    this->arrayD   = NULL;
    this->Es       = NULL;
    this->dataType = 0;
    this->seed     = _kpp_;
    this->sdata.col = NULL;
    this->sdata.index = NULL;
    this->sdata.data = NULL;
    this->SHUFFLE  = 1;
    this->_INIT_   = false;
    this->_REFINE_ = false;
    this->dstort   = 0;
    this->nLogs      = 0;
    this->kmLogs[0]  = NULL;
    this->kmLogs[1]  = NULL;
    //  this->nThrd   = HardWareSetup::getCoreNum()/4;
    this->allEg    = 0;
    strcpy(srcmatfn, "");
    strcpy(dataFn,   "");
}


bool AbstractKMeans::initMemry(const unsigned int dNum, const unsigned int clustNum)
{
    unsigned int i = 0;
    if(this->lens == NULL)
    {
        this->lens   = new double[this->count];
    }
    if(this->labels == NULL)
    {
        this->labels = new int[this->count];
    }
    memset(this->labels, 0, sizeof(int)*this->count);
    memset(this->lens,   0, sizeof(double)*this->count);

    if(this->count > this->clnumb)
    {
        this->infoMap = new CLSInfo[clustNum];
        this->dataMap = new CLData[clustNum];

        if(this->Es == NULL)
        this->Es      = new double[clustNum];
        memset(this->Es, 0, sizeof(double)*clustNum);

        if(this->Ns == NULL)
        this->Ns      = new int[clustNum];
        memset(this->Ns, 0, sizeof(int)*clustNum);

        if(this->arrayD == NULL)
        this->arrayD = new double[clustNum*this->ndim];
        memset(this->arrayD, 0, sizeof(double)*clustNum*this->ndim);
    }else{
        for(i = 0; i < this->count; i++)
        {
            this->labels[i] = i;
        }
    }

    if(kmMtd == _rbkmn_ )
    {
        if(!dataType)
        {
            this->normVects(this->data, this->ndim, this->count, this->lens);
        }
        else
            this->normVects(this->sdata, this->ndim, this->count, this->lens);
        for(i = 0; i < this->count; i++)
        {
            if(lens[i] == 0)
            {
                ///    labels[i] = -1;
            }
            else
            {
                labels[i] = 0;
            }
        }
        cout<<"Normalize input vectors .......... ";
        cout<<"yes\n";
        allEg = this->count;
    }
    else
    {
        if(!dataType)
            allEg = this->normVects(this->data, this->lens, this->ndim, this->count);
        else
            this->normVects(this->sdata, this->lens, this->ndim, this->count);
        cout<<"Normalize input vectors .......... ";
        cout<<"no\n";
    }

    return true;
}

bool AbstractKMeans::relsMemry()
{
    if(this->infoMap != NULL)
    {
        delete [] this->infoMap;
        this->infoMap = NULL;
    }

    if(this->dataMap != NULL)
    {
        delete [] this->dataMap;
        this->dataMap = NULL;
    }

    if(this->lens != NULL)
    {
        delete [] this->lens;
        this->lens = NULL;
    }

    if(this->labels != NULL)
    {
        delete [] this->labels;
        this->labels = NULL;
    }

    if(this->Es != NULL)
    {
        delete [] this->Es;
        this->Es = NULL;
    }

    if(this->arrayD != NULL)
    {
       delete [] this->arrayD;
       this->arrayD = NULL;
    }

    return true;
}

bool AbstractKMeans::setLogOn(const unsigned logSize)
{
    assert(logSize > 0);

    this->nLogs     = logSize;
    for(unsigned int i = 0; i < 128; i++)
    {
        this->kmLogs[i] = new double[logSize];
    }
    return 1;
}

void AbstractKMeans::saveLogs(const char *dstFn, const unsigned nIt)
{
    unsigned int i = 0, numb = this->nLogs > nIt?nIt:this->nLogs;

    if(this->nLogs <= 0)
    {
        return ;
    }

    ofstream *outStrm = new ofstream(dstFn, ios::out);
    for(i = 0 ; i < numb; i++)
    {
        (*outStrm)<<i<<"\t\t"<<(allEg - kmLogs[0][i])/this->count<<"\t\t"<<kmLogs[1][i]<<endl;
    }

    outStrm->close();
}

unsigned int AbstractKMeans::buildcluster(const char *srcfn, const char *dstfn, const char *_seed_, const char *lg_first,
        const char *crtrn, const int num, bool _refine_)
{
    ///cout<<srcfn<<endl;
    this->config(_seed_, lg_first, crtrn, 1);

    this->clnumb   = num;
    clock_t t1 = 0;
    double t2 = 0;

    this->_INIT_   = this->init(srcfn);

    this->_REFINE_ = _refine_;
    char logfn[512];
    Timer *mytm = new Timer();
    mytm->start();

    relsMemry();

    if(!this->_INIT_ || this->count == 0)
    {
        this->clnumb = 0;
    }
    else if(this->count <= this->clnumb)
    {
        this->infoMap = new CLSInfo[this->count];
        this->dataMap = new CLData[this->count];
        this->labels  = new int[this->count];
        this->clnumb  = this->nvclust(this->count, dstfn, 1);
    }
    else
    {
        this->initMemry(this->count, clnumb);

        this->clnumb  = this->clust(clnumb, dstfn, 1);

        t2 = (clock() - t1+0.0f)/CLOCKS_PER_SEC;
        sprintf(logfn, "%s_tkmeans_ks_log.txt", dstfn);
        ofstream *outStrm = new ofstream(logfn, ios::app);
        (*outStrm)<<this->count<<"\t"<<clnumb<<"\t"<<this->dstort<<"\t"<<t2<<endl;
        cout<<this->count<<"\t"<<clnumb<<"\t"<<this->dstort<<"\t"<<t2<<endl;

        outStrm->close();
    }
    ///sprintf(msgStr, "%s\t%s\t%s\t%s\t%d\t%lf\t", mthStr, _seed_, lg_first, crtrn, num, this->dstort);
    mytm->end(true);
    //mytm->end(msgStr, "./recordnw.txt");
    delete mytm;
    return this->clnumb;
}

double AbstractKMeans::refineSparse(const unsigned int clustNum, const unsigned int NRef)
{
    unsigned int numb   = this->count, loci = 0;
    unsigned int i = 0, t = 0, j = 0, k = 0, k0 = 0;
    int loc = 0;

    double tmpEg1 = 0, optEg = 0, tmpEg2 = 0;
    unsigned int iter = 0, r = 0, tk = 0;
    double delta = 0, delta1 = 0, mxEs = 0;
    unsigned int * rl = new unsigned int[numb];
    for(i = 0; i < clustNum; i++)
    {
        loci  = i*ndim;
        Es[i] = PQMath::dvec_norm(arrayD+loci, ndim, 2);

        /**     if(myoptz  == _I2_)
             {
                 if(Ns[i] > 1)
                 {
                     Es[i] = Es[i]*(Es[i]/Ns[i]);
                 }
                 optEg += Es[i];
             }
             else if(myoptz  == _I4_)
          /**/
        {
            optEg += Es[i];
            Es[i] *= Es[i];
        }
    }
    ///cout<<"iter 0: "<<optEg<<endl;
    for(i = 0; i < numb; i++)
    {
        rl[i] = i;
    }
    for(t = 0; t < NRef; t++)
    {
        random_shuffle(rl, rl + numb);
        for(iter = 0; iter < numb; iter++)
        {
            if(!SHUFFLE)
            {
                r   = (unsigned int)floor((rand()/(RAND_MAX+1.0f))*numb);
                r   = r>=numb?(numb-1):r;
            }
            else
                r = rl[iter];
            k0  = labels[r];

            if(Ns[k0] <= 1)
                continue;

            loci   = k0*ndim;
            if(myoptz  == _I2_)
            {


            }
            //  else if(myoptz  == _I4_)
            {
                tmpEg1 = 0;
                for(j = sdata.col[r]; j < sdata.col[r+1]; j++)
                {
                    tmpEg1 += arrayD[loci + sdata.index[j]] * sdata.data[j];
                }
                tmpEg1 = Es[k0] - 2*tmpEg1 + this->lens[r];
            }
            delta1 = 0;
            for(k = 0; k < clustNum; k++)
            {
                if(k == k0)
                    continue;

                loci   = k*ndim;
                if(myoptz  == _I2_)
                {

                }
                //   else if(myoptz  == _I4_)
                {
                    tmpEg2 = 0;
                    for(j = sdata.col[r]; j < sdata.col[r+1]; j++)
                    {
                        tmpEg2 += arrayD[loci + sdata.index[j]] * sdata.data[j];
                    }
                    tmpEg2 = Es[k] + 2*tmpEg2 + this->lens[r];
                }
                delta  = - sqrt(Es[k]) + sqrt(tmpEg2) + sqrt(tmpEg1) - sqrt(Es[k0]);

                if(delta > delta1)
                {
                    delta1 = delta;
                    tk = k;
                    mxEs = tmpEg2;
                }
            }

            if(delta1 > 0)
            {
                Ns[k0]--;
                Ns[tk]++;
                labels[r] = tk;
                Es[k0] = tmpEg1;
                Es[tk] = mxEs;
                loc = k0*ndim;
                loci = tk*ndim;
                for(j = sdata.col[r]; j < sdata.col[r+1]; j++)
                {
                    this->arrayD[loc+sdata.index[j]] -= sdata.data[j];
                    this->arrayD[loci+sdata.index[j]] += sdata.data[j];
                }
                optEg = optEg+delta1;
            }
        }///for(iter)
    }  ///for(t)

    delete [] rl;
    rl = NULL;
    return optEg;
}

unsigned int AbstractKMeans::buildcluster(float *mat, const int row, const int dim, const char *dstfn,
        const char *_seed_, const char *lg_first, const char *crtrn, const int num, bool _refine_)
{
    assert(mat);
    assert(crtrn);
    assert(lg_first);

    this->config(_seed_, lg_first, crtrn, 0);
    this->clnumb = num;
    this->_INIT_ = this->init(mat, row, dim);
    this->_REFINE_ = _refine_;

    Timer *mytm = new Timer();
    mytm->start();

    relsMemry();

    if(!this->_INIT_ || this->count == 0)
    {
        this->clnumb = 0;
    }
    else if(this->count <= this->clnumb)
    {
        this->infoMap = new CLSInfo[clnumb];
        this->dataMap = new CLData[clnumb];
        this->clnumb  = this->nvclust(clnumb, dstfn, 1);
    }
    else
    {
        this->initMemry(this->count, clnumb);
        this->clnumb  = this->clust(clnumb, dstfn, 0);
    }

    mytm->end(true);

    delete mytm;

    return this->clnumb;
}

unsigned int AbstractKMeans::nvclust(const unsigned int clnumb, const char *dstfn,  const int verbose)
{
    unsigned int i = 0, j = 0, loc = 0;
    for(i = 0; i < clnumb; i++)
    {
        this->infoMap[i].n = 1;
        this->infoMap[i].E = 1.0f;
    }

    for(i = 0; i < this->count; i++)
    {
        this->labels[i] = i;
        for(j = 0; j < this->ndim; j++)
        {
            arrayD[loc+j] = data[loc+j];
        }
        loc += ndim;
    }

    return this->count;
}

float *AbstractKMeans::getCluster(const unsigned int clabel0, unsigned int &row, unsigned int &dim)
{
    vector<int>::iterator it;
    vector<int> vects;
    unsigned int i = 0;
    int loc1 = 0, loc2 = 0, v = 0;
    unsigned int cpysize = sizeof(float)*this->ndim;

    for(i = 0; i < this->count; i++)
    {
        if((int)clabel0 == this->labels[i])
        {
            vects.push_back(i);
        }
    }
    row = vects.size();
    float *cl_data = new float[this->infoMap[clabel0].n*this->ndim];

    dim = this->ndim;
    assert(row == this->infoMap[clabel0].n);
    for(it = vects.begin(), i = 0; it != vects.end(); it++)
    {
        v = *it;
        loc1 = v*this->ndim;
        loc2 = i*this->ndim;
        memcpy(cl_data+loc2, this->data+loc1, cpysize);
        i++;
    }
    vects.clear();

    return cl_data;
}

double AbstractKMeans::refine(const unsigned int clustNum, const unsigned int NRef)
{
    unsigned int loc = 0, i = 0, j = 0, k = 0, t = 0, k0 = 0;
    unsigned int numb   = this->count, loci = 0;

    double tmpEg1 = 0, optEg = 0, tmpEg2 = 0;
    unsigned int iter = 0, r = 0, tk = 0;
    double delta = 0, delta1 = 0, mxEs = 0;
    unsigned int * rl = new unsigned int[numb];
    for(i = 0; i < clustNum; i++)
    {
        loci  = i*ndim;
        Es[i] = PQMath::dvec_norm(arrayD+loci, ndim, 2);

        if(myoptz  == _I2_)
        {
            if(Ns[i] > 1)
            {
                Es[i] = Es[i]*(Es[i]/Ns[i]);
            }
            optEg += Es[i];
        }
        else if(myoptz  == _I4_)
        {
            optEg += Es[i];
            Es[i] *= Es[i];
        }
    }
    ///cout<<"iter 0: "<<optEg<<endl;
    for(i = 0; i < numb; i++)
    {
        rl[i] = i;
    }
    for(t = 0; t < NRef; t++)
    {
        random_shuffle(rl, rl + numb);
        for(iter = 0; iter < numb; iter++)
        {
            if(!SHUFFLE)
            {
                r   = (unsigned int)floor((rand()/(RAND_MAX+1.0f))*numb);
                r   = r>=numb?(numb-1):r;
            }
            else
                r = rl[iter];
            k0  = labels[r];
            loc = r*this->ndim;

            if(Ns[k0] <= 1)
                continue;

            loci   = k0*ndim;
            if(myoptz  == _I2_)
                tmpEg1 = I2FastM(arrayD+loci, data+loc, ndim, Ns[k0]);
            else if(myoptz  == _I4_)
            {
                tmpEg1 = 0;
                for(i = 0; i < ndim; i++)
                {
                    tmpEg1 += arrayD[loci + i]* data[loc+i];
                }
                tmpEg1 = Es[k0] - 2*tmpEg1 + lens[r];
                tmpEg1 = sqrt(tmpEg1);
            }
            delta1 = 0;
            for(k = 0; k < clustNum; k++)
            {
                if(k == k0)
                    continue;

                loci   = k*ndim;
                if(myoptz  == _I2_)
                    tmpEg2 = I2FastP(arrayD+loci, data+loc, ndim, Ns[k]);
                else if(myoptz  == _I4_)
                {
                    tmpEg2 = 0;
                    for(i = 0; i < ndim; i++)
                    {
                        tmpEg2 += arrayD[loci + i] * data[loc+i];
                    }
                    tmpEg2 = Es[k] + 2*tmpEg2 +this->lens[r];
                    tmpEg2 = sqrt(tmpEg2);
                }
                delta  = tmpEg2 - sqrt(Es[k]) + tmpEg1 - sqrt(Es[k0]);

                if(delta > delta1)
                {
                    delta1 = delta;
                    tk = k;
                    mxEs = tmpEg2;
                }
            }

            if(delta1 > 0)
            {
                Ns[k0]--;
                Ns[tk]++;
                labels[r] = tk;
                Es[k0] = tmpEg1*tmpEg1;
                Es[tk] = mxEs*mxEs;
                for(j = 0; j < this->ndim; j++)
                {
                    this->arrayD[k0*ndim+j] -= data[loc+j];
                    this->arrayD[tk*ndim+j] += data[loc+j];
                }
                optEg = optEg+delta1;
                //cout<<"i a here\n";
            }
        }///for(iter)
    }  ///for(t)

    delete [] rl;
    rl = NULL;
    return optEg;
}

double AbstractKMeans::refine(int label1, int label2, int label3, const unsigned int NRef)
{
    vector<unsigned int> &cluster1= this->dataMap[label1].clustDataId;
    vector<unsigned int> &cluster2= this->dataMap[label2].clustDataId;
    vector<unsigned int> &cluster3= this->dataMap[label3].clustDataId;
    vector<int > clusterLabels;
    clusterLabels.push_back(label1);
    clusterLabels.push_back(label2);
    clusterLabels.push_back(label3);

    unsigned int loc = 0, i = 0, j = 0, k = 0, t = 0, k0 = 0;
    unsigned int numb   = cluster1.size() + cluster2.size() + cluster3.size();
    unsigned int loci = 0;

    double tmpEg1 = 0, optEg = 0, tmpEg2 = 0;
    unsigned int iter = 0, r = 0, tk = 0;
    double delta = 0, delta1 = 0, mxEs = 0;
    unsigned int * rl = new unsigned int[numb];
    for(i = 0; i < cluster1.size(); i++)
    {
        rl[i] = cluster1[i];
    }
    j = 0;
    for(i = cluster1.size(); i < cluster1.size() + cluster2.size(); i++,j++)
    {
        rl[i] = cluster2[j];
    }
    j = 0;
    for(i = cluster1.size() + cluster2.size(); i < cluster1.size() + cluster2.size() + cluster3.size(); i++, j++)
    {
        rl[i] = cluster3[j];
    }

    for(j = 0; j < 3; j++)
    {
        i = clusterLabels[j];
        loci  = i*ndim;
        Es[i] = PQMath::dvec_norm(arrayD+loci, ndim, 2);

        if(myoptz  == _I2_)
        {
            if(Ns[i] > 1)
            {
                Es[i] = Es[i]*(Es[i]/Ns[i]);
            }
            optEg += Es[i];
        }
        else if(myoptz  == _I4_)
        {
            optEg += Es[i];
            Es[i] *= Es[i];
        }
    }
    ///cout<<"iter 0: "<<optEg<<endl;

    for(t = 0; t < NRef; t++)
    {
        random_shuffle(rl, rl + numb);
        for(iter = 0; iter < numb; iter++)
        {
            if(!SHUFFLE)
            {
                r   = (unsigned int)floor((rand()/(RAND_MAX+1.0f))*numb);
                r   = r>=numb?(numb-1):r;
            }
            else
                r = rl[iter];
            k0  = labels[r];
            loc = r*this->ndim;

            if(Ns[k0] <= 1)
                continue;

            loci   = k0*ndim;
            if(myoptz  == _I2_)
                tmpEg1 = I2FastM(arrayD+loci, data+loc, ndim, Ns[k0]);
            else if(myoptz  == _I4_)
            {
                tmpEg1 = 0;
                for(i = 0; i < ndim; i++)
                {
                    tmpEg1 += arrayD[loci + i]* data[loc+i];
                }
                tmpEg1 = Es[k0] - 2*tmpEg1 + lens[r];
                tmpEg1 = sqrt(tmpEg1);
            }
            delta1 = 0;
            for(i = 0; i < 3; i++)
            {
                k = clusterLabels[i];
                if(k == k0)
                    continue;

                loci   = k*ndim;
                if(myoptz  == _I2_)
                    tmpEg2 = I2FastP(arrayD+loci, data+loc, ndim, Ns[k]);
                else if(myoptz  == _I4_)
                {
                    tmpEg2 = 0;
                    for(i = 0; i < ndim; i++)
                    {
                        tmpEg2 += arrayD[loci + i] * data[loc+i];
                    }
                    tmpEg2 = Es[k] + 2*tmpEg2 +this->lens[r];
                    tmpEg2 = sqrt(tmpEg2);
                }
                //delta  = tmpEg2 - sqrt(Es[k]) + tmpEg1 - sqrt(Es[k0]);
                delta = tmpEg2 - Es[k] + tmpEg1 - Es[k0];
                if(delta > delta1)
                {
                    delta1 = delta;
                    tk = k;
                    mxEs = tmpEg2;
                    break;
                }
            }

            if(delta1 > 0)
            {
                Ns[k0]--;
                Ns[tk]++;
                labels[r] = tk;
                Es[k0] = tmpEg1;//*tmpEg1;
                Es[tk] = mxEs;//*mxEs;
                for(j = 0; j < this->ndim; j++)
                {
                    this->arrayD[k0*ndim+j] -= data[loc+j];
                    this->arrayD[tk*ndim+j] += data[loc+j];
                }
                optEg = optEg+delta1;
                //cout<<"i a here\n";
            }
        }///for(iter)
    }  ///for(t)

    this->dataMap[clusterLabels[0]].clustDataId.clear();
    this->dataMap[clusterLabels[1]].clustDataId.clear();
    this->dataMap[clusterLabels[2]].clustDataId.clear();

    for(i = 0; i < numb; i++)
    {

        if(labels[rl[i]] == clusterLabels[0])
            this->dataMap[clusterLabels[0]].clustDataId.push_back(rl[i]);
        else if(labels[rl[i]] == clusterLabels[1])
            this->dataMap[clusterLabels[1]].clustDataId.push_back(rl[i]);
        else if(labels[rl[i]] == clusterLabels[2])
            this->dataMap[clusterLabels[2]].clustDataId.push_back(rl[i]);
        else
        {
            cout<<"cluster 3 refine erro \n";
            exit(0);
        }
    }

    delete [] rl;
    rl = NULL;
    return optEg;
}

double  AbstractKMeans::r2refine(const unsigned int clustNum, const unsigned int NRef)
{
    unsigned int i = 0, j = 0;
    int *cs = new int[clustNum];
    int c1 = 0, c2 = 0;

    for(i = 0; i < clustNum; i++)
    {
        cs[i] = i;
    }

    for(i = 0; i < NRef; i++)
    {
        random_shuffle(cs, cs+clustNum);
        for(j = 0; j < this->clnumb; j++)
        {
            c1 = rand()%clnumb;
            c2 = (c1+1)%clnumb;
            c1 = cs[c1];
            c2 = cs[c2];
            if(dataType)
                SBiOptz(c1, c2, arrayD+c1*ndim, arrayD+c2*ndim, Ns[c1], Ns[c2]);
            else
                BiOptz(c1, c2, arrayD+c1*ndim, arrayD+c2*ndim, Ns[c1], Ns[c2]);
        }
    }
    delete [] cs;
    cs = NULL;
    return 0;
}

bool AbstractKMeans::BiOptz(int clbl1, int clbl2, double *D1, double *D2, int &n1, int &n2)
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
        }
        else if(this->labels[i] == clbl2)
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
        if(iter > 64)
            break;
    }
    while(UPDATED);

    for(i = 0; i < num; i++)
    {
        if(mylabels[i] == 0)
        {
            j = crntClust[i];
            this->labels[j] = clbl1;
        }
        else
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

bool AbstractKMeans::SBiOptz(int clbl1, int clbl2, double *D1, double *D2, int &n1, int &n2)
{
    unsigned int i = 0, j = 0, loc = 0, iter = 0;
    double optEg = 0, tmpEg1 = 0, tmpEg2 = 0, delta = 0, ds1 = 0, ds2 = 0;
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
        }
        else if(this->labels[i] == clbl2)
        {
            crntClust.push_back(i);
            mylabels[j] = 1;
            j++;
        }
    }
    //cout<<num<<"\t"<<j<<endl;

    optEg = 0;
    ds1 = 0;
    ds2 = 0;
    for(i = 0; i < ndim; i++)
    {
        ds1 += D1[i]*D1[i];
        ds2 += D2[i]*D2[i];
    }


    if(n1 == 0)
        optEg += 0;
    else
        optEg += ds1/n1;

    if(n2 == 0)
        optEg += 0;
    else
        optEg += ds2/n2;

    do
    {
        UPDATED = false;
        for(i = 0; i < num ; i++)
        {
            loc = crntClust[i];

            tmpEg1 = 0;
            tmpEg2 = 0;
            for(j = sdata.col[loc]; j < sdata.col[loc+1]; j++)
            {
                tmpEg1 += D1[sdata.index[j]]*sdata.data[j];
                tmpEg2 += D2[sdata.index[j]]*sdata.data[j];
            }

            if(mylabels[i] == 0)
            {
                tmpEg1 = ds1 - tmpEg1*2 + lens[loc];
                tmpEg2 = ds2 + tmpEg2*2 + lens[loc];
                if(n1 == 1)
                    delta = tmpEg2/(n2+1) - optEg;
                else
                    delta = tmpEg2/(n2+1) + tmpEg1/(n1-1) - optEg;
            }
            else if(mylabels[i] == 1)
            {
                tmpEg1 = ds1 + tmpEg1*2 + lens[loc];
                tmpEg2 = ds2 - tmpEg2*2 + lens[loc];
                if(n2 == 1)
                    delta = tmpEg1/(n1+1) - optEg;
                else
                    delta = tmpEg2/(n2-1) + tmpEg1/(n1+1) - optEg;
            }


            if(delta > 0)
            {
                optEg += delta;
                ds1 = tmpEg1;
                ds2 = tmpEg2;
                if(mylabels[i] == 0)
                {
                    n1--;
                    n2++;
                    for(j = sdata.col[loc]; j < sdata.col[loc+1]; j++)
                    {
                        D1[sdata.index[j]] -= sdata.data[j];
                        D2[sdata.index[j]] += sdata.data[j];
                    }
                    mylabels[i] = 1;
                }
                else if(mylabels[i] == 1)
                {
                    n1++;
                    n2--;
                    for(j = sdata.col[loc]; j < sdata.col[loc+1]; j++)
                    {
                        D1[sdata.index[j]] += sdata.data[j];
                        D2[sdata.index[j]] -= sdata.data[j];
                    }
                    mylabels[i] = 0;
                }
                UPDATED = true;
            }//if(delta)
        }
        iter++;
        if(iter > 64)
            break;
    }
    while(UPDATED);

    for(i = 0; i < num; i++)
    {
        if(mylabels[i] == 0)
        {
            j = crntClust[i];
            this->labels[j] = clbl1;
        }
        else
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


double AbstractKMeans::I2FastM(const double *D1, const float *v1, const unsigned int dim, const unsigned int n0)
{
    double sumDst = 0, v = 0;
    unsigned int i = 0, n = 0;
    if(n0 == 1)
        return 0;

    n = n0 - 1;
    for(i = 0; i < dim; i++)
    {
        v = D1[i]-v1[i];
        sumDst += v*v;
    }
    sumDst = sumDst/n;

    return sumDst;
}

double AbstractKMeans::I2FastP(const double *D1, const float *v1, const unsigned int dim, const unsigned int n0)
{
    double sumDst = 0, v = 0;
    unsigned int i = 0, n = 0;
    n = n0 + 1;
    for(i = 0; i < dim; i++)
    {
        v = D1[i]+v1[i];
        sumDst += v*v;
    }
    sumDst = sumDst/n;

    return sumDst;
}

double AbstractKMeans::getI2(const double *Ds, const unsigned int k, const unsigned int dim, const int *tmpns)
{
    double sumDst = 0, tmpE = 0;
    unsigned int i = 0, j = 0, loc = 0;

    for(j = 0; j < k; j++)
    {
        if(tmpns[j] <= 0)
            continue;

        loc = dim*j;
        tmpE = 0;

        for(i = 0; i < dim; i++)
        {
            tmpE += Ds[loc+i]*Ds[loc+i];
        }
        sumDst += tmpE/tmpns[j];
    }
    return sumDst;
}


double AbstractKMeans::getI4(const double *Ds, const unsigned int k, const unsigned int dim, const int *tmpns)
{
    double sumDst = 0, tmpE = 0;
    unsigned int i = 0, j = 0, loc = 0;

    for(j = 0; j < k; j++)
    {
        if(tmpns[j] <= 0)
            continue;

        loc = dim*j;
        tmpE = 0;

        for(i = 0; i < dim; i++)
        {
            tmpE += Ds[loc+i]*Ds[loc+i];
        }
        sumDst += sqrt(tmpE);
    }
    return sumDst;
}

double AbstractKMeans::getE2(const double *Ds, const unsigned int k, const unsigned int dim, const int *tmpns)
{
    double sumDst = 0, tmpE = 0;
    unsigned int i = 0, j = 0, loc = 0, loc1 = 0;

    for(j = 0; j < k; j++)
    {
        if(tmpns[j] <= 0)
            continue;

        loc1 = dim*j;
        for(unsigned int j2 = j + 1; j2 < k; j2++)
        {
            tmpE = 0;
            if(tmpns[j2] <= 0)
                continue;
            loc = dim*j2;
            for(i = 0; i < dim; i++)
            {
                tmpE += (Ds[loc+i]/tmpns[j2] -Ds[loc1+i]/tmpns[j])*(Ds[loc+i]*tmpns[j] -Ds[loc1+i]*tmpns[j2]);
            }
            sumDst += tmpE*tmpns[j]*tmpns[j2];
        }

    }
    return sumDst;
}

double AbstractKMeans::getI1(const double *Ds, const unsigned int k, const unsigned int dim, const int *tmpns)
{
    double i1 = 0, tmpE;
    unsigned int i, j, loc = 0;
    int clustNum = this->clnumb;
    double *clustE = new double[clustNum];

    memset(clustE, 0, sizeof(double)*clustNum);
    for(i = 0; i < this->count; i++)
    {
        clustE[this->labels[i]] += this->lens[i];
    }

    for(j = 0; j < k; j++)
    {
        if(tmpns[j] <= 0)
            continue;

        loc = dim*j;
        tmpE = 0;

        for(i = 0; i < dim; i++)
        {
            tmpE += Ds[loc+i]*Ds[loc+i];
        }
        if(tmpns[j] == 1)
        {
            tmpE = 0;
            clustE[j] = 0;
        }
        else
        {
            tmpE = tmpE/(tmpns[j] - 1);
            clustE[j] = clustE[j]*tmpns[j]/(tmpns[j] - 1);
        }
        i1 += tmpE - clustE[j];
    }
    delete [] clustE;
    return i1;
}

double AbstractKMeans::getE1(const double *Ds, const unsigned int k, const unsigned int dim, const int *tmpns)
{
    double e1 = 0, tmpE = 0;
    unsigned int i, j, loc = 0, loc1 = 0;
    int clustNum = this->clnumb;
    double *clustE = new double[clustNum];

    memset(clustE, 0, sizeof(double)*clustNum);
    for(i = 0; i < this->count; i++)
    {
        clustE[this->labels[i]] += this->lens[i];
    }

    for(j = 0; j < k; j++)
    {
        if(tmpns[j] <= 0)
            continue;

        loc1 = dim*j;
        for(unsigned int j2 = j + 1; j2 < k; j2++)
        {
            tmpE = 0;
            if(tmpns[j2] <= 0)
                continue;
            loc = dim*j2;
            for(i = 0; i < dim; i++)
            {
                tmpE += Ds[loc+i]*Ds[loc1+i];
            }
            e1 += (clustE[j]*tmpns[j] + clustE[j2]*tmpns[j2] - 2*tmpE)*tmpns[j]*tmpns[j2];
        }

    }
    delete [] clustE;
    return e1;
}

double AbstractKMeans::calAVGDist(const double *Ds, const unsigned int clustNum, CLSInfo *infos)
{
    double *dcts = new double[clustNum], tmpE = 0, sumEg = 0, tmpEg = 0;
    double sumDcts = 0;
    double *C    = new double[this->ndim];
    memset(dcts, 0, sizeof(double)*clustNum);
    unsigned int i, j = 0, loc = 0;
    int c = 0;
    for(i = 0; i < this->count; i++)
    {
        c   = labels[i];
        if(kmMtd == _rbkmn_)
        {
            dcts[c] += lens[i]*lens[i];
            sumDcts += lens[i]*lens[i];
        }
        else
        {
            dcts[c] += lens[i];
            sumDcts += lens[i];
        }
    }
    for(i = 0; i < clustNum; i++)
    {
        infos[i].n = this->Ns[i];

        if(this->Ns[i] < 1)
        {
            infos[i].E = 0;
            continue;
        }
        loc = i*this->ndim;
        tmpEg = 0;
        for(j = 0; j < this->ndim; j++)
        {
            C[j] = Ds[loc+j]/this->Ns[i];
            tmpEg += Ds[loc+j]*Ds[loc+j];
        }
        sumEg += tmpEg/this->Ns[i];
        tmpE   = PQMath::dvec_norm(C, this->ndim, 2);
        tmpE   = dcts[i] - infos[i].n*tmpE*tmpE;
        if(infos[i].n > 1)
        {
            infos[i].E = 2*(tmpE/(infos[i].n - 1));
        }
        else
        {
            infos[i].E = 0;
        }
    }

    this->dstort = (sumDcts - sumEg)/this->count;
    ///cout<<"E: "<<sumDcts<<"\t"<<sumEg<<"\t"<<"Distortion: "<<this->dstort<<endl;
    delete [] dcts;
    delete [] C;
    dcts = C = NULL;
    return sumEg;
}

void AbstractKMeans::printClusters(const char *dstdir)
{
    unsigned int clabel = 0, i,  j, loc, counter = this->count;
    CLSInfo &crntinfo    = infoMap[0];
    FILE *fp = NULL;
    char matfn[1024];
    assert(dstdir);

    for(j = 0; j < this->clnumb; j++)
    {
        crntinfo = infoMap[j];
        sprintf(matfn, "%s%d.mat", dstdir, j);
        fp = fopen(matfn, "w");
        fprintf(fp,"%d %d\n", crntinfo.n, this->ndim);
        fclose(fp);
    }

    for(i = 0; i < this->count; i++)
    {
        clabel = labels[i];
        loc = i*this->ndim;

        sprintf(matfn, "%s%d.mat", dstdir, clabel);
        fp = fopen(matfn, "a");
        for(j = 0; j < this->ndim; j++)
        {
            fprintf(fp,"%f ", data[loc+j]);
        }
        fprintf(fp,"\n");
        fclose(fp);
        counter--;
        cout<<"\r\r"<<setw(8)<<counter;
    }
    cout<<endl;
}

void AbstractKMeans::printCluster(const char *dstfn)
{
    FILE *fp = fopen(dstfn, "w");
    if(fp == NULL)
    {
        cout<<"File '"<<dstfn<<"' cannot open for write!\n";
        return;
    }
    unsigned int i;

    for(i = 0; i < this->count; i++)
    {
        fprintf(fp, "%d\n", this->labels[i]);
    }

    fclose(fp);
    return ;
}

int AbstractKMeans::idx2clabel(const int i)
{
    if(this->labels == NULL)
        return 0;

    if(i < 0 ||i >= (int)this->count || this->clnumb == 0)
    {
        return 0;
    }
    else
    {
        return this->labels[i];
    }
}

double AbstractKMeans::normVects(float *vects, double *lens,
                                 const unsigned int d0, const unsigned int n0)
{
    assert(vects);
    assert(lens);
    unsigned int i = 0, j = 0, loc = 0;
    double len = 0, sumEg1 = 0;
    for(i = 0; i < n0; i++)
    {
        len = 0;
        for(j = 0; j < d0; j++)
        {
            len += vects[loc+j]*vects[loc+j];
        }
        lens[i] = len;
        sumEg1  += lens[i];
        loc += d0;
    }

    return sumEg1;
}

void AbstractKMeans::normVects(const SparseMatrix &svects, double *lens,
                               const unsigned int d0, const unsigned int n0)
{
    assert(svects.data);
    assert(svects.index);
    assert(svects.col);
    assert(lens);
    unsigned int i = 0, j = 0;
    double len = 0;
    for(i = 0; i < n0; i++)
    {
        len = 0;
        for(j = svects.col[i]; j < svects.col[i+1]; j++)
        {
            len += svects.data[j]*svects.data[j];
        }
        lens[i] = len;
    }
}


void AbstractKMeans::normVects(SparseMatrix &svects, const unsigned int d0,
                               const unsigned int n0, double *lens)
{
    assert(svects.data);
    assert(svects.index);
    assert(svects.col);
    assert(lens);
    unsigned int i = 0, j = 0, loc = 0;
    float len = 0;
    for(i = 0; i < n0; i++)
    {
        len = 0;
        for(j = svects.col[i]; j < svects.col[i+1]; j++)
        {
            len += svects.data[j]*svects.data[j];
        }
        lens[i] = len;
        len     = sqrt(len);
        lens[i] = 1;
        if(len > AbstractKMeans::smallVal0)
            for(j = svects.col[i]; j < svects.col[i+1]; j++)
            {
                svects.data[j] = svects.data[j]/len;
            }
        else
        {
            lens[i] = 0;
        }

        loc += d0;
    }
    return ;
}

void AbstractKMeans::normVects(float *vects, const unsigned int d0,
                               const unsigned int n0, double *lens)
{
    assert(vects);
    assert(lens);
    unsigned int i = 0, j = 0, loc = 0;
    float len = 0;
    for(i = 0; i < n0; i++)
    {
        len = 0;
        for(j = 0; j < d0; j++)
        {
            len += vects[loc+j]*vects[loc+j];
        }
        lens[i] = len;
        len     = sqrt(len);
        lens[i] = 1;
        if(len > AbstractKMeans::smallVal0)
        {
            for(j = 0; j < d0; j++)
            {
                vects[loc+j] = vects[loc+j]/len;
            }
        }
        else
        {
            lens[i] = 0;
        }

        loc += d0;
    }

    return ;
}


double AbstractKMeans::calAvgDistort(const unsigned int nclust, int lp)
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

    if(lp == 2)
    {
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
    }
    else if(lp == 1)
    {
        for(i = 0; i < this->count; i++)
        {
            c = this->labels[i];
            cloc = c*this->ndim;
            dloc = i*this->ndim;
            for(j = 0; j < this->ndim; j++)
            {
                distort = abs(centers[cloc + j] - this->data[dloc + j]);
                avgDistort += distort;
            }
        }
    }else if(lp == 3)
    {
        for(i = 0; i < this->count; i++)
        {
            c = this->labels[i];
            cloc = c*this->ndim;
            dloc = i*this->ndim;
            for(j = 0; j < this->ndim; j++)
            {
                distort = abs(centers[cloc + j] - this->data[dloc + j]);
                distort = pow(distort, 2.2);
                avgDistort += distort;
            }
        }
    }
    else
    {
        cout<<lp << ": lp-norm undefined!!\n";
    }
    ///caculate avg distortion
    //cout<<"before: "<<avgDistort<<endl;
    avgDistort /= this->count;
    //cout<<"hi: "<<avgDistort<<"\t"<<this->count<<endl;

    delete[] centers;
    centers = NULL;

    return avgDistort;
}

void AbstractKMeans::saveCenters(const char *dstfn, bool append, bool _norm_)
{
    unsigned int clabel = 0, j, loc;
    ofstream *outStrm  = NULL;
    if(append)
    {
        outStrm = new ofstream(dstfn, ios::app);
    }
    else
    {
        outStrm = new ofstream(dstfn, ios::out);
        (*outStrm)<<this->clnumb<<" "<<this->ndim<<endl;;
    }

    for(clabel = 0; clabel < this->clnumb; clabel++)
    {
        loc  = clabel*this->ndim;
        for(j = 0; j < this->ndim; j++)
        {
            this->arrayD[loc+j] = this->arrayD[loc+j]/this->infoMap[clabel].n;
            ///cout<<this->arrayD[loc+j]<<" ";
        }

        /**/
        if(_norm_)
        {
            PQMath::l2_norm(this->arrayD+loc, ndim);
        }
        /**/

        for(j = 0; j < this->ndim; j++)
        {
            (*outStrm)<<this->arrayD[loc+j]<<" ";
        }
        (*outStrm)<<endl;
    }

    outStrm->close();

    return ;
}

bool AbstractKMeans::save_clust(const char *dstfn)
{
    if(dstfn == NULL || !strcmp(dstfn,""))
        return false;

    ofstream *outStrm = new ofstream(dstfn, ios::out);
    if(!outStrm->is_open())
    {
        cout<<"\nDestine file '"<<dstfn<<"' cannot open for write!\n";
        return false;
    }

    for(unsigned int i = 0; i < this->count; i++)
    {
        (*outStrm)<<this->labels[i]<<endl;
    }
    outStrm->close();
    return true;
}


void AbstractKMeans::resetlabels(unsigned int &n)
{
    unsigned int i;
    vector<unsigned > realid;
    map<int, unsigned> lbmap;
    unsigned int *lb = new unsigned[this->clnumb];
    memset(lb, 0, sizeof(unsigned)*(this->clnumb));
    for(i = 0; i < this->count; i++)
    {
        lb[labels[i]] = 1;
    }
    for(i = 0; i < this->clnumb; i++)
    {
        if(lb[i] != 0)
            realid.push_back(i);
    }
    for(i = 0; i < realid.size(); i++)
    {
        lbmap.insert(pair<unsigned, unsigned>(realid[i], i));
    }
    for(i = 0; i < this->count; i++)
    {
        labels[i] = lbmap[labels[i]];
    }
    n = lbmap.size();

    delete [] lb;
    lb = NULL;
}

AbstractKMeans::~AbstractKMeans()
{
    if(this->_REFER_)
    {
        this->data = NULL;
        this->sdata.data = NULL;
        this->sdata.index = NULL;
        this->sdata.col = NULL;
    }
    else
    {
        if(this->data != NULL)
        {
            delete [] this->data;
            this->data = NULL;
        }
        if(this->sdata.index != NULL)
        {
            delete [] this->sdata.index;
            delete [] this->sdata.col;
            this->sdata.index = NULL;
            this->sdata.col = NULL;

        }
        if(this->sdata.data != NULL)
        {
            delete [] this->sdata.data;
            this->sdata.data = NULL;
        }
    }
    relsMemry();

    if(this->nLogs > 0)
    {
        delete [] this->kmLogs[0];
        delete [] this->kmLogs[1];
        this->kmLogs[0] = NULL;
        this->kmLogs[1] = NULL;
        this->nLogs     = 0;
    }

}
