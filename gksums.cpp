#include "gksums.h"

#include <unordered_map>
#include <iostream>
#include <cassert>
#include <cstring>
#include <fstream>
#include <queue>
#include <time.h>

#include "vstring.h"
#include "ioagent.h"
#include "cleaner.h"
#include "pqmath.h"
#include "timer.h"

//#define _DB_

using namespace std;
const int GKSums::k0  = 40;
const int GKSums::BD0 = 100;
const int GKSums::nIter1 = 4;
const int GKSums::nIter2 = 4;    ///set interation times
const int GKSums::nTrail = 10;    ///set repeat times


GKSums::GKSums()
{
    this->_INIT_    = false;
    this->dstSum    = NULL;
    this->checked   = NULL;
    this->dstSum    = NULL;
    this->nCmps     = 0;
    this->tCmps     = 0;
    this->errSum    = 0;
    this->verbose   = 1;
    this->bstLabels = NULL;
    this->labels    = NULL;
    this->entropies = NULL;

    this->tcs1 = NULL;
    this->tcs2 = NULL;
}

bool GKSums::init(const char *srcfn)
{
    cout<<"Method ........................... gksum (graph-based)"<<endl;
    cout<<"Loading matrix ................... ";
    assert(srcfn);
    strcpy(srcmatfn, srcfn);
    this->refresh();

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

    this->checked = new char[this->count];
    memset(this->checked, 0, sizeof(char)*this->count);

    this->entropies = new float[this->count];
    memset(this->entropies, 1.0f, this->count*sizeof(float));

    this->knnGraph.resize(this->count);
    this->initKnnGraph(GKSums::k0);
    this->nbGraph.resize(this->count);

    this->_INIT_  = true;
    this->_REFER_ = false;
    return true;
}

bool GKSums::init(float *mat, const int row, const int dim)
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

    this->checked = new char[this->count];
    memset(this->checked, 0, sizeof(char)*this->count);

    this->entropies = new float[this->count];
    memset(this->entropies, 1.0f, this->count*sizeof(float));

    this->_INIT_  = true;
    this->_REFER_ = true;
    return true;
}

bool GKSums::config(const char *_seed_, const char *crtrn, const char *lg_first, int verbose)
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

bool GKSums::refresh()
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

    if(this->arrayD != NULL)
    {
        delete [] this->arrayD;
        this->arrayD = NULL;
    }

    if(this->checked != NULL)
    {
        delete [] this->checked;
        this->checked = NULL;
    }

    if(this->bstLabels != NULL)
    {
        delete [] this->bstLabels;
        this->bstLabels = NULL;
    }

    if(this->entropies != NULL)
    {
        delete [] this->entropies;
        this->entropies = NULL;
    }

    return true;
}

void GKSums::getRndSds(std::mt19937 &rng, const unsigned int k, unsigned int *addr, const unsigned int N)
{
    unsigned int i = 0;
    if (N == k)
    {
        for (i = 0; i < k; ++i)
        {
            addr[i] = i;
        }
        return;
    }
    for (i = 0; i < k; ++i)
    {
        addr[i] = rng() % (N - k);
    }
    sort(addr, addr + k);
    for (i = 1; i < k; ++i)
    {
        if (addr[i] <= addr[i-1])
        {
            addr[i] = addr[i-1] + 1;
        }
    }
    unsigned off = rng() % N;
    for (i = 0; i < k; ++i)
    {
        addr[i] = (addr[i] + off) % N;
    }

}

double GKSums::calAvgDistort(const unsigned int nclust)
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

int GKSums::loadkNNGraph(const char *graphFn)
{
    ifstream *inStrm = new ifstream(graphFn, ios::in);
    unsigned int row = 0, col = 0, iRow = 0, idx = 0, j = 0;
    if(!inStrm->is_open())
    {
        cout<<"Fail to open '"<<graphFn<<"' for read!\n";
        exit(0);
    }
    (*inStrm)>>row;
    (*inStrm)>>col;
    assert(row == count);
    this->knnGraph.resize(row);

    while(!inStrm->eof())
    {
        (*inStrm)>>idx;
        (*inStrm)>>col;
        this->knnGraph[iRow].resize(col);
        for(j = 0; j < col; j++)
        {
            MiniNN &itm = this->knnGraph[iRow][j];
            (*inStrm)>>itm.idx;
            (*inStrm)>>itm.val;
            //itm.val = 0;
        }
        //cout<<iRow<<"\t"<<this->knnGraph[iRow].size()<<endl;
        iRow++;
    }
    inStrm->close();
    //cout<<this->knnGraph.size();
    //exit(0);
    return 0;
}

int GKSums::saveKNNGraph(const char *dstFn)
{
    ofstream *outStrm = new ofstream(dstFn, ios::out);
    unsigned int i = 0, j =  0;
    (*outStrm)<<this->knnGraph.size()<<" "<<this->k0<<endl;
    for(i = 0; i < this->knnGraph.size(); i++)
    {
        (*outStrm)<<i<<" "<<k0;
        for(j = 0; j < k0; j++)
        {
            (*outStrm)<<" "<<knnGraph[i][j].idx;
        }
        (*outStrm)<<endl;
    }
    outStrm->close();
    return 0;
}

int GKSums::initKnnGraph(const unsigned int k)
{
    unsigned i = 0, j = 0, L = 0;
    unsigned *seeds  = new unsigned[512];
    std::mt19937 rng(time(NULL));
    cout<<"Initialize nn graph .............. ";
    for( i = 0; i < count; i++)
    {
        knnGraph[i].resize(k);
        this->getRndSds(rng, k+1, seeds, count);
        L = 0;
        for (j = 0; j < k+1;  j++)
        {
            MiniNN &nb = knnGraph[i][L];

            if(seeds[j] == i) continue;

            nb.idx  = seeds[j];
            assert(nb.idx != i);
            nb.val  = PQMath::l2f(data, i, data, nb.idx, this->ndim);
            nb.nw   = 1;
            this->nCmps++;
            L++;
            if(L >= k)
                break;
        }
        sort(knnGraph[i].begin(), knnGraph[i].end(), MiniNN::LLcomparer);
        ///cout<<i<<"\t"<<knnGraph[i].size()<<endl;
    }
    delete [] seeds;
    seeds = NULL;
    cout<<"done \n";

    return 0;
}

int GKSums::updateLst(const unsigned int id, const unsigned int nb, float dst)
{
    int i =0, topk =0, j = 0, l;
    i = topk = k0;
    vector<MiniNN> & addr = this->knnGraph[id];

    if(addr[topk-1].val < dst)
        return 0;

    while(i > 0)
    {
        j = i-1;
        if(addr[j].val <= dst) break;
        i = j;
    }
    l = i;
    while(l > 0)
    {
        j = l - 1;
        if(addr[j].val < dst) break;
        if(addr[j].idx == nb) return 0;
        l = j;
    }

    j = topk - 1;

    while(j > i)
    {
        addr[j].idx = addr[j-1].idx;
        addr[j].val = addr[j-1].val;
        addr[j].nw  = addr[j-1].nw;
        --j;
    }

    addr[j].idx = nb;
    addr[j].val = dst;
    addr[j].nw  = 1;

    return j;
}


int GKSums::getNBGraph(const unsigned int smplNum)
{
    float maxDst = 0;
    for(unsigned irow = 0; irow < this->nbGraph.size(); irow++)
    {
        auto &newlist  = this->nbGraph[irow].newnb;
        auto &oldlist  = this->nbGraph[irow].oldnb;

        unsigned nnSz = k0;
        for(unsigned idim = 0; idim < nnSz; idim++)
        {
            ///Neighborx &nb = knnGraph[irow].pool[idim];

            MiniNN &nb = this->knnGraph[irow][idim];

            auto &nhood_o = nbGraph[nb.idx];

            maxDst = knnGraph[nb.idx][nnSz - 1].val;
            if(nb.nw)
            {
                newlist.push_back(nb.idx);

                if( maxDst < nb.val )
                {
                    ///LockGuard guard(nhood_o.lock);
                    nhood_o.rnewnb.push_back(irow);
                }
                nb.nw = false;
            }
            else
            {
                oldlist.push_back(nb.idx);
                if(maxDst < nb.val)
                {
                    ///LockGuard guard(nhood_o.lock);
                    nhood_o.roldnb.push_back(irow);
                }
            }
        }///for(idim)
    }///for(irow)


    for( unsigned irow = 0; irow < this->nbGraph.size(); irow ++ )
    {
        vector<unsigned> &newlist  = nbGraph[irow].newnb;
        vector<unsigned> &oldlist  = nbGraph[irow].oldnb;
        vector<unsigned> &rnewlist = nbGraph[irow].rnewnb;
        vector<unsigned> &roldlist = nbGraph[irow].roldnb;

        random_shuffle(rnewlist.begin(), rnewlist.end());
        if(rnewlist.size() > smplNum)
        {
            rnewlist.resize(smplNum);
        }

        random_shuffle(roldlist.begin(), roldlist.end());

        if(roldlist.size() > smplNum)
        {
            roldlist.resize(smplNum);
        }

        newlist.insert(newlist.end(), rnewlist.begin(), rnewlist.end());
        oldlist.insert(oldlist.end(), roldlist.begin(), roldlist.end());
    }///for(irow)
    if(0)
    cout<<"\nDistance Sum is: "<<errSum<<"\t"<<this->nCmps<<endl;
    return 0;
}

unsigned GKSums::nnDescent()
{
    unsigned int i = 0, j = 0, k = 0, a = 0, b = 0;
    unsigned ccmps = 0, ignored = 0;
    float dst = 0;
    ///cout<<"my graph size: "<<this->nbGraph.size()<<endl;
    cout<<"perform NNDescent ...... \n";
    for(i = 0; i < this->nbGraph.size(); i++)
    {
        if(entropies[i] < this->entrp0)
         {
             ignored++;
             continue;
         }

        auto & oldnb = this->nbGraph[i].oldnb;
        auto & newnb = this->nbGraph[i].newnb;
        for(j = 0; j < oldnb.size(); j++)
        {
            a = oldnb[j];
            for(k = 0; k < newnb.size(); k++)
            {
                b = newnb[k];
                dst = PQMath::l2f(data, a, data, b, ndim);
                this->nCmps++;
                ccmps++;

                updateLst(a, b, dst);
                updateLst(b, a, dst);
            }
        }

        if(newnb.size() == 0)
            continue;

        for(j = 0; j < newnb.size()-1; j++)
        {
            a = newnb[j];
            for(k = j+1; k < newnb.size(); k++)
            {
                b = newnb[k];
                dst = PQMath::l2f(data, a, data, b, ndim);
                updateLst(a, b, dst);
                updateLst(b, a, dst);
                this->nCmps++;
                ccmps++;
            }
        }

    }///i

    cout<<"ignored:\t"<<ignored<<endl;

    return ccmps;
}

int GKSums::initNon(const unsigned int nclust)
{
    unsigned int i = 0, j = 0, c = 0, loc = 0;
    if(this->arrayD == NULL)
    {
       this->arrayD = new double[nclust*this->ndim];
    }

    if(this->labels == NULL)
    {
        this->labels = new int[this->count];
    }

    if(this->Ns == NULL)
    {
        this->Ns = new int[nclust];
    }
    memset(this->arrayD, 0, sizeof(double)*nclust*this->ndim);
    memset(this->labels, 0, sizeof(int)*this->count);
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
    }
    cout<<"done\n";
    return 0;
}


int GKSums::clust(const unsigned int nclust, const char *dstFn, const int verbose)
{
    double avgDist = 0, minDist = RAND_MAX + 0.0f;
    unsigned int i = 0, t = 0, ccmps = 0;
    float rate = 0.0f;
    if(nclust >= k0)
    {
        this->entrp0 = log(40)/log(2);
    }else{
        this->entrp0 = log(nclust)/log(2);
    }
    for(i = 0; i < this->count; i++)
    {
        this->entropies[i] = this->entrp0;
    }
    this->entrp0 = this->entrp0*0.5;

    this->dstSum = new double[nclust];
    memset(this->dstSum, 0, sizeof(double)*nclust);

    this->distortions1 = new double[nIter1+1];
    memset(this->distortions1, 0, sizeof(double)*(nIter1+1));
    this->distortions2 = new double[nIter2+1];
    memset(this->distortions2, 0, sizeof(double)*(nIter2+1));
    this->tcs1 = new double[nIter1+1];
    memset(this->tcs1, 0, sizeof(double)*(nIter1+1));
    this->tcs2 = new double[nIter2+1];
    memset(this->tcs2, 0, sizeof(double)*(nIter2+1));

    this->bg = clock();
    cout<<"cluster numb ..................... "<<nclust<<endl;
    cout<<"Optimize ......................... \n";
    ///this->tCmps += this->nCmps;
    i = 0;
    this->nCmps = 0;

    for(i = 0; i < 5; i++)
    {
        getNBGraph(100);
        ccmps = nnDescent();
        Cleaner::clearNbs(this->nbGraph);
    }
    ///cout<<"i am here\n";
    this->initNon(nclust);
    optzIn(nIter1, nclust);
    for(i = 0; i < 8; i++)
    {
        optzI2(nIter2, nclust);
        getNBGraph(100);
        ccmps = nnDescent();
        Cleaner::clearNbs(this->nbGraph);
    }

    avgDist = calAvgDistort(nclust);

    Cleaner::clearkNNGraph(knnGraph);
    cout<<"Run " << t << "\tDistortion: "<<avgDist<<endl;

    this->save_clust(dstFn);

    cout<<"\n\t\tdone\n";

    rate = (this->nCmps*2.0)/(this->count*(this->count-1.0));
    cout<<"Averge distortion ................ "<<minDist<<endl;
    cout<<"Scanning rate .................... "<<rate<<endl;

    return 0;
}


int GKSums::numInKnn(const unsigned int i, const int clabel)
{
    vector<MiniNN>::iterator lit;
    int num = 0;

    vector<MiniNN> & crntList = this->knnGraph[i];
    for(lit = crntList.begin(); lit != crntList.end(); lit++)
    {
        MiniNN & itm = *lit;
        if(this->labels[itm.idx] == clabel)
            num++;
    }

    return num;
}

double GKSums::calDist(double *D, const unsigned int x, const unsigned int n)
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

double GKSums::cosDst(double *D, const unsigned int x, bool _plus_)
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

int GKSums::optzIn(const unsigned int nIter, const int nclust)
{
    unsigned int i = 0, iter = 0, cloc = 0, dloc = 0, sloc = 0;
    double rvsNumA = 0, rvsNumB = 0, rvsErr = 0, rvsMax = 0;
    unsigned int s = 0, c = 0, c0 = 0, nwC = 0, j = 0;
    int numA = 0, numB = 0, maxNum = 0;
    double tcost = 0, avgDistort = 0;

    set<unsigned int>::iterator sit;
    vector<MiniNN>::iterator vit;
    set<unsigned int> others;
    int *seq = NULL;
    seq = new int[this->count];

    for(i = 0; i < this->count; i++)
    {
        seq[i] = i;
    }

    this->ed = clock();
    tcost = (double)(this->ed - this->bg)/(0.0 + CLOCKS_PER_SEC);
    this->tcs1[0] = tcost;

    avgDistort = calAvgDistort(nclust);
    this->distortions1[0] = avgDistort;
    cout<<"Start Integral clustering ........... \n";
    double sum = 0;
    for(iter = 0; iter < nIter; iter++)///set nIter to 30 normally
    {
        random_shuffle(seq, seq+this->count);
        sum  = 0;
        for(i = 0; i < this->count; i++)///one round
        {
            s = seq[i];
            if(this->checked[s] == 1)
                continue;

            c0 = this->labels[s]; ///only one element, no movement
            if(this->Ns[c0] == 1)
                continue;

            vector<MiniNN> &crntLst = knnGraph[s];
            this->checked[s] = 1;
            for(vit = crntLst.begin(); vit != crntLst.end(); vit++)
            {
                MiniNN &crntItm = *vit;
                c  = this->labels[crntItm.idx];

                if(c != c0)
                {
                   others.insert(c);
                }
            }

            numA    = numInKnn(s, c0);
            rvsNumA = (double)numA / Ns[c0];
            rvsMax  = 0;
            ///search for the best movement
            for(sit = others.begin(); sit != others.end(); sit++)
            {
                numB    = numInKnn(s, *sit);
                rvsNumB = (double)numB / Ns[*sit];
                rvsErr  = rvsNumB - rvsNumA;

                if(rvsErr > rvsMax)
                {
                    rvsMax = rvsErr;
                    nwC    = *sit;
                }
            }

            sum += maxNum;
            if(rvsMax > 0) ///carry out the movement
            {
                this->labels[s] = nwC;
                this->Ns[nwC]  += 1;
                this->Ns[c0]   -= 1;
                cloc = c0*ndim;
                dloc = nwC*ndim;
                sloc = s*ndim;
                for(j = 0; j < this->ndim; j++)
                {
                   this->arrayD[cloc+j] -= this->data[sloc + j];
                   this->arrayD[dloc+j] += this->data[sloc + j];
                }
            }
            others.clear();
        }///end-for(i)

        this->ed = clock();
        tcost = (double)(this->ed - this->bg)/CLOCKS_PER_SEC;
        ///this->tcs1[iter+1] = tcost;

        avgDistort = calAvgDistort(nclust);
        ///this->distortions1[iter+1] = avgDistort;

        if(1)
        {
            cout<<nclust<<"\t"<<"iter"<<iter+1<<": "<<avgDistort<<"\t\ttime cost: "<<tcost<<"\tsec(s)"<<endl;
        }
        memset(this->checked, 0, this->count);
    }///end-for(iter)

    delete [] seq;
    seq = NULL;
    cout<<"end"<<endl;

    return 0;

}

int GKSums::optzI2(const unsigned int nIter, const int nclust)
{
    unsigned int i = 0, iter = 0, j = 0, cloc = 0, sloc = 0, dloc = 0;
    unsigned int s = 0, c = 0, c0 = 0, nwC = 0, n1 = 0, n2 = 0;
    double distA = 0, distB = 0, err = 0, minErr = 0;

    unordered_map<int, short>::iterator mit;
    unordered_map<int, short> groups;
    set<unsigned int>::iterator sit;
    vector<MiniNN>::iterator vit;
    set<unsigned int> nblabels;
    double avgDistort = 0;
    float entrp = 0, szknn = 0, a = 0, b = 1.0/log(2);
    clock_t tcost = 0;
    int label = 0;

    int *seq = NULL;
    memset(this->Ns, 0, sizeof(int)*nclust);
    seq = new int[this->count];
    for(i = 0; i < this->count; i++)
    {
        c = this->labels[i];
        this->Ns[c] += 1;
        seq[i] = i;
    }

    avgDistort = calAvgDistort(nclust);
    for(iter = 0; iter < nIter; iter++)///set nIter to 30 normally
    {
        random_shuffle(seq, seq+this->count);

        for(i = 0; i < this->count; i++)///one round
        {
            s = seq[i];
            if(this->checked[s] == 1)
                continue;

            c0 = this->labels[s]; ///only one element, no movement
            if(this->Ns[c0] == 1)
                continue;

            this->checked[s] = 1;

            cloc = c0*ndim;
            n1 = this->Ns[c0]-1;

            /**/
            distA = calDist(this->arrayD+cloc, s, n1+1);
            distA = distA/(n1*n1);
            ///distA = cosDst(this->arrayD+cloc, s, false);
            szknn = knnGraph[s].size();
            for(vit = knnGraph[s].begin(); vit != knnGraph[s].end(); vit++)
            {
                MiniNN &crntNN = *vit;
                label = this->labels[crntNN.idx];
                if(groups.find(label) == groups.end())
                {
                    groups.insert(pair<int, short>(label, 1));
                }else{
                    groups[label] = groups[label] + 1;
                }
                if(label != c0)
                nblabels.insert(label);
            }

            minErr = 0;
            ///search for the best movement
            for(sit = nblabels.begin(); sit != nblabels.end(); sit++)
            {
                label = *sit;
                if(label != (int)c0)
                {
                    cloc  = label*ndim;
                    n2    = this->Ns[label] + 1;
                    /**/
                    distB = calDist(arrayD+cloc, s, n2-1);
                    distB = distB/(n2*n2);/**/
                    ///distB = cosDst(this->arrayD+cloc, s, true);

                    err = distB - distA;
                    if(err < minErr)
                    {
                        minErr = err;
                        nwC    = label;
                    }
                }
            }
            entrp = 0.0f;
            for(mit = groups.begin(); mit != groups.end(); mit++)
            {
                 a = mit->second/szknn;
                 entrp -= a*log(a)*b;
                 ///cout<<mit->first<<"\t"<<mit->second<<"\t"<<szknn<<"\t"<<b<<endl;
            }
            entropies[s] = entrp;

            nblabels.clear();
            groups.clear();

            if(minErr >= 0.0f)
            {
                continue;
            }

            ///carry out movement
             this->labels[s] = nwC;
             this->Ns[nwC]  += 1;
             this->Ns[c0]   -= 1;
             cloc = c0*ndim;
             dloc = nwC*ndim;
             sloc = s*ndim;
             for(j = 0; j < this->ndim; j++)
             {
                 this->arrayD[cloc+j] -= this->data[sloc + j];
                 this->arrayD[dloc+j] += this->data[sloc + j];
             }

        }///end-for(i)
        //exit(0);

        this->ed = clock();
        tcost = (double)(this->ed - this->bg)/(0.0 + CLOCKS_PER_SEC);

        avgDistort = calAvgDistort(nclust);

        if(1)
          cout<<"iter"<<iter+1<<": "<<avgDistort<<"\t\ttime cost: "<<tcost<<"\tsec(s)"<<endl;

        memset(this->checked, 0, this->count);
    }///end-for(iter)

    delete [] seq;
    seq = NULL;

    return avgDistort;
}


void GKSums::saveCenters(const char *dstfn, bool append)
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

int GKSums::fetchCenters(float *centers)
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
            centers[idxi + j] = this->arrayD[clabel*this->ndim+j]/this->Ns[clabel];
        }
        idxi += this->ndim;
    }
    return rCNum;
}

void GKSums::saveDistort(const char* dstFn)
{
    ofstream *outStrm = new ofstream(dstFn, ios::out);
    unsigned int i = 0;

    for(i = 0; i <= this->nIter1; i++)
    {
        (*outStrm)<<i<<"\t"<<this->distortions1[i]<<"\t\t"<<this->tcs1[i]<<endl;
    }
    (*outStrm)<<"\n";
    for(i = 1; i <= this->nIter2; i++)
    {
         (*outStrm)<<i<<"\t"<<this->distortions2[i]<<"\t\t"<<this->tcs2[i]<<endl;
    }

    outStrm->close();
}

GKSums::~GKSums()
{
    if(this->data != NULL)
    {
        delete [] this->data;
        this->data = NULL;
    }

    if(this->checked != NULL)
    {
        delete [] this->checked;
        this->checked = NULL;
    }

    if(this->labels != NULL)
    {
        delete [] this->labels;
        this->labels = NULL;
    }

    if(this->bstLabels != NULL)
    {
        delete [] this->bstLabels;
        this->bstLabels = NULL;
    }

    if(this->knnGraph.size() > 0)
    {
        Cleaner::clearkNNGraph(this->knnGraph);
    }

    if(this->dstSum != NULL)
    {
        delete [] dstSum;
        this->dstSum = NULL;
    }


    if(this->distortions1 != NULL)
    {
        delete [] this->distortions1;
        this->distortions1 = NULL;
    }

    if(this->distortions2 != NULL)
    {
        delete [] this->distortions2;
        this->distortions2 = NULL;
    }

    if(this->entropies != NULL)
    {
        delete [] this->entropies;
        this->entropies = NULL;
    }

    if(this->tcs1 != NULL)
    {
         delete [] this->tcs1;
         this->tcs1 = NULL;
    }

    if(this->tcs2 != NULL)
    {
         delete [] this->tcs2;
         this->tcs2 = NULL;
    }
}

void GKSums::test()
{

    const char *srcFn1 = "/home/hye/datasets/clust/sift100k/sift_learn.txt";
    const char *knng1 = "/home/hye/datasets/clust/sift100k/sift100k_40knn.txt";
    const char *distortFn1 = "/home/hye/result/sift100k_10+60_tcs_distort.txt";

    const char *srcFn2 = "/home/hye/datasets/clust/sift1m/sift_base.txt";
    const char *knng2 = "/home/hye/datasets/clust/sift1m/sift1m_40knn.txt";
    const char *distortFn2 = "/home/hye/result/sift1m_ksum_distort.txt";

    const char *srcFn3 = "/home/hye/datasets/glove1m/glove1m_base.txt";
    const char *knng3 = "/home/hye/datasets/glove1m/glove1m_40knn.txt";
    const char *distortFn3 = "/home/hye/result/glove1m_10+60_tcs_distort.txt";

    const char *srcFn4 = "/home/hye/datasets/clust/susy/SUSY.csv";
    const char *knng4 = "";
    const char *distortFn4 = "/home/hye/datasets/clust/susy/SUSY_60t_distort.txt";

    const char *srcFn5 = "/home/hye/datasets/clust/msd/msd-rh.txt";
    const char *knng5 = "/home/hye/datasets/clust/msd/msd-rh_40knn.txt";
    const char *distortFn5 = "/home/hye/datasets/clust/msd/msd-rh_ksums_distort.txt";


    unsigned int cnum = 10000;

    GKSums *mykm = new GKSums();
    Timer *mytm = new Timer();
    mytm->start();
    mykm->init(srcFn5);

    mykm->clust(cnum,"", 0);
    mytm->end(true);

    mykm->saveDistort(distortFn5);
}

void GKSums::wltest()
{
    const char *srcFn1 = "/home/wlzhao/datasets/bignn/sift1m/sift_learn.txt";
    const char *knng1  = "/home/wlzhao/datasets/bignn/sift1m/sift100k_k=40v17.txt";
    const char *distortFn1 = "/home/wlzhao/datasets/clust/distort/sift100k_80_gksum_distort_v4.txt";

    const char *srcFn2 = "/home/wlzhao/datasets/bignn/sift1m/sift_base.txt";
    const char *distortFn2 = "/home/wlzhao/datasets/clust/distort/sift_80_gksum_distort.txt";
    const char *dstFn2 = "/home/wlzhao/datasets/clust/rslt/sift_clust.txt";

    const char *srcFn3 = "/home/wlzhao/datasets/bignn/glove/glove1m_base.txt";
    const char *distortFn3 = "/home/wlzhao/datasets/clust/distort/glove_80_gksum_distort.txt";
    const char *dstFn3 = "/home/wlzhao/datasets/clust/rslt/glove_clust.txt";

    const char *srcFn4 = "/home/wlzhao/datasets/clust/msd-rh.txt";
    const char *distortFn4 = "/home/wlzhao/datasets/clust/distort/msd-rh_80_gksum_distort.txt";
    const char *dstFn4 = "/home/wlzhao/datasets/clust/rslt/msd_clust.txt";

    const char *srcFn5 = "/home/wlzhao/datasets/clust/SUSY.csv";
    const char *distortFn5 = "/home/wlzhao/datasets/clust/distort/susy_80_gksum_distort.txt";
    const char *dstFn5 = "/home/wlzhao/datasets/clust/rslt/susy_clust.txt";

    unsigned int cnum = 10000;

    GKSums *mykm = new GKSums();
    Timer *mytm = new Timer();
    mytm->start();
    mykm->init(srcFn3);
    mykm->clust(cnum, "", 0);
    mykm->saveDistort(distortFn3);
    mytm->end(true);

}
