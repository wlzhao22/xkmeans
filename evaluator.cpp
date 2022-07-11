#include "evaluator.h"
#include "cleaner.h"
#include "pqmath.h"
#include "ioagent.h"
#include "vstring.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <limits>


using namespace std;

double Evaluator::distortion(const char* srcfn, const char* clusterfn)
{
    float *data;
    unsigned row, col;
    if(VString::endWith(srcfn, ".txt"))
    {
        data = IOAgent::loadMatrix(srcfn, row, col);
    }
    else if(VString::endWith(srcfn, ".fvecs"))
    {
        data = IOAgent::load_fvecs(srcfn, col, row);
    }
    else if(VString::endWith(srcfn, ".mat"))
    {
        data = IOAgent::loadDat(srcfn, row, col);
    }
    else if(VString::endWith(srcfn, ".itm"))
    {
        data = IOAgent::loadItms(srcfn, "fsvtab", row, col);
    }
    else
    {
        cout<<"Unrecognizable input file format!!!\n";
        data = NULL;
        exit(0);
    }

    map<unsigned int, set<int> *> clusters;
    double en;
    unsigned int i, j, k;
    en = 0;

    for(i = 0; i < row; i++)
    {
        for(j = 0; j < col; j++)
        en += data[i*col+j]*data[i*col+j];
    }
    Evaluator::loadClust(clusterfn, clusters);
    double opt;
    unsigned clnumb = clusters.size();
    double *Ds = new double[clnumb*col];
    memset(Ds, 0, sizeof(double)*clnumb*col);
    opt = 0;

    for(i = 0; i < clusters.size(); i++)
    {
        for(auto it = clusters[i]->begin(); it != clusters[i]->end(); it++)
        {
            for(k = 0; k < col; k++)
            {
                Ds[i*col+k] += data[(*it)*col+k];
            }
        }
    }
    for(i = 0; i < clusters.size(); i++)
    {
        for(j = 0; j < col; j++)
        {
            if(clusters[i]->size() != 0)
            opt += Ds[i*col+j]*Ds[i*col+j]/clusters[i]->size();
        }
    }
    double distortion;
    distortion = en - opt;
    distortion /= row;
    cout<<"en "<<en<<" opt "<<opt<<" distortion "<<distortion<<endl;
    delete [] Ds;
    delete [] data;
    Cleaner::clearClust(clusters);
    return distortion;
}

void Evaluator::jaccard(map<unsigned int, set<int> *> &clust, map<unsigned int, set<int> *> &grdTrh, const char *fn)
{
    unsigned int m = 0, n = 0, i = 0, j = 0;
    m = clust.size();
    n = grdTrh.size();

    float *val = new float[m*n];

    memset(val, 0, sizeof(float)*m*n);

    float tmp;
    vector<int> interset(numeric_limits<int>::max());

    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            if(clust[i+1]->size() != 0)
            {
                auto it = std::set_intersection(clust[i+1]->begin(), clust[i+1]->end(), grdTrh[j+1]->begin(), grdTrh[j+1]->end(),interset.begin());
                tmp = it - interset.begin();
                tmp = tmp/grdTrh[j+1]->size();
            }
            else
                tmp = 0;
            val[i*n+j] = tmp;
        }
    }
    ofstream os(fn);

    if(!os.is_open())
    {
        cout<<"can not open file "<<fn<<" for jaccard result\n";
        exit(0);
    }

    os<<"\t";
    for(i = 0; i < n; i++)
    {
        os<<setw(9)<<grdTrh[i+1]->size()<<"\t";
    }
    os<<"\n";
    for(i = 0; i < m; i++)
    {
        os<<clust[i+1]->size()<<"\t";
        for(j = 0; j < n; j++)
        {
            os<<setw(9)<<val[i*n+j]<<"\t";
        }
        os<<endl;
    }

    os.close();

    delete [] val;
}

void Evaluator::jaccard(const char *setfn, const char* grdfn, const char *dstfn)
{
    map<unsigned int, set<int> *> clust, grdTrh;

    IOAgent::loadClust(setfn, clust);
    IOAgent::loadClust(grdfn, grdTrh);

    Evaluator::jaccard(clust, grdTrh, dstfn);

    Cleaner::clearClust(clust);
    Cleaner::clearClust(grdTrh);
}

void Evaluator::loadClust(const char *srcFn, map<unsigned int, set<int> *> &clusts)
{
    string tmp;
    vector<string> clust;
    set<string> temp;
    set<string>::iterator it;
    map<string, unsigned int> mp;
    int i;

    if(clusts.size())
    {
        cout<<"clear clusts *******"<<endl;
        Cleaner::clearClust(clusts);
    }

    ifstream inStr(srcFn);

    if(!inStr.is_open())
    {
        cout<<"Evaluator::loadClust open file "<<srcFn<<" fail"<<endl;
        exit(0);
    }
    while(inStr>>tmp)
    {
        clust.push_back(tmp);
    }

    inStr.close();


    int s = clust.size();
    for(int i = 0; i < s; i++)
    {
        temp.insert(clust[i]);
    }

    i = 0;
    for(it = temp.begin(); it != temp.end(); it++)
    {
        set<int > * tmpclust = new set<int>;
        mp.insert(pair<string,unsigned int>(*it, i));
        clusts.insert(pair<unsigned int, set<int> *>(i, tmpclust));
        i++;

    }

    for(int i = 0; i < s; i++)
    {
        clusts[mp[clust[i]]]->insert(i);
    }

    clust.clear();

    return ;
}

float Evaluator::entropy(map<unsigned int, set<int> *> &clusts, map<unsigned int, set<int> *> &grdTrth)
{
    float E = 0;
    float tmpE = 0, tmp = 0;
    int q = grdTrth.size();
    float lgq = log(q);
    float *nr = new float[q];
    int n = 0;
    float nrs;
    for(unsigned int i = 0; i < clusts.size(); i++)
    {
        n += clusts[i]->size();
    }
    for(unsigned int j = 0; j < clusts.size(); j++)
    {
        for(int i = 0; i < q; i++)
        {
            nr[i] = Evaluator::interSect(*clusts[j], *grdTrth[i]);
        }
        nrs = clusts[j]->size();
        tmpE = 0;
        for(int i = 0; i < q; i++)
        {
            if(nr[i] != 0)
            {
                tmp = nr[i]/nrs;
                tmpE += tmp*log(tmp);
            }
        }
        E += -(tmpE/lgq)*(nrs/n);
    }
    delete [] nr;
    return E;
}

float Evaluator::entropy(set<int> &clust, map<unsigned int, set<int> *> &grdTrth)
{
    float E = 0;
    float tmpE = 0, tmp = 0;
    int q = grdTrth.size();
    float lgq = log(q);
    float *nr = new float[q];
    float nrs;

    for(int i = 0; i < q; i++)
    {
        nr[i] = Evaluator::interSect(clust, *grdTrth[i]);
    }
    nrs = clust.size();
    tmpE = 0;
    for(int i = 0; i < q; i++)
    {
        if(nr[i] != 0)
        {
            tmp = nr[i]/nrs;
            tmpE += tmp*log(tmp);
        }
    }
    E = -tmpE/lgq;
    delete [] nr;
    return E;
}

float Evaluator::purity(map<unsigned int, set<int> *> &clusts, map<unsigned int, set<int> *> &grdTrth)
{
    float P = 0;
    float tmp = 0;
    int q = grdTrth.size();
    float *nr = new float[q];
    int n = 0;
    for(unsigned int i = 0; i < clusts.size(); i++)
    {
        n += clusts[i]->size();
    }
    for(unsigned int j = 0; j < clusts.size(); j++)
    {
        for(int i = 0; i < q; i++)
        {
            nr[i] = Evaluator::interSect(*clusts[j], *grdTrth[i]);
        }
        tmp = 0;
        for(int i = 0; i < q; i++)
        {
            if(nr[i] > tmp)
                tmp = nr[i];
        }
        P += tmp;
    }
    P = P/n;
    delete [] nr;
    return P;
}

float Evaluator::purity(set<int> &clust, map<unsigned int, set<int> *> &grdTrth)
{
    float P = 0;
    float tmp = 0;
    int q = grdTrth.size();
    int n = clust.size();
    float *nr = new float[q];
    for(int i = 0; i < q; i++)
    {
        nr[i] = Evaluator::interSect(clust, *grdTrth[i]);
    }
    tmp = 0;
    for(int i = 0; i < q; i++)
    {
        if(nr[i] > tmp)
            tmp = nr[i];
    }

    P = tmp/n;
    delete [] nr;
    return P;
}

int Evaluator::interSect(const set<int> &s1, const set<int> &s2)
{
    int n = 0;
    set<int>::iterator it;
    for(it = s1.begin(); it != s1.end(); it++)
    {
        if(s2.find(*it) != s2.end())
            n++;
    }
    return n;
}

float Evaluator::funScore(const char *clustFn, const char *datFn, const char *fun)
{
    float score = 0;
    double *Ds;
    unsigned int k;
    unsigned int row;
    unsigned int dim;
    int *tmpns;

//load clusts;
    map<unsigned int, set<int> *>clusts;
    float *dat;
    Evaluator::loadClust(clustFn, clusts);
    dat = IOAgent::loadDat(datFn, row, dim);
    k = clusts.size();
    tmpns = new int[k];
    Ds = new double[k*dim];
//caculate Ds;
    for(unsigned int i = 0; i < k; i++)
    {
        tmpns[i] = clusts[i]->size();
        int loc = i*dim;

        for(unsigned int j = 0; j < dim; j++)
        {
            Ds[loc + j] = 0;

            for(set<int>::iterator k = clusts[i]->begin(); k != clusts[i]->end(); k++)
            {
                Ds[loc + j] += dat[j + (*k)*dim];
            }
        }

    }
    if(!strcmp(fun, "i2"))
    {
        score = PQMath::getI2(Ds, k, dim, tmpns);
    }
    else if(!strcmp(fun, "e2"))
    {
        score = PQMath::getE2(Ds, k, dim, tmpns);
    }
    else if(!strcmp(fun, "i1"))
    {
        score = PQMath::getI1(Ds, k, dim, tmpns);
    }
    else if(!strcmp(fun, "e1"))
    {
        score = PQMath::getE1(Ds, k, dim, tmpns);
    }
    delete [] tmpns;
    delete [] Ds;
    delete [] dat;
    Cleaner::clearClust(clusts);


    return score;
}

float Evaluator::entropy(const char *srcFn, const char *grdFn)
{
    float e;
    map<unsigned int, set<int> *> ourclusts, grdclusts;
    Evaluator::loadClust(srcFn, ourclusts);
    Evaluator::loadClust(grdFn, grdclusts);
    e = Evaluator::entropy(ourclusts, grdclusts);
    Cleaner::clearClust(ourclusts);
    Cleaner::clearClust(grdclusts);
    return e;
}

float Evaluator::purity(const char *srcFn, const char *grdFn)
{
    float p;
    map<unsigned int, set<int> *> ourclusts, grdclusts;
    Evaluator::loadClust(srcFn, ourclusts);
    Evaluator::loadClust(grdFn, grdclusts);
    p = Evaluator::purity(ourclusts, grdclusts);
    Cleaner::clearClust(ourclusts);
    Cleaner::clearClust(grdclusts);
    return p;
}

void Evaluator::getSEP(const char *clustFn, const char *grdFn, const char *datFn, float &score, float &e, float &p, const char *fun)
{
    score = Evaluator::funScore(clustFn, datFn, fun);
    e = Evaluator::entropy(clustFn, grdFn);
    p = Evaluator::purity(clustFn, grdFn);
}



void Evaluator::test()
{
    /**
    const char * ourFn = "tr41_xkb_kpp_10.txt";
    const char * grdFn = "15rclass/tr41.mat.rclass";
    float e = Evaluator::entropy(ourFn, grdFn);
    cout<<"e "<<e<<endl;
    float p = Evaluator::purity(ourFn, grdFn);
    cout<<"p "<<p<<endl;
    /**/
    //const char * clustFn = "15rclass/tr41.mat.rclass";
    /**  char *clustFn = "tr41dst.txt";
      const char *datFn = "norm15dat/tr41.mat";
      const char * grdFn = "15rclass/tr41.mat.rclass";
      const char *fun = "i2";
      //clustFn = "15rclass/tr41.mat.rclass";
      float score, e, p;
      Evaluator::getSEP(clustFn, grdFn, datFn, score, e, p, fun);
      cout<<clustFn<<" score = "<<score<<" e = "<<e<<" p = "<<p<<"\n";
    /**
    const char *srcfn = "/home/chenghaod/dataset/vlad/ox5kbowhe.outs";
    const char *grdfn = "/home/chenghaod/dataset/vlad/ox5k_clust0.txt";
    const char *dstfn = "/home/chenghaod/dataset/vlad/ox5kdst2";
    const char *srcfn1 = "/home/chenghaod/dataset/vlad/vlad_ox5k_dst5.outs";

    Evaluator::jaccard(srcfn, grdfn, dstfn);

    /**/

    const char *srcfn = "/home/chenghaod/data/sift/raw/sift_base.fvecs";
    const char *clusterfn = "/home/chenghaod/dataset/sift1m_pq/sift1m_pq8.out";
    Evaluator::distortion(srcfn, clusterfn);
}
