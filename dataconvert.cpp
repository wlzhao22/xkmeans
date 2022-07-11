#include "dataconvert.h"
#include "ioagent.h"
#include "vstring.h"
#include "cleaner.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>
#include <cmath>
#include <map>

using namespace std;

void DataConvert::getIdf(float *mat, unsigned int row, unsigned col, float *idf)
{
    double N = static_cast<double>(row);
    int *df = new int[col];
    int id = 0;
    for(unsigned int i = 0; i < col; i++)
    {
        df[i] = 0;
        id = i;
        idf[i] = 0;
        for(unsigned int j = 0; j < row; j++)
        {
            if(mat[id] != 0)
                df[i]++;
            id += col;
        }
    }
    for(unsigned int i = 0; i < col; i++)
    {

        if(df[i] != 0)
        {
            idf[i] = log(N/df[i]);
        }
    }
    delete []df;
}

void DataConvert::tfidfMat(const float *mat, const float *idf, unsigned int row, unsigned col, float *nmat)
{
    int rowId = 0;

    for(unsigned int i = 0; i < row; i++)
    {
        rowId = i*col;
        for(unsigned int j = 0; j < col; j++)
        {
            nmat[rowId + j] = mat[rowId + j]*idf[j];
        }
    }

}

void DataConvert::normMat(float *mat, unsigned int row, unsigned col, float *nmat)
{
    float ln = 0;
    int rowId;

    for(unsigned int i = 0; i < row; i++)
    {
        rowId = i*col;
        ln = 0;
        for(unsigned int j = 0; j < col; j++)
        {
            ln += mat[rowId + j]*mat[rowId + j];
        }
        ln = sqrt(ln);
        for(unsigned int j = 0; j < col; j++)
        {
            nmat[rowId + j] = mat[rowId + j]/ln;
        }
    }
}

void DataConvert::yingzhaoForm(const char *dataFn, const char *srcDir, const char *dstDir)
{
    string srcdir = srcDir;
    string dstdir = dstDir;
    string srcFn  = dataFn;
    string normFn = dstdir + srcFn;
    srcFn = srcdir + srcFn;

    float *mat    = NULL;
    float *idf    = NULL;
    float *idfmat = NULL;
    float *nmat   = NULL;
    unsigned int row = 0, col = 0;

    if(VString::endWith(srcFn.c_str(), ".txt"))
    {
        mat = IOAgent::loadMatrix(srcFn.c_str(), row, col);
    }
    else if(VString::endWith(srcFn.c_str(), ".fvecs"))
    {
        mat = IOAgent::load_fvecs(srcFn.c_str(), col, row);
    }
    else if(VString::endWith(srcFn.c_str(), ".mat"))
    {
        mat = IOAgent::loadDat(srcFn.c_str(), row, col);
    }
    else
    {
        cout<<"Unrecognizable input file format!!!\n";
        mat = NULL;
        exit(0);
    }
    idf = new float[col];
    DataConvert::getIdf(mat, row, col, idf);
    /**/
    //IOAgent::saveDat(idfFn.c_str(), 1, col, idf);
    idfmat = new float[row*col];
    DataConvert::tfidfMat(mat, idf, row, col, idfmat);
    ///IOAgent::saveDat(tfidfFn.c_str(), row, col, idfmat);
    nmat = new float[row*col];
    DataConvert::normMat(idfmat, row, col, nmat);
    IOAgent::saveDat(normFn.c_str(), row, col, nmat);
    /**/
    delete [] mat;
    delete [] nmat;
    delete [] idf;
    delete [] idfmat;
    idf = idfmat = nmat = mat = NULL;
}

void DataConvert::iggg2ii(const char *igfn, const char * ggfn, const char *iifn)
{
    ifstream ig(igfn);
    if(!ig.is_open())
    {
        cout<<"can not open file "<<igfn<<" for read\n";
        exit(0);
    }
    ifstream gg(ggfn);
    if(!gg.is_open())
    {
        cout<<"can not open file "<<ggfn<<" for read\n";
        exit(0);
    }
    ofstream ii(iifn);
    if(!ii.is_open())
    {
        cout<<"can not open write "<<iifn<<" for read\n";
        exit(0);
    }


    map<int, string> i2g;
    map<string, int> g2i;
    map<string, string> g2g;
    ///  multimap<string, string> g2g;
    vector< pair<int, int> > i2i;
    char tmps[1024];
    char tmps1[1024];

    int tmpi;


    while(ig>>tmpi)
    {
        ig>>tmps;

        i2g.insert(pair<int, string>(tmpi, tmps));
    }

    ig.close();

    while(gg>>tmps)
    {
        gg>>tmps1;
        g2g.insert(pair<string, string>(tmps, tmps1));
    }


    gg.close();

    for(map<int, string>::iterator it = i2g.begin(); it != i2g.end(); it++)
    {
        g2i.insert(pair<string, int>(it->second, it->first));
    }



    for(map<int, string>::iterator it = i2g.begin(); it != i2g.end(); it++)
    {
        i2i.push_back(pair<int, int >(it->first, g2i[g2g[it->second]]));
    }

    for(unsigned int i = 0; i < i2i.size(); i++)
    {
        ii<<i2i[i].first<<"\t"<<i2i[i].second<<"\n";
    }

    ii.close();

    Cleaner::clearMap(i2g);
    Cleaner::clearMap(g2i);
    Cleaner::clearMap(g2g);

}

void DataConvert::label2set(vector<int> labels, map<int,set<int> *> &clusts)
{
    int i = 0, n = 0;
    vector<int> clust;
    set<int> temp;
    set<int>::iterator it;
    map<int, int> mp;
    n = labels.size();
    for(int i = 0; i < n; i++)
    {
        temp.insert(labels[i]);
    }

    i = 0;
    for(it = temp.begin(); it != temp.end(); it++)
    {
        set<int > * tmpclust = new set<int>;
        clusts.insert(pair<int, set<int> *>(i, tmpclust));
        i++;

    }

    for(int i = 0; i < n; i++)
    {
        clusts[labels[i]]->insert(i);
    }

}

void DataConvert::label2set(vector<pair<int, int> >labels, map<int, set<int > *> &sets)
{
    int i = 0, n = 0;
    n = labels.size();
    vector<set<int> > initset;
    set<int > tmpset;
    for(i = 0; i < n; i++)
    {
        tmpset.clear();
        tmpset.insert(labels[i].first);
        tmpset.insert(labels[i].second);
        initset.push_back(tmpset);
    }

}

void DataConvert::label2label(const char  *fn, const char* igfn, const char *dstfn)
{
    ifstream is(fn);
    ifstream ig(igfn);
    ofstream os(dstfn);

    int i, n, s;

    vector<set<int> *> cluster;

    map<int, string> igmap;
    map<string, int> gimap;
    set<int>::iterator it;

    char tmps[1024];

    int tmpi;

    set<int> *pset;

    while(ig>>tmpi)
    {
        ig>>tmps;
        igmap.insert(pair<int, string>(tmpi, tmps));
    }

    ig.close();

    for(map<int, string>::iterator it = igmap.begin(); it != igmap.end(); it++)
    {
        gimap.insert(pair<string, int>(it->second, it->first));
    }

    while(is>>tmpi)
    {
        pset = new set<int>;
        for(i = 0; i < tmpi; i++)
        {
            is>>tmps;
            pset->insert(gimap[tmps]);
        }
        cluster.push_back(pset);
    }

    is.close();

    n = cluster.size();

    for(i = 0; i < n; i++)
    {
        s = cluster[i]->size();
        os<<s<<"\t";
        for(it = cluster[i]->begin(); it != cluster[i]->end(); it++)
        {
            os<<*it<<"\t";
        }
        os<<"\n";
    }

    for(i = n-1; i >= 0; i--)
    {
        cluster[i]->clear();
    }

    os.close();

}

void DataConvert::similarityBoost(SparseMatrix &sdata, unsigned int row, unsigned int ndim, float boostRate)
{
    unsigned int i, j, id, loc;

    unsigned int *data = new unsigned[row*ndim];

    memset(data, 0, sizeof(unsigned)*row*ndim);

    for(i = 0; i < row; i++)
    {
        loc = i*ndim;
        for(j = sdata.col[i]; j < sdata.col[i+1]; j++)
        {
            data[loc + sdata.index[j]] = 1;
        }
    }

    unsigned int *boostMap = new unsigned[ndim];
    for(i = 0; i < ndim; i++)
    {
        boostMap[i] = i;
    }
    random_shuffle(boostMap, boostMap + ndim);

    /**/
    id = 0;
    if(boostRate <=1)
        for(i = 0; i < row; i++)
        {
            loc = i*ndim;
            for(j = sdata.col[i]; j < sdata.col[i+1]; j++)
            {
                if(rand() < boostRate*RAND_MAX)
                    data[loc + boostMap[sdata.index[j]]] = 1;
            }
        }
    else
    {
        while(id < boostRate)
        {
            random_shuffle(boostMap, boostMap + ndim);
            for(i = 0; i < row; i++)
            {
                loc = i*ndim;
                for(j = sdata.col[i]; j < sdata.col[i+1]; j++)
                {
                    if(rand() < boostRate*RAND_MAX)
                        data[loc + boostMap[sdata.index[j]]] = 1;
                }
            }
            id++;
        }
    }
    /**/
    /**
        for(i = 0; i < ndim; i++)
        {
            if(rand() < RAND_MAX*boostRate)
            for(j = 0; j < row; j++)
            {
                if(data[j*ndim + i] == 1)
                data[j*ndim + boostMap[i]] == 1;
            }

        }
    **/
    delete [] sdata.col;
    delete [] sdata.data;
    delete [] sdata.index;

    id = 0;
    for(i = 0; i < row; i++)
    {
        loc = i*ndim;
        for(j = 0; j < ndim; j++)
            if(data[loc + j] == 1)
                id++;
    }

    sdata.col = new unsigned[row+1];
    sdata.index = new unsigned[id+1];
    sdata.data = new float[id];
    id = 0;
    for(i = 0; i < row; i++)
    {
        loc = i*ndim;
        sdata.col[i] = id;
        for(j = 0; j < ndim; j++)
            if(data[loc + j] == 1)
            {
                sdata.data[id] = 1;
                sdata.index[id] = j;
                id++;
            }
    }
    sdata.col[row] = id;
    delete [] data;
}

void DataConvert::rawcsv2mat(const char *srcFn, const char *dstFn)
{
    unsigned int count = 0, i = 0;
    unsigned int row = 0,  col = 0;
    string txt;
    vector<float> vals;
    ifstream *inStrm = new ifstream(srcFn, ios::in);
    cout<<srcFn<<endl;
    while(inStrm->good())
    {
        getline(*inStrm, txt, '\n');
        if(col == 0)
        {
            col = VString::countsof(txt, ',');
        }
        if(txt.length() > 1)
            count++;
    }

    if(count == 0)
    {
        return ;
    }

    row = count;

    inStrm->clear();
    inStrm->seekg(0, ios::beg);
    cout<<"Row: "<<row<<" col: "<<col<<endl;
    ofstream *outStrm = new ofstream(dstFn, ios::out);
    (*outStrm)<<100000<<" "<<col+1<<endl;
    count = 0;
    while(inStrm->good())
    {
        getline(*inStrm, txt, '\n');
        if(txt.size() > 1)
        {
            VString::str2float(txt.c_str(), ',', vals);
            for(i = 0; i < vals.size(); i++)
            {
                (*outStrm)<<vals[i]<<" ";
            }
            (*outStrm)<<endl;
            vals.clear();
            txt.clear();
        }
        count++;
        if(count == 100000)
        break;
    }
    outStrm->close();
    inStrm->close();
    return ;
}

void DataConvert::xcsv2mat(const char *srcFn, const char *dstFn)
{
    unsigned int count = 0, i = 0;
    unsigned int row = 0,  col = 0;
    string txt;
    vector<float> vals;
    float *maxVals = NULL;
    float *minVals = NULL;
    ifstream *inStrm = new ifstream(srcFn, ios::in);
    cout<<srcFn<<endl;
    while(inStrm->good())
    {
        getline(*inStrm, txt, '\n');
        if(col == 0)
        {
            col = VString::countsof(txt, ',');
        }
        if(txt.length() > 1)
            count++;
    }

    if(count == 0)
    {
        return ;
    }

    row = count;

    inStrm->clear();
    inStrm->seekg(0, ios::beg);

    cout<<"Row: "<<row<<" col: "<<col<<endl;
    maxVals = new float[col];
    minVals = new float[col];
    for(i = 0; i < col; i++)
    {
       maxVals[i] = 0;
       minVals[i] = RAND_MAX;
    }

    count = 0;
    while(inStrm->good())
    {
        getline(*inStrm, txt, '\n');
        if(txt.size() > 1)
        {
            VString::str2float(txt.c_str(), ',', vals);
            for(i = 1; i < vals.size(); i++)
            {
                maxVals[i-1] = (maxVals[i-1] < vals[i])?vals[i]:maxVals[i-1];
                minVals[i-1] = (minVals[i-1] > vals[i])?vals[i]:minVals[i-1];
            }
            vals.clear();
            txt.clear();
        }
        count++;
    }
    for(i = 0; i < col; i++)
    {
       maxVals[i] = maxVals[i] - minVals[i];
       if(maxVals[i] == 0)
       {
           maxVals[i] = 1;
       }
    }

    inStrm->clear();
    inStrm->seekg(0, ios::beg);
    for(i = 0; i < col; i++)
    {
         cout<<"[ "<<minVals[i]<<", "<<maxVals[i]<<" ]"<<endl;
    }

    ofstream *outStrm = new ofstream(dstFn, ios::out);
    (*outStrm)<<count<<" "<<col<<endl;
    count = 0;
    while(inStrm->good())
    {
        getline(*inStrm, txt, '\n');
        if(txt.size() > 1)
        {
            VString::str2float(txt.c_str(), ',', vals);
            for(i = 1; i < vals.size(); i++)
            {
                (*outStrm)<<(vals[i] - minVals[i-1])/maxVals[i-1]<<" ";
            }
            (*outStrm)<<endl;
            vals.clear();
            txt.clear();
        }
        count++;
        if(count >= 100000)
        break;
    }
    outStrm->close();
    inStrm->close();
    delete [] maxVals;
    delete [] minVals;
    maxVals = minVals = NULL;

    return ;
}

void DataConvert::test()
{
    const char *srcFn1 = "/home/wlzhao/datasets/clust/USCensus1990.data.txt";
    const char *dstFn1 = "/home/wlzhao/datasets/clust/USCensus1990_2mxxx.txt";
    const char *srcFn2 = "/home/wlzhao/datasets/clust/SUSY.csv";
    const char *dstFn2 = "/home/wlzhao/datasets/clust/susy_mat2m.txt";
    const char *srcs[15]= {"classic.mat", "fbis.mat", "hitech.mat", "k1a.mat", "k1b.mat", "la12.mat",
    "new3.mat", "ohscal.mat", "re0.mat", "re1.mat", "reviews.mat", "sports.mat", "tr31.mat", "tr41.mat", "wap.mat"};
    const char *srcDir = "/home/wlzhao/datasets/umn/umn/";
    const char *dstDir = "/home/wlzhao/datasets/umn/norm/";

    ///DataConvert::xcsv2mat(srcFn1, dstFn1);
    for(unsigned int i = 0; i < 15; i++)
    DataConvert::yingzhaoForm(srcs[i], srcDir, dstDir);
}


