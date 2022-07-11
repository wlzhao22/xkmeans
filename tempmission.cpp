#include "tempmission.h"

#include "abstractkmeans.h"
#include "minikmeans.h"
#include "dataconvert.h"
#include "evaluator.h"
#include "lvectquant.h"
#include "xbkmeans.h"
#include "rbkmeans.h"
#include "xtkmeans.h"
#include "vstring.h"
#include "ioagent.h"
#include "kmeans.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>

using namespace std;

void TempMission::datEvaluate(const char *datFn, const char *grdFn, float &e, float &p)
{
    char ourFn[1024];
    strcpy(ourFn, datFn);

    map<unsigned int, set<int> *> ourclusts, grdclusts;

    Evaluator::loadClust(ourFn, ourclusts);
    Evaluator::loadClust(grdFn, grdclusts);
    e = Evaluator::entropy(ourclusts, grdclusts);
    p = Evaluator::purity(ourclusts, grdclusts);
}

void TempMission::datEvaluate()
{
    vector<string> filenames;
    const char * srcdir = "dst15dat/";
    string datdir = srcdir;
    string grddir = "15rclass/";
    string grdend = ".mat.rclass";

    IOAgent::readFileName(srcdir, filenames);

    int s = filenames.size();
    vector<string>::iterator it = filenames.begin();
    for(int i = s - 1; i >= 0; i--)
    {
        if(!VString::endWith(filenames[i].c_str(),"norm"))
            filenames.erase(it + i);
    }

    s = filenames.size();
    for(int i = s - 1; i >= 0; i--)
    {
        if(filenames[i].find("large") == filenames[i].npos)
            filenames.erase(it + i);
    }

    s = filenames.size();
    for(int i = s - 1; i >= 0; i--)
    {
        if(filenames[i].find("i1") == filenames[i].npos)
            filenames.erase(it + i);
    }


    s =filenames.size();

    vector<string> datname;

    datname.resize(s);
    char tmpchar[1024];

    float *e = new float[s];
    float *p = new float[s];
    string datFn, grdFn;

    for(int i = 0; i < s; i++)
    {
        VString::parseFName(tmpchar, filenames[i].c_str());
        datname[i] = tmpchar;
    }
    set<string> fname;
    for(int i = 0; i < s; i++)
    {
        fname.insert(datname[i]);
    }

    for(int i = 0; i < s; i++)
    {
        datFn = datdir + filenames[i];
        grdFn = grddir + datname[i] + grdend;
        TempMission::datEvaluate(datFn.c_str(), grdFn.c_str(), e[i], p[i]);
    }

    map<string, int> dat5, dat10, dat15, dat20;
    for(int i = 0; i < s; i++)
    {
        strcpy(tmpchar, datname[i].c_str());
        if(filenames[i].find("_5") != filenames[i].npos)
            dat5.insert(pair<string,int>(tmpchar, i));
        else if(filenames[i].find("10") != filenames[i].npos)
            dat10.insert(pair<string,int>(tmpchar, i));
        else if(filenames[i].find("15") != filenames[i].npos)
            dat15.insert(pair<string,int>(tmpchar, i));
        else if(filenames[i].find("20") != filenames[i].npos)
            dat20.insert(pair<string,int>(tmpchar, i));
        else
        {
            cout<<" erro "<<filenames[i]<<endl;
            exit(0);
        }
    }

    const char * exfile = "result1.txt";
    ofstream oStr(exfile);
    if(!oStr.is_open())
    {
        cout<<"\ncant not open file "<<exfile<<" for write"<<endl;
        exit(0);
    }
    oStr<<setw(8)<<"dataset"<<"\t"<<setw(8)<<"k=5"<<"\t"<<setw(8)<<"k=10"<<"\t"<<setw(8)<<"k=15"<<"\t"<<setw(8)<<"k=20"<<"\t";
    oStr<<setw(8)<<"p=5"<<"\t"<<setw(8)<<"p=10"<<"\t"<<setw(8)<<"p=15"<<"\t"<<setw(8)<<"p=20"<<"\t"<<endl;


    set<string>::iterator sit;
    string id ;
    for(sit = fname.begin(); sit != fname.end(); sit++)
    {
        id = *sit;
        oStr<<setw(8)<<id<<"\t";
        if(dat5.find(id) != dat5.end())
            oStr<<setw(8)<<e[dat5[id]]<<"\t";
        else
            oStr<<setw(8)<<" "<<"\t";
        if(dat10.find(id) != dat10.end())
            oStr<<setw(8)<<e[dat10[id]]<<"\t";
        else
            oStr<<setw(8)<<" "<<"\t";
        if(dat15.find(id) != dat15.end())
            oStr<<setw(8)<<e[dat15[id]]<<"\t";
        else
            oStr<<setw(8)<<" "<<"\t";
        if(dat20.find(id) != dat20.end())
            oStr<<setw(8)<<e[dat20[id]]<<"\t";
        else
            oStr<<setw(8)<<" "<<"\t";
        if(dat5.find(id) != dat5.end())
            oStr<<setw(8)<<p[dat5[id]]<<"\t";
        else
            oStr<<setw(8)<<" "<<"\t";
        if(dat10.find(id) != dat10.end())
            oStr<<setw(8)<<p[dat10[id]]<<"\t";
        else
            oStr<<setw(8)<<" "<<"\t";
        if(dat15.find(id) != dat15.end())
            oStr<<setw(8)<<p[dat15[id]]<<"\t";
        else
            oStr<<setw(8)<<" "<<"\t";
        if(dat20.find(id) != dat20.end())
            oStr<<setw(8)<<p[dat20[id]]<<"\t";
        oStr<<endl;
    }

}
void TempMission::SEPLine(const char *filename, int t, int k)
{
    char normFn[1024];
    char grdFn[1024];
    char dstFn[1024];
    char result[1024];
    char crt[16] = "i2";
    char prr[16] = "large";
    char sed[16] = "non";
    char m[16] = "lvq";
    char rf[16] = "aug";
    char tmprst[1024] = "sepresult";
    char srcMatFn[1024];
    float score[t];
    float e[t];
    float p[t];


    sprintf(srcMatFn, "%s%s", "15dat/", filename);
    sprintf(result, "%s_%s_%s_%s_%s_%s.txt", tmprst, crt, prr, sed, m, rf);
    sprintf(dstFn, "%s%s_%d_%s_%s_%s_%s_%s_%s", "dst", srcMatFn, k, sed, prr, crt, m, rf, ".txt");

    sprintf(normFn, "%s%s", "norm15dat/", filename);
    sprintf(grdFn, "%s%s%s", "15rclass/", filename, ".rclass");

    AbstractKMeans *mykm = NULL;

    for(int i = 0; i < t; i++)
    {
        if(!strcmp(m, "tkm"))
        {
            mykm = new KMeans();
        }
        else if(!strcmp(m, "xbk"))
        {
            mykm = new XBKMeans();
        }
        else if(!strcmp(m, "xtk"))
        {
            mykm = new XTKMeans();
        }
        else if(!strcmp(m, "min"))
        {
            mykm = new MiniKMeans();
        }
        else if(!strcmp(m, "rb"))
        {
            mykm = new RBKMeans();
        }
        else if(!strcmp(m, "lvq"))
        {
            mykm = new LVectQuant();
        }
        cout<<i<<" times\n";
        mykm->buildcluster(normFn, dstFn, sed, prr, crt, k, false);
        Evaluator::getSEP(dstFn, grdFn, normFn, score[i], e[i], p[i], crt);
        cout<<"score "<<"\t"<<"e\t"<<"p\n";
        cout<<score[i]<<"\t"<<e[i]<<"\t"<<p[i]<<"\n";
        delete mykm;
    }

    ofstream oStr(result);
    if(!oStr.is_open())
    {
        cout<<"file "<<result<<" can not open for save sep\n";
        exit(0);
    }
    oStr<<setw(8)<<"score "<<"\t"<<setw(8)<<"e"<<"\t"<<setw(8)<<"p\n";
    for(int i = 0; i < t; i++)
    {
        oStr<<setw(8)<<score[i]<<"\t"<<setw(8)<<e[i]<<"\t"<<setw(8)<<p[i]<<"\n";
    }


    oStr.close();
}
void TempMission::matlabTest()
{
    vector<string> filenames;
    const char * srcdir = "../../matlabsrc/result";
    cout<<"readding srcdir"<<endl;
    IOAgent::readFileName(srcdir,filenames);
    int s =filenames.size();

    char clustFn[1024];
    const char * datFn = "norm15dat/tr41.mat";
    const char * grdFn = "15rclass/tr41.mat.rclass";
    const char * result = "matlabresult.txt";
    const char * crt = "i2";
    float *score = new float[s];
    float *e = new float[s];
    float *p = new float[s];
    for(int i = 0; i < s; i++)
    {
        sprintf(clustFn, "%s%s", "../../matlabsrc/result/", filenames[i].c_str());
        Evaluator::getSEP(clustFn, grdFn, datFn, score[i], e[i], p[i], crt);
        cout<<filenames[i]<<" score = "<<score[i]<<" e = "<<e[i]<<" p = "<<p[i]<<"\n";
    }

    ofstream os(result);
    if(!os.is_open())
    {
        cout<<"file "<<result<<" can not open for write\n";
        exit(0);
    }
    os<<setw(8)<<"score "<<"\t"<<setw(8)<<"e"<<"\t"<<setw(8)<<"p\n";
    for(int i = 0; i < s; i++)
    {
        os<<score[i]<<"\t"<<e[i]<<"\t"<<p[i]<<"\n";
    }
    os.close();
    delete [] score;
    delete [] e;
    delete [] p;
    filenames.clear();
}

void TempMission::yzscore()
{
    vector<string> filenames;
    vector<string> trimedfilenames;
    const char * srcdir = "/home/chenghaod/Downloads/datasets/15dat";
    cout<<"readding srcdir"<<endl;
    IOAgent::readFileName(srcdir,filenames);
    int s =filenames.size();
    vector<string>::iterator it = filenames.begin();
    for(int i = s - 1; i >= 0; i--)
    {
        if(filenames[i].find("clustering") == filenames[i].npos)
            filenames.erase(it + i);
    }


    s = filenames.size();
    trimedfilenames.resize(s);

    char tmp[1024];

    for(int i = 0; i < s; i++)
    {
        int ed = VString::firstindexof(filenames[i].c_str(), '.');
        memset(tmp, 0, 1024);
        memcpy(tmp, filenames[i].c_str(), ed*sizeof(char));
        trimedfilenames[i] = tmp;
        cout<<filenames[i]<<" s"<<"\n";
        cout<<tmp<<" t"<<"\n";
    }


    char clustFn[1024];
    char datFn[1024];
    char grdFn[1024];
    const char *result = "yzresult.txt";
    float *score = new float[s];
    float *e = new float[s];
    float *p = new float[s];
    for(int i = 0; i < s; i++)
    {
        sprintf(clustFn, "%s/%s", srcdir, filenames[i].c_str());
        sprintf(datFn, "%s/%s%s", "norm15dat", trimedfilenames[i].c_str(), ".mat");
        sprintf(grdFn, "%s/%s%s", "15rclass", trimedfilenames[i].c_str(), ".mat.rclass");

        Evaluator::getSEP(clustFn, grdFn, datFn, score[i], e[i], p[i], "i2");
        cout<<filenames[i]<<" score = "<<score[i]<<" e = "<<e[i]<<" p = "<<p[i]<<"\n";
    }

    ofstream os(result);
    if(!os.is_open())
    {
        cout<<"file "<<result<<" can not open for write\n";
        exit(0);
    }

    os<<setw(16)<<"filename"<<"\t"<<setw(8)<<"score "<<"\t"<<setw(8)<<"e"<<"\t"<<setw(8)<<"p\n";

    for(int i = 0; i < s; i++)
    {
        os<<setw(16)<<filenames[i]<<"\t"<<score[i]<<"\t"<<e[i]<<"\t"<<p[i]<<"\n";
    }
    os.close();
    delete [] score;
    delete [] e;
    delete [] p;
    filenames.clear();

}


