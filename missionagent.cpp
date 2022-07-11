#include "missionagent.h"

#include "abstractkmeans.h"
#include "dataconvert.h"
#include "minikmeans.h"
#include "lvectquant.h"
#include "evaluator.h"
#include "xbkmeans.h"
#include "xtkmeans.h"
#include "rbkmeans.h"
#include "hartigan.h"
#include "xtksums.h"
#include "xbksums.h"
#include "ioagent.h"
#include "cleaner.h"
#include "vstring.h"
#include "kmeans.h"
#include "ksums.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>


using namespace std;

bool MissionAgent::buildClust(map<string, const char*>   &arguments)
{
    const char *paras[] = {"-i", "-d", "-k", "-m"};
    char seed[16], crt[16], prr[16];
    unsigned int i = 0;
    bool _refine_  = false;
    for(i = 0; i < 4; i++)
    {
        if(arguments.find(paras[i]) == arguments.end())
        {
            cout<<"Required parameter '"<<paras[i]<<"' is missing!\n";
            return false;
        }
    }
    AbstractKMeans *mykmn = NULL;

    if(!strcmp(arguments["-m"], "tkm"))
    {
        mykmn = new KMeans();
    }
    else if(!strcmp(arguments["-m"], "xbk"))
    {
        mykmn = new XBKMeans();
    }

    if(arguments.find("-s") != arguments.end())
    {
        strcpy(seed, arguments["-s"]);
    }
    else
    {
        strcpy(seed, "kpp");
    }

    if(arguments.find("-r") != arguments.end())
    {
        _refine_ = true;
    }


    if(!strcmp(seed, "rnd") || !strcmp(seed, "kpp")||!strcmp(seed, "non"))
    {
        cout<<" Unknown seeding option '"<<crt<<"'!\n";
        cout<<" Valid options are 'rnd', 'kpp' or 'non'!\n";
        exit(0);
    }

    if(arguments.find("-p") != arguments.end())
    {
        strcpy(prr, arguments["-p"]);
    }
    else
    {
        strcpy(prr, "large");
    }

    if(!strcmp(prr, "large") || !strcmp(prr, "best"))
    {
        cout<<" Unknown seeding option '"<<crt<<"'!\n";
        cout<<" Valid options are 'large' or 'best'!\n";
        exit(0);
    }

    if(arguments.find("-c") != arguments.end())
    {
        strcpy(crt, arguments["-c"]);
    }
    else
    {
        strcpy(crt, "i2");
    }

    int clustNum = atoi(arguments["-k"]);

    if(clustNum < 0 || clustNum > 2147483648)
    {
        cout<<" Invalid cluster number!\n";
        cout<<" Suggested range of cluster number: 1 < k < matrix size\n";
        exit(0);
    }

    mykmn->buildcluster(arguments["-i"], arguments["-d"], seed, prr, crt, clustNum, _refine_);
    delete mykmn;

    return true;
}

void MissionAgent::yingzhaoSolution(const char *fileName,const char *mth, float *e, float *p, const char *crt, const char *prr, const char *sed, const char *rf)
{

    cout<<"file "<<fileName<<endl<<endl;
    string name ;
    char dataFn[1024];
    const char *record = "/home/wlzhao/datasets/umn/record_hart_Im.txt";
    string dir = "/home/wlzhao/datasets/umn/";
    string dstDir = "/home/wlzhao/datasets/umn/rslt/";
    name = "15dat/";
    name = name + fileName;
    strcpy(dataFn, name.c_str());

    //name = "norm/";
    //name = dir + name + fileName;

    //if(!VString::existFile(name.c_str()))
    //DataCovert::yingzhaoForm(dataFn);
    name = "norm";
    name = dir + name + dataFn;
    string srcMatFn = name;


    name = "15rclass/";
    name = dir + name + fileName + ".rclass";
    char grdFn[1024];
    strcpy(grdFn, name.c_str());

    AbstractKMeans *mykm = NULL;

    const int a[4] = {5,10,15,20};
    int k = 0, i = 0;
    for(i = 0; i < 4; i++)
    {
        k = a[i];

        cout<<"\n\nk = "<<k<<endl;
        char dstFn[1024];
        char t[16];
        sprintf(t, "%d", k);
        if(!strcmp(mth, "tkm"))
        {
            mykm = new KMeans();
        }
        else if(!strcmp(mth, "xbk"))
        {
            mykm = new XBKMeans();
        }
        else if(!strcmp(mth, "xtk"))
        {
            mykm = new XTKMeans();
        }
        else if(!strcmp(mth, "min"))
        {
            mykm = new MiniKMeans();
        }
        else if(!strcmp(mth, "rbk"))
        {
            mykm = new RBKMeans();
        }
        else if(!strcmp(mth, "lvq"))
        {
            mykm = new LVectQuant();
        }
        else if(!strcmp(mth, "xtksum"))
        {
            mykm = new XTKSums();
        }
        else if(!strcmp(mth, "xbksum"))
        {
            mykm = new XBKSums();
        }else if(!strcmp(mth, "ksum"))
        {
            mykm = new KSums();
        }else if(!strcmp(mth, "hart"))
        {
            mykm = new Hartigan();
        }

        sprintf(dstFn, "%s_dst_%s_%s_%d.txt", dstDir.c_str(), fileName, mth, a[i]);

        mykm->buildcluster(srcMatFn.c_str(), dstFn, sed, prr, crt, k, false);

        delete mykm;

        map<unsigned int, set<int> *> ourclusts, grdclusts;
        Evaluator::loadClust(dstFn, ourclusts);
        Evaluator::loadClust(grdFn, grdclusts);
        e[i] = Evaluator::entropy(dstFn, grdFn);
        p[i] = Evaluator::purity(dstFn, grdFn);
        ofstream *outStrm = new ofstream(record, ios::app);
        ///cout<<"e \t"<<" p"<<endl;
        (*outStrm)<<mth<<"\t"<<fileName<<"\t\t"<<a[i]<<"\t"<<e[i]<<"\t"<<p[i]<<endl;
        ///cout<<mth<<"\t"<<fileName<<"\t\t"<<a[i]<<"\t"<<e[i]<<"\t"<<p[i]<<endl;
        outStrm->close();
        Cleaner::clearClust(ourclusts);
        Cleaner::clearClust(grdclusts);
    }

}

void MissionAgent::experiment()
{
    vector<string> filenames;
    const char * srcdir = "15dat";
    cout<<"reading srcdir"<<endl;
    IOAgent::readFileName(srcdir, filenames);
    int s =filenames.size();

    vector<string>::iterator it = filenames.begin();
    for(int i = s - 1; i >= 0; i--)
    {
        if(!VString::endWith(filenames[i].c_str(),".mat"))
            filenames.erase(it + i);
    }

    s = filenames.size();

    float *e = new float[s*4];
    float *p = new float[s*4];
    int id = 0;
    char crt[16] = "i2";
    char prr[16] = "large";
    char sed[16] = "non";
    char m[16] = "lvq";
    char rf[16] = "oug";
    char tmpresult[15024];
    char tmprst[1024] = "tmpresult";

    sprintf(tmpresult, "%s_%s_%s_%s_%s_%s.txt", tmprst, crt, prr, sed, m, rf);

    for(int i = 0; i < s; i++)
    {
        cout<<" FILE "<<i<<endl;
        MissionAgent::yingzhaoSolution(filenames[i].c_str(), m, e + id, p + id, crt, prr, sed, rf);
        IOAgent::addRslt(filenames[i].c_str(), e + id, p + id, tmpresult);
        id += 4;
    }

    char result[1024];
    sprintf(result, "%s_%s_%s_%s_%s_%s%s", "result", crt, prr, sed, m, rf ,".txt");
    IOAgent::saveEx(filenames, e, p, result);

    delete [] e;
    delete [] p;

}







