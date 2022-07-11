#include <iostream>

#include "missionagent.h"
#include "dataconvert.h"
#include "minikmeans.h"
#include "lvectquant.h"
#include "evaluator.h"
#include "randseed.h"
#include "rbkmeans.h"
#include "xbkmeans.h"
#include "xtkmeans.h"
#include "vstring.h"
#include "ioagent.h"
#include "kmeans.h"
#include "speedtest.h"
#include "xtksums.h"
#include "xbksums.h"
#include "gksums.h"
#include "ksums.h"
#include "seqkmeans.h"

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <map>

using namespace std;

void help()
{
    const char *version = "1.38";

    cout<<" xkmeans -i mat -k num -m xbk|tkm|xtk|rbk [-s kpp|rnd|non] [-c i1|i2|e1|e2] [-p best|large] [-r] -d fn\n";
    cout<<"    -i  mat\tinput matrix\n";
    cout<<"    -k  num\tnumber of clusters\n";
    cout<<"    -m  xbk|.\tkmean method option\n";
    cout<<"    -s  kpp|.\tseed option (optional)(kpp by default)\n";
    cout<<"    -r  \twith refinement (without refinement by default)\n";
    cout<<"    -c  i1|.\toptimization objective function (i2 by default)\n";
    cout<<"    -p  best|.\tbisecting stradegy (large by default)\n";
    cout<<"    -d  fn\tdestine visual words file\n\n";

    cout<<"  Author:\tWan-Lei Zhao, Cheng-Hao Deng and Hui Ye, email to wlzhao@xmu.edu.cn\n";
    cout<<"  Copyrights:\tAll rights are reserved by the atuhors\n";
    cout<<"  Version:\t"<<version<<endl;
}

void test()
{
    ///RandSeed::test();
    ///KMeans::test();
    ///IOAgent::test();
    ///RRClust::test();
    ///DataConvert::test();
    ///XBKMeans::test();
    ///XTKMeans::test();
    ///PQTrainer::test();
    ///RBKMeans::test();
    ///LVectQuant::test();
    ///MiniKMeans::test();
    ///MiniKMeans::pctest();
    ///Evaluator::test();
    ///MissionAgent::experiment();
    /**
    const char *files[15] = {"classic.mat", "fbis.mat", "hitech.mat", "k1a.mat", "k1b.mat", "la12.mat",
    "new3.mat", "ohscal.mat", "re0.mat", "re1.mat", "reviews.mat", "sports.mat", "tr31.mat", "tr41.mat", "wap.mat"};
    float e[4] = {0}, p[4] = {0};
    for(unsigned i = 0; i < 15; i++)
    {
       ///MissionAgent::yingzhaoSolution(files[i], "xtksum", e, p, "i2", "large", "rnd", "-r");
       MissionAgent::yingzhaoSolution(files[i], "hart", e, p, "i2", "large", "rnd", "-r");
    }
    /**/
    ///SpeedTest::test();
    ///KMeans::test();
    ///KSums::test();
    ///KSums::wltest();
    ///GKSums::wltest();
    ///XBKSums::test();
    ///XTKSums::wltest();
}


int main(int argc, const char *argv[])
{
    ///srand(time(NULL));
    test();
    exit(0);

    if(argc < 9)
    {
        help();
        exit(0);
    }

    int i = 0;
    map<string, const char*> arguments;

    for(i = 1; i < argc; i += 2)
    {
        arguments.insert(pair<string, const char*>(argv[i], argv[i+1]));
    }

    MissionAgent::buildClust(arguments);

    arguments.clear();
    return 0;
}
