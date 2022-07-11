#ifndef CLEANER_H
#define CLEANER_H

#include "clsinfo.h"
#include "nnitem.h"


#include <vector>
#include <thread>
#include <list>
#include <map>
#include <set>

/**
In charge of memory recycling, acts as
a rubbish collector

@author:    Wanlei Zhao
@email:     stonescx@gmail.com
@date:      30-Jul-2012
@Institute: INRIA
**/

using namespace std;

class Cleaner
{


public:

    static void clearNNList(list<NNItem*> &nnItems);

    static void clearMatrix(vector<vector<unsigned int> > &matrix);

    static void clear_k2iMap(map<string, unsigned int>      &refTab);
    static void clear_k2iMap(map<string, unsigned char>     &refTab);
    static void clear_i2kMap(map<unsigned int, const char*> &i2kMap);

    static void freeItemMaps(map<unsigned int, map<string, const char*> *> &itmMaps);

    static void freeParaMap(map<string, const char*> &paras);
    static void freeStrVect(vector<const char*> &strVect);
    static void clearThread(vector<std::thread*> &thrds);

    static void clearClust(map<unsigned int, set<int> *> &clust);

    static bool clearVector(vector<NNItem*> &vects);
    static void clearInfoMap(map<unsigned int, CLSInfo*> &infoMap);

    static void clearMap(map<int, string> &is);
    static void clearMap(map<string, string> &ss);
    static void clearMap(map<string, int> &si);

    static void clearkNNGraph(vector<vector<MiniNN> > &knnGraph);
    static void clearNbs(vector<NbHood> &nbGraph);
    static void clearDs(vector<vector<double>> &Ds);
};

#endif
