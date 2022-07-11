#include "cleaner.h"

#include <cassert>

void Cleaner::clearMatrix(vector<vector<unsigned int> > &matrix)
{
    vector<vector<unsigned int> >::iterator vit;
    for(vit = matrix.begin(); vit != matrix.end(); vit++)
    {
        vector<unsigned int> &crntVect = *vit;
        crntVect.clear();
    }
    matrix.clear();

}

void Cleaner::freeItemMaps(map<unsigned int, map<string, const char*> *> &itmMaps)
{
    map<unsigned int, map<string,const char*>* >::iterator it;
    map<string, const char*>::iterator iter;
    map<string, const char*> *crntMap;
    string keystr;
    for(it = itmMaps.begin(); it != itmMaps.end(); it++)
    {
        crntMap = it->second;
        for(iter = crntMap->begin(); iter != crntMap->end(); iter++)
        {
            keystr = iter->first;
            /**i am not clear whether we should clear following
            strings by ourselves**/
            /**
            const char *val = iter->second;
            delete [] val;
            keystr.clear();
            **/
        }
        crntMap->clear();
        delete crntMap;
    }
    itmMaps.clear();
}

void Cleaner::clear_i2kMap(map<unsigned int, const char*> &i2kMap)
{
    map<unsigned int,const char*>::iterator mit;

    const char *val;
    for(mit = i2kMap.begin(); mit != i2kMap.end(); mit++)
    {
        val = mit->second;
        delete [] val;
    }
    i2kMap.erase(i2kMap.begin(),i2kMap.end());

}

void Cleaner::freeParaMap(map<string, const char*> &paras)
{
    map<string,const char*>::iterator mit;
    string crntstr;
    const char *val;
    for(mit = paras.begin(); mit != paras.end(); mit++)
    {
        val = mit->second;
        crntstr = mit->first;
        crntstr.erase(crntstr.begin(),crntstr.end());
        delete [] val;
    }
    paras.erase(paras.begin(),paras.end());

}


void Cleaner::freeStrVect(vector<const char*> &str_vect)
{
    const char *crntChrs;
    vector<const char*>::iterator it;
    for(it = str_vect.begin(); it != str_vect.end(); it++)
    {
        crntChrs = *it;
        delete [] crntChrs;
    }
    str_vect.clear();
    return ;
}

void Cleaner::clear_k2iMap(map<string, unsigned int> &refTab)
{
    map<string, unsigned int>::iterator mit;
    string key;

    for(mit = refTab.begin(); mit != refTab.end(); mit++)
    {
        key = mit->first;
        key.clear();
    }

    refTab.clear();
}

void Cleaner::clear_k2iMap(map<string, unsigned char> &refTab)
{
    map<string, unsigned char>::iterator mit;
    string key;

    for(mit = refTab.begin(); mit != refTab.end(); mit++)
    {
        key = mit->first;
        key.clear();
    }

    refTab.clear();
}

void Cleaner::clearNNList(list<NNItem*> &nnItems)
{
    NNItem *crnt_itm = NULL;
    list<NNItem*>::iterator itm_iter;

    for(itm_iter = nnItems.begin(); itm_iter != nnItems.end(); itm_iter++)
    {
        crnt_itm = *itm_iter;
        delete crnt_itm;
    }
    nnItems.erase(nnItems.begin(), nnItems.end());
}

bool Cleaner::clearVector(vector<NNItem*> &vects)
{
    vector<NNItem*>::iterator vit;
    NNItem *crntItm;
    for(vit = vects.begin(); vit != vects.end(); vit++)
    {
        crntItm = *vit;
        delete crntItm;
    }
    vects.clear();
    return true;
}

void Cleaner::clearInfoMap(map<unsigned int, CLSInfo*> &infoMap)
{
    map<unsigned int, CLSInfo*>::iterator mit;
    CLSInfo *crnt_info;
    for(mit = infoMap.begin(); mit != infoMap.end(); mit++)
    {
        crnt_info = mit->second;
        delete crnt_info;
    }
    infoMap.clear();
    return ;
}


void Cleaner::clearThread(vector<std::thread*> &thrds)
{
    vector<std::thread*>::iterator it;
    std::thread* thrd;
    for(it = thrds.begin(); it != thrds.end(); it++)
    {
        thrd = *it;
        delete thrd;
    }
    thrds.clear();
}

void Cleaner::clearClust(map<unsigned int, set<int> *> &clust)
{
    for(auto it = clust.begin(); it != clust.end(); it++)
    {
        it->second->clear();
    }
    clust.clear();

}

void Cleaner::clearMap(map<int, string> &is)
{
    map<int, string>::iterator it;

    for(it = is.begin(); it != is.end(); it++)
    {
        (it->second).clear();
    }
    is.clear();
}

void Cleaner::clearMap(map<string, string> &is)
{
    map<string, string>::iterator it;
    string key;
    for(it = is.begin(); it != is.end(); it++)
    {
        key = it->first;
        key.clear();
        (it->second).clear();
    }
    is.clear();
}

void Cleaner::clearMap(map<string, int > &is)
{
    map<string, int>::iterator it;
    string key;
    for(it = is.begin(); it != is.end(); it++)
    {
        key = it->first;
        key.clear();
    }
    is.clear();
}

void Cleaner::clearkNNGraph(vector<vector<MiniNN> > &knnGraph)
{
    unsigned int i = 0;

    for(i = 0; i < knnGraph.size(); i++)
    {
        vector<MiniNN> & crntLst = knnGraph[i];
        crntLst.clear();
    }
    knnGraph.clear();
}

void Cleaner::clearNbs(vector<NbHood> &nbGraph)
{
    unsigned int i = 0;

    for(i = 0; i < nbGraph.size(); i++)
    {
        nbGraph[i].newnb.clear();
        nbGraph[i].oldnb.clear();
        nbGraph[i].rnewnb.clear();
        nbGraph[i].roldnb.clear();
    }
}

void Cleaner::clearDs(vector<vector<double>> &Ds)
{
    unsigned int i = 0;

    for(i = 0; i < Ds.size(); i++)
    {
        vector<double>& D = Ds[i];
        D.clear();
    }
    Ds.clear();
}


