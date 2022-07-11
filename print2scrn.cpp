#include "print2scrn.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>

using namespace std;

void Print2Scrn::printCLInfo(const CLSInfo *info, const unsigned int n)
{
    cout<<"-------------------------------\n";
    for(unsigned int i = 0; i < n; i++)
    {
        cout<<i<<"\t"<<info[i].n<<endl;
    }
    cout<<"================================\n";
}

void Print2Scrn::printCLInfo(const CLSInfo *info, const unsigned int n,
                                 const char *dstfn)
{
    ofstream *outStrm = new ofstream(dstfn, ios::out);
    for(unsigned int i = 0; i < n; i++)
    {
        (*outStrm)<<i<<"\t"<<info[i].n<<endl;
    }
    outStrm->close();
    return ;
}

void Print2Scrn::printvect(vector<NNItem*> &vects)
{
    vector<NNItem*>::iterator vit;
    NNItem* crntItm = NULL;
    cout<<"Solution------------------\n";
    cout<<"  "<<vects.size()<<" ways of clustering\n";
    cout<<"--------------------------\n";
    cout<<" Clabel\t  Size\tISim\t|\n";
    cout<<"--------------------------\n";
    for(vit = vects.begin(); vit != vects.end(); vit++)
    {
        crntItm = *vit;
        cout<<setw(4)<<crntItm->index<<"\t"<<setw(4)<<crntItm->size<<"\t+"<<setprecision(4)<<crntItm->val<<"|"<<endl;
    }
    cout<<"--------------------------\n";
}

void Print2Scrn::printvect(double *v, const unsigned int n)
{
    cout<<"----------------------------"<<endl;
    for(unsigned int i = 0; i < n; i++)
    {
        cout<<v[i]<<" ";
    }
    cout<<"\n----------------------------"<<endl;
}


void Print2Scrn::printvect(vector<NNItem*> &vects, const char *dst_info_fn)
{
    assert(dst_info_fn);
    ofstream *outStrm = new ofstream(dst_info_fn, ios::out);
    vector<NNItem*>::iterator vit;
    NNItem* crntItm;
    (*outStrm)<<"Clabel\tSize\tISim/IDst\n";
    for(vit = vects.begin(); vit != vects.end(); vit++)
    {
        crntItm = *vit;
        (*outStrm)<<crntItm->index<<"\t"<<crntItm->size<<"\t"<<setprecision(4)<<crntItm->val<<endl;
    }
    (*outStrm).close();
}
