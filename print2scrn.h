#ifndef PRINT2SCRN_H
#define PRINT2SCRN_H

#include "clsinfo.h"
#include "nnitem.h"

#include <vector>

using namespace std;

class Print2Scrn
{
public:
    Print2Scrn() {}
    virtual ~Print2Scrn() {}

    static void   printvect(vector<NNItem*> &vects);
    static void   printvect(vector<NNItem*> &vects, const char *dst_info_fn);
    static void   printvect(double *v, const unsigned int n);
    static void   printCLInfo(const CLSInfo *info, const unsigned int n);
    static void   printCLInfo(const CLSInfo *info, const unsigned int n, const char *dstfn);
};

#endif // PRINT2SCRN_H
