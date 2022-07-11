#ifndef CLSINDEX_H
#define CLSINDEX_H

#include <vector>

using std::vector;

struct CLSInfo
{
public:
    CLSInfo()
    {
        n = 0;
        E = 0.0;
    }
    unsigned int n;
    double E;
};

struct CLData : public CLSInfo
{
    vector<unsigned int > clustDataId;
};


#endif
