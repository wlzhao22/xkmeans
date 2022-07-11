#include "randseed.h"
#include "pqmath.h"

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cmath>
#include <ctime>

using namespace std;

bool RandSeed::randTwo(unsigned int &r1, unsigned int &r2, vector<int> &vect)
{
    const unsigned int bound   = vect.size();
    unsigned int t = 0;
    r1 = (unsigned int)floor((rand()/(RAND_MAX+1.0f))*bound);
    r1 = (r1 >= bound)?(bound-1):r1;
    r2 = (unsigned int)floor((rand()/(RAND_MAX+1.0f))*bound);
    r2 = (r2 >= bound)?(bound-1):r2;

    while(r1 == r2 && t < 5000)
    {
        r2 = (unsigned int)floor((rand()/(RAND_MAX+1.0f))*bound);
        r2 = (r2 >= bound)?(bound-1):r2;
        t++;
    }
    assert(r1 != r2);
    r1 = vect[r1];
    r2 = vect[r2];
    return true;
}

bool RandSeed::randTwo(unsigned int &r1, unsigned int &r2, const unsigned int bound)
{
    unsigned int t = 0;
    r1 = (unsigned int)floor((rand()/(RAND_MAX+1.0f))*bound);
    r1 = (r1 >= bound)?(bound-1):r1;
    r2 = (unsigned int)floor((rand()/(RAND_MAX+1.0f))*bound);
    r2 = (r2 >= bound)?(bound-1):r2;

    while(r1 == r2 && t < 5000)
    {
        r2 = (unsigned int)floor((rand()/(RAND_MAX+1.0f))*bound);
        r2 = (r2 >= bound)?(bound-1):r2;
        t++;
    }
    assert(r1 != r2);
    return true;
}

unsigned int RandSeed::randOne(const unsigned int bound)
{
    unsigned int sed0 = time(NULL);
    unsigned int r1   = (unsigned int)floor((rand_r(&sed0)/(RAND_MAX+1.0f))*bound);
    ///unsigned int r1 = floor((rand()/(RAND_MAX+1.0f))*bound);
    cout<<"r1: "<<sed0<<"\t"<<r1<<endl;

    return r1;
}

void RandSeed::test()
{
    unsigned int i = time(NULL);
    unsigned int k = RandSeed::randOne(30);
    cout<<i<<"\t"<<k<<endl;
}
