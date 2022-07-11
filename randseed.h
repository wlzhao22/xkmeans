#ifndef RANDSEED_H
#define RANDSEED_H

/***
*
*
*@author: Wan-Lei Zhao
*@date:   Apr.-1-2016
*
*
*produce random series of numbers
*
**/
#include <vector>

using namespace std;

class RandSeed
{
public:
    RandSeed(){}
    virtual ~RandSeed(){}
    static bool randTwo(unsigned int &r1, unsigned int &r2, vector<int> &vect);
    static bool randTwo(unsigned int &r1, unsigned int &r2, const unsigned int bound);
    static unsigned int  randOne(const unsigned int bound);
    static void test();
};

#endif // RANDSEED_H
