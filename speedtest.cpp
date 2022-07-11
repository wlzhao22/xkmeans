#include "speedtest.h"

#include "timer.h"

#include <iostream>
#include <cstdlib>
#include <ctime>
using namespace std;

void SpeedTest::generateData()
{
    unsigned int i = 0;
    ddata = new double[num];
    fdata = new float[num];
    for(i = 0; i < num; i++)
    {
        ddata[i] = double(rand())/5*12;
        fdata[i] = float(rand())/5*12;
    }
}

void SpeedTest::testmultiDouble()
{
    unsigned int i = 0, j = 0;
    double tmp = 0;
    Timer *mytime = new Timer;
    mytime->start();
    for(j = 0; j < num; j++)
        for(i = 0; i < num; i++)
        {
            tmp += ddata[i]*ddata[i];
        }
    cout<<tmp<<" tp\n";
    cout<<"double\t***\n";
    mytime->end(true);
    delete mytime;
}


void SpeedTest::testmultiFloat()
{
    unsigned int i = 0, j = 0;
    double tmp = 0;
    Timer *mytime = new Timer;
    mytime->start();
    for(j = 0; j < num; j++)
        for(i = 0; i < num; i++)
        {
            tmp += fdata[i]*fdata[i];
        }
    cout<<tmp<<" tp\n";
    cout<<"float\t***\n";
    mytime->end(true);
    delete mytime;
}

void SpeedTest::testmultiFloatDouble()
{
    unsigned int i = 0, j = 0;
    double tmp = 0;
    Timer *mytime = new Timer;
    mytime->start();
    for(j = 0; j < num; j++)
        for(i = 0; i < num; i++)
        {
            tmp += fdata[i]*ddata[i];
        }
    cout<<tmp<<" tp\n";
    cout<<"float\tdouble***\n";


    mytime->end(true);
    delete mytime;
}

SpeedTest::~SpeedTest()
{
    delete [] ddata;
    delete [] fdata;
    ddata  = NULL;
    fdata = NULL;
}

void SpeedTest::test()
{
    srand(time(0));
    SpeedTest spt;
    spt.generateData();
    spt.testmultiFloat();
    spt.testmultiDouble();
    spt.testmultiFloatDouble();


}
