#ifndef SPEED_TEST
#define SPEED_TEST

class SpeedTest
{
    float *fdata;
    double *ddata;
    unsigned int num = 100000;
public :

    void generateData();
    void testmultiDouble();
    void testmultiFloat();
    void testmultiFloatDouble();
    static void test();
    virtual   ~SpeedTest();
};


#endif
