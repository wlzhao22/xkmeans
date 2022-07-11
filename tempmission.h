#ifndef TEMPMISSION_H
#define TEMPMISSION_H




class TempMission
{

public :

    static void datEvaluate(const char* datFn, const char *grdFn, float &e, float &p);
    static void datEvaluate();
    static void SEPLine(const char* filename, int t, int k);
    static void matlabTest();
    static void yzscore();
};








#endif
