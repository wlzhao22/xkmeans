#ifndef MISSIONAGENT_H
#define MISSIONAGENT_H

/**
* in charge of distributing the tasks
*
*
*@author:  Wan-Lei Zhao
*@date:    Apr.-12-2016
*
*
**/

#include <cstring>
#include <string>
#include <map>

using namespace std;

class MissionAgent
{
    public:
        MissionAgent(){}
        virtual ~MissionAgent(){}
        static bool buildClust(map<string,  const char*> &arguments);
        static void yingzhaoSolution(const char *fileName,const char *m, float *e, float *p, const char *crt, const char *prr, const char *sed, const char *rf);
        static void experiment();

};

#endif
