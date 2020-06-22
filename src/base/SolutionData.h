
#ifndef incl_SolutionData_h
#define incl_SolutionData_h

#include "headersBasic.h"
#include "headersEigen.h"


class  SolutionData
{
  public:
    int  vecsize;

    VectorXd  solnDot, solnDotPrev, solnDotCur, solnExtrap;
    VectorXd  soln, solnPrev, solnPrev2, solnPrev3, solnPrev4, solnCur, solnInit, solnApplied;
    VectorXd  td;

    vector<int>  node_map_get_old, node_map_get_new;

    SolutionData();

    ~SolutionData(){}

    void  initialise(int size1);

    void  setZero();

    void  setTimeParam(int tis, double rhoInf, double dt);

    void  timeUpdate();

    void  updateIterStep();

    void  reset();
};





#endif


