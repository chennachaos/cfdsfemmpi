
#ifndef incl_SolutionData_h
#define incl_SolutionData_h

#include "headersBasic.h"
#include "headersEigen.h"


class  SolutionData
{
  public:
    int  vecsize, tis;
    double  rhoInf, dt;

    VectorXd  solnDot, solnDotPrev, solnDotCur, solnExtrap;
    VectorXd  soln, solnPrev, solnPrev2, solnPrev3, solnPrev4, solnCur, solnInit, solnApplied;
    VectorXd  td;

    vector<int>  node_map_get_old, node_map_get_new;
    vector<int>  dof_map_get_old, dof_map_get_new;


    SolutionData();

    ~SolutionData(){}

    void setTimeIncrementType(int ttt)
    {  tis = ttt; }

    void setSpectralRadius(double ttt)
    {  rhoInf = ttt; }

    void  initialise(int size1);

    void  setZero();

    void  setTimeParam();

    void  timeUpdate();

    void  updateIterStep();

    void  reset();
};





#endif


