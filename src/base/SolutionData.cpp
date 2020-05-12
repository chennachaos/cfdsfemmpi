
#include "SolutionData.h"
#include "headersBasic.h"
#include "UtilitiesGeneral.h"


SolutionData::SolutionData()
{
  td.resize(100);

  tis = 0;
  rhoInf = 0.0;
  dt  = 1.0;

  return;
}



void SolutionData::initialise(int size1)
{
  vecsize = size1;

  soln.resize(vecsize);
  soln.setZero();

  solnInit  = soln;
  solnPrev  = soln;
  solnPrev2 = soln;
  solnPrev3 = soln;
  solnPrev4 = soln;

  solnDot     = soln;
  solnDotPrev = soln;
  solnDotCur  = soln;
  solnApplied = soln;

  return;
}





void SolutionData::setZero()
{
  soln.setZero();

  return;
}





void  SolutionData::setTimeParam()
{
  //cout << tis << '\t' << rhoInf << '\t' << dt << endl;
  SetTimeParametersFluid(tis, rhoInf, dt, td);
  //printVector(td);

   return;
}



void  SolutionData::timeUpdate()
{
  solnPrev4  = solnPrev3;
  solnPrev3  = solnPrev2;
  solnPrev2  = solnPrev;
  solnPrev   = soln;

  solnDotPrev  = solnDot;

  solnExtrap = solnPrev;
  //solnExtrap = 2.0*solnPrev - solnPrev2;

  return;
}




void SolutionData::updateIterStep()
{
  solnDot    = td[9]*soln + td[10]*solnPrev + td[11]*solnPrev2 + td[12]*solnPrev3 + td[13]*solnPrev4 + td[15]*solnDotPrev ;

  solnCur    = td[2]*soln    + (1.0-td[2])*solnPrev;
  solnDotCur = td[1]*solnDot + (1.0-td[1])*solnDotPrev;

  return;
}




void  SolutionData::reset()
{
  solnPrev  = solnPrev2;
  solnPrev2 = solnPrev3;
  solnPrev3 = solnPrev4;

  solnDot   = solnDotPrev;

  return;
}





