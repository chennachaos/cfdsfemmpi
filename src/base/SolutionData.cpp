
#include "SolutionData.h"
#include "headersBasic.h"


SolutionData::SolutionData()
{
  td.resize(100);

  tis = 0;
  rho = 0.0;
  dt  = 1.0;

  return;
}



void SolutionData::initialise()
{
  int ind=nNode_Velo*ndim;

  velo.resize(ind);
  velo.setZero();

  veloPrev  = velo;
  veloCur   = velo;
  veloPrev2 = velo;
  veloPrev3 = velo;

  veloDot     = velo;
  veloDotPrev = veloDot;
  acceCur  = veloDot;

  pres.resize(nNode_Pres);
  pres.setZero();

  presPrev  = pres;
  presPrev2 = pres;
  presPrev3 = pres;
  presCur   = pres;

  presDot     = velo;
  presDotPrev = pres;
  presDotCur  = pres;

  //velo.setZero();      veloPrev.setZero();      veloCur.setZero();
  //veloDot.setZero();      veloDotPrev.setZero();      acceCur.setZero();
  //pres.setZero();      presPrev.setZero();      presCur.setZero();
  //presDot.setZero();   presDotPrev.setZero();   presDotCur.setZero();

  return;
}





void  SolutionData::setTimeParam()
{
  //cout << PHYSICS_TYPE << endl;
  //cout << tis << '\t' << rho << '\t' << mpapTime.dt << endl;  
  //printVector(td);

   // solve for the initial profile

   return;
}



void  SolutionData::timeUpdate()
{
/*
        veloPrev3  = veloPrev2;
        veloPrev2  = veloPrev;
        veloPrev   = velo;

        veloDotPrev  = veloDot;

        presPrev3  = presPrev2;
        presPrev2  = presPrev;
        presPrev   = pres;
        presDotPrev  = presDot;



  // store the variables 

  var1Prev4 = var1Prev3;
  var1Prev3 = var1Prev2;
  var1Prev2 = var1Prev;

  var1Prev = var1;
  var2Prev = var2;
  var3Prev = var3;
  var4Prev = var4;
  
  var1Extrap = var1Prev;
  //var1Extrap = 2.0*var1Prev - var1Prev2;
  //var1Extrap = 3.0*var1Prev - 3.0*var1Prev2 + var1Prev3;

  var1DotPrev = var1Dot;
  var1DotDotPrev = var1DotDot;
  
  dDotPrev = dDot;

  var2DotPrev = var2Dot;
  var2DotDotPrev = var2DotDot;

  bodyForcePrev = bodyForce;

  // initialise the variables

  forcePrev5 = forcePrev4;
  forcePrev4 = forcePrev3;
  forcePrev3 = forcePrev2;
  forcePrev2 = forcePrev;
  forcePrev  = force;
*/
  return;
}




void SolutionData::updateIterStep()
{
/*
  char fct[] = "SolutionData::updateIterStep";

  if( PHYSICS_TYPE == PHYSICS_TYPE_SOLID )
  {
    //if(STAGGERED)
    //{
      // displacement as the primary variable
      //cout << "  SOLID SOLID SOLID " << endl;

      var1Dot     = td[10]*var1 + td[11]*var1Prev + td[12]*var1DotPrev + td[13]*var1DotDotPrev + td[14]*dDotPrev;
      var1DotDot  = td[15]*var1 + td[16]*var1Prev + td[17]*var1DotPrev + td[18]*var1DotDotPrev + td[19]*dDotPrev;
      dDot        = td[20]*var1 + td[21]*var1Prev + td[22]*var1DotPrev + td[23]*var1DotDotPrev + td[24]*dDotPrev;
    //}
    //else
    //{
      // velocity as the primary variable

      //var1       = td[40]*var1Dot + td[41]*var1Prev + td[42]*var1DotPrev + td[43]*var1DotDotPrev + td[44]*dDotPrev;
      //var1DotDot = td[45]*var1Dot + td[46]*var1Prev + td[47]*var1DotPrev + td[48]*var1DotDotPrev + td[49]*dDotPrev;
      //// ddot_{n+1} for modified state-space formulation
      //dDot       = td[50]*var1Dot + td[51]*var1Prev + td[52]*var1DotPrev + td[53]*var1DotDotPrev + td[54]*dDotPrev;
    //}

    // compute Current values

    var1Cur       = td[2]*var1       + (1.0-td[2])*var1Prev;
    var1DotCur    = td[2]*var1Dot    + (1.0-td[2])*var1DotPrev;
    var1DotDotCur = td[1]*var1DotDot + (1.0-td[1])*var1DotDotPrev;

    //printVector(var1Cur);
    //printVector(var1DotCur);
    //printVector(var1DotDotCur);

      //cout << " ggggggggggggg " << endl;
  }
  else
  {
    //cout << " uuuuuuuuuuuu  " << endl;
    //printVector(td);
    var1Dot    = td[9]*var1 + td[10]*var1Prev + td[11]*var1Prev2 + td[12]*var1Prev3 + td[13]*var1Prev4 + td[15]*var1DotPrev ;
    var2Dot    = td[9]*var2 + td[10]*var2Prev + td[15]*var2DotPrev ;

    var1Cur    = td[2]*var1    + (1.0-td[2])*var1Prev; // velocity
    var2Cur    = td[2]*var2    + (1.0-td[2])*var2Prev; // pressure
    var3Cur    = td[2]*var3    + (1.0-td[2])*var3Prev; // Lagrange parameters
    var4Cur    = td[2]*var4    + (1.0-td[2])*var4Prev; // Solid dof (for contact elements)

    var1DotCur = td[1]*var1Dot + (1.0-td[1])*var1DotPrev;
    var2DotCur = td[1]*var2Dot + (1.0-td[1])*var2DotPrev;
    //cout << " uuuuuuuuuuuu  " << endl;
  }
*/
  return;
}


void  SolutionData::reset()
{
/*
  var1 = var1Prev;
  var2 = var2Prev;
  var3 = var3Prev;
  var4 = var4Prev;

  var1Prev  = var1Prev2;
  var1Prev2 = var1Prev3;
  var1Prev3 = var1Prev4;

  var1Dot = var1DotPrev ;
  var1DotDot = var1DotDotPrev ;
  
  dDot = dDotPrev;

  var2Dot = var2DotPrev ;
  var2DotDot = var2DotDotPrev ;

  bodyForce = bodyForcePrev;
  bodyForcePrev = bodyForcePrev2;

  force = forcePrev;
  forcePrev = forcePrev2;
  forcePrev2 = forcePrev3;
  forcePrev3 = forcePrev4;
  forcePrev4 = forcePrev5;
*/
  return;
}




void  SolutionData::applyDirichletBCs(double fact)
{
  int  ii, jj, n1, n2;
  double  gamm;

        //Add specified Dirichlet BC
        for(ii=0; ii<nDBC; ii++)
        {
          n1 = DirichletBCs[ii][0];
          n2 = DirichletBCs[ii][1];

          if( n2 < 2)
          {
            jj = n1*2+n2;

            velo(jj) = DirichletBCs[ii][2]*fact;

            //xx = node_coords[n1][0];
            //yy = node_coords[n1][1];
            //if(n2 == 0)
              //velo(jj) = -cos(xx)*sin(yy)*sin(2.0*timeNow) ;
            //else
              //velo(jj) =  sin(xx)*cos(yy)*sin(2.0*timeNow) ;

            veloDot(jj) = (velo(jj)-veloPrev(jj))/(gamm*dt) - (1.0-gamm)*veloDotPrev(jj)/gamm;
          }
        }
  return;
}


void  SolutionData::addNodalForces(double fact)
{
  int  ii, jj, n1, n2;
  double  gamm;
        // Add specified nodal force 
        //for(ii=0; ii<nFBC; ii++)
        //{
          //n1 = ForceBCs[ii][0];
          //n2 = ForceBCs[ii][1];

          //row  = n1*2+n2;
          //fact = ForceBCs[ii][2];

          //rhsVecVelo(row) = rhsVecVelo(row) + fact;
        //}
  return;
}






