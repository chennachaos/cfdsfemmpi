
#ifndef KimMoinFlow_h
#define KimMoinFlow_h


#include <iostream>
//#include <limits.h>
//#include <float.h>
//#include <math.h>
//#include <cmath>

#include "headersBasic.h"

//using namespace Eigen;
//using namespace std;

using std::cout;
using std::endl;



//class  KimMoinFlow : public Function

/*
class  KimMoinFlow
{
  private:
    double rho, mu, a, c;

  public:
    KimMoinFlow(double r, double m, double c1)
    {
      rho = r;
      mu  = m;
      c   = c1;
      a   = 2.0;
    }

    virtual ~KimMoinFlow(){}

    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double tt=0.0)
    {
      double cx = cos(a*PI*xx);
      double sx = sin(a*PI*xx);
      double cy = cos(a*PI*yy);
      double sy = sin(a*PI*yy);
      double gu = exp(-2.0*tt);
      double gp = exp(-4.0*tt);

      //if(dir == 0)
        //return  ( 2.0*cx*sy*gu*(rho-mu*a*a*PI*PI) + a*PI*sx*cx*gp*(1.0-rho) );
      //else //if(dir == 1)
        //return  (-2.0*sx*cy*gu*(rho-mu*a*a*PI*PI) + a*PI*sy*cy*gp*(1.0-rho) );

      if(dir == 0)
        return  ( 2.0*cx*sy*gu*(rho-mu*a*a*PI*PI) + a*PI*sx*cx*gp );
      else //if(dir == 1)
        return  (-2.0*sx*cy*gu*(rho-mu*a*a*PI*PI) + a*PI*sy*cy*gp );
    }

    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double tt=0.0)
    {
      if(dir == 0)
        return ( -cos(a*PI*xx)*sin(a*PI*yy)*exp(-2.0*tt) );
      else if(dir == 1)
        return (  sin(a*PI*xx)*cos(a*PI*yy)*exp(-2.0*tt) );
      else //if(dir == 2)
        return  ( -0.25*(cos(2.0*a*PI*xx)+cos(2.0*a*PI*yy))*exp(-4.0*tt) );
    }

    virtual  void computeDerivatives(double xx, double yy, double tt, double* dv)
    {
      double c1 = exp(-2.0*tt);
      
      dv[0] =  a*PI*sin(a*PI*xx) * sin(a*PI*yy)*c1; // du/dx
      dv[2] = -a*PI*cos(a*PI*xx) * cos(a*PI*yy)*c1; // du/dy
      dv[1] =  a*PI*cos(a*PI*xx) * cos(a*PI*yy)*c1; // dv/dx
      dv[3] = -a*PI*sin(a*PI*xx) * sin(a*PI*yy)*c1; // dv/dy

      return;
    }
};
*/


/*
class  KimMoinFlow
{
  // Unsteady Stokes

  private:
    double rho, mu;

  public:
    KimMoinFlow(double r, double m, double c1)
    {
      rho = r;
      mu  = m;
    }

    virtual ~KimMoinFlow(){}

    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double tt=0.0)
    {
      double cx = cos(xx);
      double sx = sin(xx);
      double cy = cos(yy);
      double sy = sin(yy);

      if(dir == 0)
        return  (-2.0*cx*sy*(rho*cos(2.0*tt)+mu*sin(2.0*tt)) + sx*cx*sin(2.0*tt)*sin(2.0*tt) );
      else //if(dir == 1)
        return  ( 2.0*sx*cy*(rho*cos(2.0*tt)+mu*sin(2.0*tt)) + sy*cy*sin(2.0*tt)*sin(2.0*tt) );
    }

    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double tt=0.0)
    {
      if(dir == 0)
        return ( -cos(xx)*sin(yy)*sin(2.0*tt) );
      else if(dir == 1)
        return (  sin(xx)*cos(yy)*sin(2.0*tt) );
      else //if(dir == 2)
        return  ( -0.25*(cos(2.0*xx)+cos(2.0*yy))*sin(2.0*tt)*sin(2.0*tt) );
    }

    virtual  void computeDerivatives(double xx, double yy, double tt, double* dv)
    {
      double c1 = sin(2.0*tt);
      
      dv[0] =  sin(xx)*sin(yy)*c1; // du/dx
      dv[2] = -cos(xx)*cos(yy)*c1; // du/dy
      dv[1] =  cos(xx)*cos(yy)*c1; // dv/dx
      dv[3] = -sin(xx)*sin(yy)*c1; // dv/dy

      return;
    }
};
*/


/*
class  KimMoinFlow
{
  // For Steady Stokes  
  
  private:
    double rho, mu, a, c;

  public:
    KimMoinFlow(double r, double m, double c1)
    {
      rho = r;
      mu  = m;
    }

    virtual ~KimMoinFlow(){}

    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double tt=0.0)
    {
      double cx = cos(xx);
      double sx = sin(xx);
      double cy = cos(yy);
      double sy = sin(yy);

      if(dir == 0)
        return  (-2.0*mu*cx*sy + sx*cx );
      else //if(dir == 1)
        return  ( 2.0*mu*sx*cy + sy*cy );
    }

    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double tt=0.0)
    {
      if(dir == 0)
        return ( -cos(xx)*sin(yy) );
      else if(dir == 1)
        return (  sin(xx)*cos(yy) );
      else //if(dir == 2)
        return  ( -0.25*(cos(2.0*xx)+cos(2.0*yy)) );
    }

    virtual  void computeDerivatives(double xx, double yy, double tt, double* dv)
    {
      dv[0] =  sin(xx)*sin(yy); // du/dx
      dv[2] = -cos(xx)*cos(yy); // du/dy
      dv[1] =  cos(xx)*cos(yy); // dv/dx
      dv[3] = -sin(xx)*sin(yy); // dv/dy

      return;
    }
};
*/


/*
class  KimMoinFlow
{
  // For Steady Navier-Stokes  
  
  private:
    double rho, mu;

  public:
    KimMoinFlow(double r, double m, double c1)
    {
      rho = r;
      mu  = m;
    }

    virtual ~KimMoinFlow(){}

    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double tt=0.0)
    {
      double cx = cos(xx);
      double sx = sin(xx);
      double cy = cos(yy);
      double sy = sin(yy);

      if(dir == 0)
        return  (-rho*sx*cx - 2.0*mu*cx*sy + sx*cx );
      else //if(dir == 1)
        return  (-rho*sy*cy + 2.0*mu*sx*cy + sy*cy );
    }

    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double tt=0.0)
    {
      if(dir == 0)
        return ( -cos(xx)*sin(yy) );
      else if(dir == 1)
        return (  sin(xx)*cos(yy) );
      else //if(dir == 2)
        return  ( -0.25*(cos(2.0*xx)+cos(2.0*yy)) );
    }

    virtual  void computeDerivatives(double xx, double yy, double tt, double* dv)
    {
      dv[0] =  sin(xx)*sin(yy); // du/dx
      dv[2] = -cos(xx)*cos(yy); // du/dy
      dv[1] =  cos(xx)*cos(yy); // dv/dx
      dv[3] = -sin(xx)*sin(yy); // dv/dy

      return;
    }
};
*/


class  KimMoinFlow
{
  // Unsteady Navier-Stokes

  private:
    double rho, mu;

  public:
    KimMoinFlow(double r, double m, double c1)
    {
      rho = r;
      mu  = m;
    }

    virtual ~KimMoinFlow(){}

    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double tt=0.0)
    {
      double cx = cos(xx);
      double sx = sin(xx);
      double cy = cos(yy);
      double sy = sin(yy);

      //if(dir == 0)
        //return  (-2.0*cx*sy*(rho*cos(2.0*tt)+mu*sin(2.0*tt)) - rho*sx*cx*sin(2.0*tt)*sin(2.0*tt) + sx*cx*sin(2.0*tt)*sin(2.0*tt) );
      //else //if(dir == 1)
        //return  ( 2.0*sx*cy*(rho*cos(2.0*tt)+mu*sin(2.0*tt)) - rho*sy*cy*sin(2.0*tt)*sin(2.0*tt) + sy*cy*sin(2.0*tt)*sin(2.0*tt) );

      if(dir == 0)
        return ( sx*cx*sin(2.0*tt)*sin(2.0*tt) - 2.0*rho*cos(2.0*tt)*cx*sy - 2.0*mu*sin(2.0*tt)*cx*sy );
      else //if(dir == 1)
        return ( sy*cy*sin(2.0*tt)*sin(2.0*tt) + 2.0*rho*cos(2.0*tt)*sx*sy + 2.0*mu*sin(2.0*tt)*sx*cy );
    }

    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double tt=0.0)
    {
      if(dir == 0)
        return ( -cos(xx)*sin(yy)*sin(2.0*tt) );
      else if(dir == 1)
        return (  sin(xx)*cos(yy)*sin(2.0*tt) );
      else //if(dir == 2)
        return  ( -0.25*(cos(2.0*xx)+cos(2.0*yy))*sin(2.0*tt)*sin(2.0*tt) );
    }

    virtual  void computeDerivatives(double xx, double yy, double tt, double* dv)
    {
      double c1 = sin(2.0*tt);

      dv[0] =  sin(xx)*sin(yy)*c1; // du/dx
      dv[2] = -cos(xx)*cos(yy)*c1; // du/dy
      dv[1] =  cos(xx)*cos(yy)*c1; // dv/dx
      dv[3] = -sin(xx)*sin(yy)*c1; // dv/dy

      return;
    }
};



class  Kovasznay
{
  private:
    double  mu, Re, lamda, p0;

  public:
    Kovasznay()
    {
      Re = 40.0;
      mu = 1.0/Re;
      lamda = 0.5*Re - sqrt(0.25*Re*Re + 4.0*PI*PI);
      p0 = 0.5*exp(-lamda);
    }

    Kovasznay(double Re1)
    {
      Re = Re1;

      lamda = 0.5*Re - sqrt(0.25*Re*Re + 4.0*PI*PI);

      p0 = 0.5*exp(-lamda);
    }

    virtual ~Kovasznay(){}

    void SetPressure(double pp)
    {
      p0 = pp;
    }
    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      // for Navier-Stokes the body force is zero
      return  0.0;

      //if(dir == 0)
        //return  (-mu*(4.0*PI*PI-lamda*lamda)*exp(lamda*xx)*cos(2.0*PI*yy)-lamda*exp(2.0*lamda*xx));
      //else
        //return  (-mu*(lamda/2.0/PI)*(-4.0*PI*PI+lamda*lamda)*exp(lamda*xx)*sin(2.0*PI*yy));
    }

    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      if(dir == 0)
        return ( 1.0 - exp(lamda*xx) * cos(2.0*PI*yy) );
      else if(dir == 1)
        return ( (lamda/2.0/PI) * exp(lamda*xx) * sin(2.0*PI*yy) );
      else //if(dir == 2)
        return  (p0 - 0.5*exp(2.0*lamda*xx));
    }

    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      double  fact = exp(lamda*xx), fact2 = 2.0*PI*yy;;

      dv[0] = - lamda * fact * cos(fact2);
      dv[2] = 2.0*PI*fact* sin(fact2);
      dv[1] = (lamda*lamda/2.0/PI) * fact * sin(fact2);
      dv[3] = lamda*fact* cos(fact2);

      return;
    }
};




class  Stokes2DEx1
{
  public:

    Stokes2DEx1() {}

    virtual ~Stokes2DEx1() {}

    double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      if(dir == 0)
        return  ( (12.0*(1.0-2.0*yy)*(xx-2.0)*pow(xx,3.0)) + ((12.0-48.0*yy+72.0*yy*yy-48.0*yy*yy*yy)*pow(xx,2.0)) +  ((-2.0+24.0*yy-72*yy*yy+48*yy*yy*yy)*xx) + (1.0-4.0*yy+12.0*yy*yy-8.0*yy*yy*yy) );
      else //if(dir == 1)
        return  ( ((8.0-48.0*yy+48.0*yy*yy)*pow(xx,3.0)) + ((-12.0+72.0*yy-72.0*yy*yy)*pow(xx,2.0)) + ((4.0-24.0*yy+48*yy*yy-48*yy*yy*yy+24.0*yy*yy*yy*yy)*xx) + (-12.0*yy*yy+24.0*yy*yy*yy-12.0*yy*yy*yy*yy) );
    }

    double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      if(dir == 0)
        return  ( pow((xx - xx*xx),2.0)*(2.0*yy - 6.0*yy*yy + 4.0*yy*yy*yy) );
      else if(dir == 1)
        return  ( -pow((yy - yy*yy),2.0)*(2.0*xx - 6.0*xx*xx + 4.0*xx*xx*xx) );
      else //if(dir == 2)
        return  (xx*(1.0-xx));
    }

    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      dv[0] = (2.0*xx-6.0*xx*xx+4.0*xx*xx*xx)*(2.0*yy-6.0*yy*yy+4.0*yy*yy*yy);
      dv[2] = (xx*xx-2.0*xx*xx*xx+xx*xx*xx*xx)*(2.0-12.0*yy+12.0*yy*yy);
      dv[1] = -(2.0-12.0*xx+12.0*xx*xx)*(yy*yy-2.0*yy*yy*yy+yy*yy*yy*yy);
      dv[3] = -(2.0*xx-6.0*xx*xx+4.0*xx*xx*xx)*(2.0*yy-6.0*yy*yy+4.0*yy*yy*yy);

      return ;
    }
};



class  Stokes2DEx2
{
  public:

    double  nu, a, b;

    Stokes2DEx2(){nu=1.0;a=2.0*PI; b=3.0*PI;}

    virtual ~Stokes2DEx2() {}

    double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      if(dir == 0)
        return  ( 2.0*nu*a*a*sin(a*xx)*cos(a*yy) + b*cos(b*xx)*cos(b*yy));
      else //if(dir == 1)
        return  (-2.0*nu*a*a*cos(a*xx)*sin(a*yy) - b*sin(b*xx)*sin(b*yy));
    }

    double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      if(dir == 0)
        return  ( sin(a*xx)*cos(a*yy) );
      else if(dir == 1)
        return  ( -cos(a*xx)*sin(a*yy) );
      else //if(dir == 2)
        return  (sin(b*xx)*cos(b*yy));
    }

    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      dv[0] =  a*cos(a*xx)*cos(a*yy); // du/dx
      dv[2] = -a*sin(a*xx)*sin(a*yy); // du/dy
      dv[1] =  a*sin(a*xx)*sin(a*yy); // dv/dx
      dv[3] = -a*cos(a*xx)*cos(a*yy); // dv/dy

      return ;
    }
};



class  PoissonEx3
{
  public:
    PoissonEx3() {}

    virtual ~PoissonEx3() {}

    virtual double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  (cosh(PI*yy) - sinh(PI*yy)/tanh(PI))*sin(PI*xx);
    }
    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  0.0;
    }
    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      dv[0] = (cosh(PI*yy) - sinh(PI*yy)/tanh(PI))*PI*cos(PI*xx);
      dv[1] = (PI*sinh(PI*yy) - PI*cosh(PI*yy)/tanh(PI))*sin(PI*xx);

      return ;
    }
};




class  PoissonAnnulus
{
  public:
    PoissonAnnulus() {}

    virtual ~PoissonAnnulus() {}

    virtual double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  ( (2.0*(xx*xx + yy*yy - 1.0)*yy) / (3.0*(xx*xx+yy*yy)) );
    }
    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return 0.0;
    }
    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      double r2 = xx*xx+yy*yy;
      dv[0] = (4.0*xx*yy)/ (3.0*r2*r2);
      dv[1] = (2.0/(3.0*r2*r2))*(r2*r2 - xx*xx + yy*yy);

      return ;
    }
};




class  KimMoinFlowUnsteadyNavierStokes
{
  // Unsteady Navier-Stokes

  private:
    double rho, mu;

  public:
    KimMoinFlowUnsteadyNavierStokes(double r, double m, double c1)
    {
      rho = r;
      mu  = m;
    }

    virtual ~KimMoinFlowUnsteadyNavierStokes(){}

    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0, double tt=0.0)
    {
      double cx = cos(xx);
      double sx = sin(xx);
      double cy = cos(yy);
      double sy = sin(yy);
      double s2t = sin(2.0*tt);
      double c2t = cos(2.0*tt);

      if(dir == 0)
        return  (-2.0*cx*sy*(rho*c2t+mu*s2t) - rho*sx*cx*s2t*s2t + sx*cx*s2t*s2t );
      else //if(dir == 1)
        return  ( 2.0*sx*cy*(rho*c2t+mu*s2t) - rho*sy*cy*s2t*s2t + sy*cy*s2t*s2t );
    }

    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0, double tt=0.0)
    {
      if(dir == 0)
        return ( -cos(xx)*sin(yy)*sin(2.0*tt) );
      else if(dir == 1)
        return (  sin(xx)*cos(yy)*sin(2.0*tt) );
      else //if(dir == 2)
        return  ( -0.25*(cos(2.0*xx)+cos(2.0*yy))*sin(2.0*tt)*sin(2.0*tt) );
    }

    virtual  void computeDerivatives(double xx, double yy, double zz, double tt, double* dv)
    {
      double s2t = sin(2.0*tt);

      dv[0] =  sin(xx)*sin(yy)*s2t; // du/dx
      dv[2] = -cos(xx)*cos(yy)*s2t; // du/dy
      dv[1] =  cos(xx)*cos(yy)*s2t; // dv/dx
      dv[3] = -sin(xx)*sin(yy)*s2t; // dv/dy

      return;
    }
};




#endif


