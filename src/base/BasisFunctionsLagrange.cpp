
#include "BasisFunctionsLagrange.h"
#include "UtilitiesGeneral.h"
#include "headersBasic.h"
#include "headersEigen.h"
#include "ElementBase.h"

#include <math.h>
#include <cmath>
#include <stdio.h>
#include <vector>


using namespace std;



int computeBasisFunctions2D(bool flag, int ELEM_TYPE, int degree, double* param, double* xNode, double* yNode, double*  N, double*  dN_dx, double* dN_dy, double&  Jac)
{
    int  ii=0, jj=0, count=0, nlbf=0;

    if( ELEM_TYPE == ELEM_SHAPE_TRIA )
    {
      if(degree == 1)
      {
        count = 3;
        nlbf = count;
      }
      else if(degree == 2)
      {
        count = 6;
        nlbf = count;
      }
    }
    else if( ELEM_TYPE == ELEM_SHAPE_QUAD )
    {
      count = degree + 1;
      nlbf = count*count;
    }
    else
    {
      cerr << " Invalid Element Type in computeBasisFunctions2D... " << endl;
      exit(-2);
    }

    vector<double>  N1(count), N2(count), dN1(count), dN2(count);
    vector<double>  dN_du1(nlbf), dN_du2(nlbf);
    double  xx, yy, detinv, B[2][2], Binv[2][2] ;

    if( ELEM_TYPE == ELEM_SHAPE_TRIA )
    {
      LagrangeBasisFunsTria(degree, param[0], param[1], &N[0], &dN_du1[0], &dN_du2[0]);
    }
    else if( ELEM_TYPE == ELEM_SHAPE_QUAD )
    {
      LagrangeBasisFunsQuad(degree, param[0], param[1], &N[0], &dN_du1[0], &dN_du2[0]);
    }
    else
    {
      cerr << " Invalid Element Type in computeBasisFunctions2D... " << endl;
      exit(-2);
    }

    // Gradient of mapping from parameter space to physical space

    B[0][0] = B[1][0] = B[0][1] = B[1][1] = 0.0;
    for(ii=0; ii<nlbf; ii++)
    {
      xx = xNode[ii];
      yy = yNode[ii];

      B[0][0] +=  (xx * dN_du1[ii]) ;
      B[1][0] +=  (xx * dN_du2[ii]) ;
      B[0][1] +=  (yy * dN_du1[ii]) ;
      B[1][1] +=  (yy * dN_du2[ii]) ;
    }

    Jac  = B[0][0]*B[1][1] - B[0][1]*B[1][0];

    //printf("Bmat  \t%20.18f\t%20.18f\t%20.18f\t%20.18f\t%20.18f \n\n", B[0][0], B[0][1], B[1][0], B[1][1], Jac);

    detinv = 1.0/Jac ;

    Binv[0][0] =  B[1][1] * detinv;
    Binv[0][1] = -B[0][1] * detinv;
    Binv[1][0] = -B[1][0] * detinv;
    Binv[1][1] =  B[0][0] * detinv;

    // Compute derivatives of basis functions w.r.t physical coordinates
    for(ii=0; ii<nlbf; ii++)
    {
      dN_dx[ii] = dN_du1[ii] * Binv[0][0] + dN_du2[ii] * Binv[0][1];
      dN_dy[ii] = dN_du1[ii] * Binv[1][0] + dN_du2[ii] * Binv[1][1];
    }

    return 0;
}



int computeBasisFunctions3D(bool flag, int ELEM_TYPE, int degree, double* param, double* xNode, double* yNode, double* zNode, double*  N, double*  dN_dx, double* dN_dy, double* dN_dz, double&  Jac)
{
    int  ii=0, jj=0, count=0, nlbf=0;

    if( ELEM_TYPE == ELEM_SHAPE_TETRA )
    {
      if(degree == 1)
      {
        count = 4;
        nlbf = count;
      }
      else if(degree == 2)
      {
        count = 10;
        nlbf = count;
      }
    }
    else if(ELEM_TYPE == ELEM_SHAPE_HEXA)
    {
      count = degree + 1;
      nlbf = count*count*count;
    }
    else
    {
      cerr << " Invalid Element Type in computeBasisFunctions3D... " << endl;
      exit(-2);
    }

    vector<double>  N1(count), N2(count), dN1(count), dN2(count), dN3(count);
    vector<double>  dN_du1(nlbf), dN_du2(nlbf), dN_du3(nlbf);
    double  xx, yy, zz, detinv;
    MatrixXd B(3,3), Binv(3,3) ;

    if( ELEM_TYPE == ELEM_SHAPE_TETRA)
    {
      LagrangeBasisFunsTetra(degree, param[0], param[1], param[2], &N[0], &dN_du1[0], &dN_du2[0], &dN_du3[0]);
    }
    else if(ELEM_TYPE == ELEM_SHAPE_HEXA)
    {
      LagrangeBasisFunsHexa(degree, param[0], param[1], param[2], &N[0], &dN_du1[0], &dN_du2[0], &dN_du3[0]);
    }
    else
    {
      cerr << " Invalid Element Type in computeBasisFunctions3D... " << endl;
      exit(-2);
    }

    // Gradient of mapping from parameter space to physical space

    B.setZero();
    for(ii=0; ii<nlbf; ii++)
    {
      xx = xNode[ii];
      yy = yNode[ii];
      zz = zNode[ii];

      B(0,0) +=  (xx * dN_du1[ii]) ;
      B(1,0) +=  (xx * dN_du2[ii]) ;
      B(2,0) +=  (xx * dN_du3[ii]) ;

      B(0,1) +=  (yy * dN_du1[ii]) ;
      B(1,1) +=  (yy * dN_du2[ii]) ;
      B(2,1) +=  (yy * dN_du3[ii]) ;

      B(0,2) +=  (zz * dN_du1[ii]) ;
      B(1,2) +=  (zz * dN_du2[ii]) ;
      B(2,2) +=  (zz * dN_du3[ii]) ;
    }

    //printMatrix(B);

    Jac  = B.determinant();
    Binv = B.inverse();

    //printMatrix(Binv);

    // Compute derivatives of basis functions w.r.t physical coordinates
    for(ii=0; ii<nlbf; ii++)
    {
      dN_dx[ii] = dN_du1[ii] * Binv(0,0) + dN_du2[ii] * Binv(0,1) + dN_du3[ii] * Binv(0,2);
      dN_dy[ii] = dN_du1[ii] * Binv(1,0) + dN_du2[ii] * Binv(1,1) + dN_du3[ii] * Binv(1,2);
      dN_dz[ii] = dN_du1[ii] * Binv(2,0) + dN_du2[ii] * Binv(2,1) + dN_du3[ii] * Binv(2,2);
    }

    return 0;
}




void Lagrange_BasisFuns1D(int p, double xi, double* N, double* dN_dxi, double* d2N_dxi2)
{
  double  fact1, fact2, val1, val2, val3, val4;

  switch(p)
  {
      case 0:

          N[0] = 1.0;

          dN_dxi[0] = 0.0;

          d2N_dxi2[0] = 0.0;

      break;

      case 1:

          N[0] = 0.5*(1.0 - xi);
          N[1] = 0.5*(1.0 + xi);

          dN_dxi[0] = -0.5;
          dN_dxi[1] =  0.5;

          d2N_dxi2[0] = 0.0;
          d2N_dxi2[1] = 0.0;

      break;

      case 2:

          val1 = xi*xi;

          N[0] = -0.5*(xi - val1);
          N[1] =  1.0 - val1;
          N[2] =  0.5*(xi + val1);

          val1 = 2.0*xi;

          dN_dxi[0] = -0.5*(1.0 - val1);
          dN_dxi[1] = -val1;
          dN_dxi[2] =  0.5*(1.0 + val1);

          d2N_dxi2[0] =  1.0;
          d2N_dxi2[1] = -2.0;
          d2N_dxi2[2] =  1.0;

      break;

      case 3:

          fact1 = 9.0/16.0;
          fact2 = 27.0/16.0;
          val1  = xi*xi;

          N[0] = -fact1 * (1 - xi)   * (1/9 - val1);
          N[1] =  fact2 * (1 - val1) * (1/3 - xi);
          N[2] =  fact2 * (1 - val1) * (1/3 + xi);
          N[3] = -fact1 * (1 + xi)   * ( 1/9 - val1);

          val2 = 3.0*val1;

          dN_dxi[0] = -fact1*(-1/9 - 2.0*xi   +  val2);
          dN_dxi[1] =  fact2*(-1   - 2.0/3*xi +  val2);
          dN_dxi[2] =  fact2*(1    - 2.0/3*xi -  val2);
          dN_dxi[3] = -fact1*(1/9  - 2.0*xi   -  val2);

          val2 = 6.0*xi;

          d2N_dxi2[0] = -fact1 * (-2   + val2);
          d2N_dxi2[1] =  fact2 * (-2/3 + val2);
          d2N_dxi2[2] =  fact2 * (-2/3 - val2);
          d2N_dxi2[3] = -fact1 * (-2   - val2);

      break;

      case 4:

          fact1 = 2.0/3.0;
          fact2 = 8.0/3.0;
          val1 = xi*xi;
          val2 = val1*xi;
          val3 = val2*xi;

          N[0] =  fact1 * (0.25*xi  - 0.25*val1  -  val2     + val3);
          N[1] = -fact2 * (0.5 *xi  - val1       -  0.5*val2 + val3);
          N[2] =    4.0 * (0.25     - 1.25*val1  -  0        + val3);
          N[3] =  fact2 * (0.5*xi   + val1       -  0.5*val2 - val3);
          N[4] = -fact1 * (0.25*xi  + 0.25*val1  -  val2     - val3);

          val4 = 4.0*val2;

          dN_dxi[0] =  fact1 * (0.25 - 0.5*xi  - 3.0*val1  + val4);
          dN_dxi[1] = -fact2 * (0.5  - 2.0*xi  - 1.5*val1  + val4);
          dN_dxi[2] =    4.0 * (  0  - 2.5*xi  -   0.0     + val4);
          dN_dxi[3] =  fact2 * (0.5  + 2.0*xi  - 1.5*val1  - val4);
          dN_dxi[4] = -fact1 * (0.25 + 0.5*xi  - 3.0*val1  - val4);

          val4 = 12.0*val1;

          d2N_dxi2[0] =  fact1 * (-0.5  -  6.0*xi  +  val4);
          d2N_dxi2[1] = -fact2 * (-2.0  -  3.0*xi  +  val4);
          d2N_dxi2[2] =    4.0 * (-2.5  -  0.0     +  val4);
          d2N_dxi2[3] =  fact2 * ( 2.0  -  3.0*xi  -  val4);
          d2N_dxi2[4] = -fact1 * ( 0.5  -  6.0*xi  -  val4);

      break;

      default:

          printf("no basis functions defined for this degree = %5d \n", p);

      break;
  }

   return;
}



void Lagrange_BasisFuns1D(int p, double xi, double* N, double* dN_dxi)
{
  double  fact1, fact2, val1, val2, val3, val4;

  switch(p)
  {
      case 0:

          N[0] = 1.0;

          dN_dxi[0] = 0.0;

      break;

      case 1:

          N[0] = 0.5*(1.0 - xi);
          N[1] = 0.5*(1.0 + xi);

          dN_dxi[0] = -0.5;
          dN_dxi[1] =  0.5;

      break;

      case 2:

          val1 = xi*xi;

          N[0] = -0.5 * (xi - val1);
          N[1] = 1.0 - val1;
          N[2] = 0.5 *(xi + val1);

          val1 = 2.0*xi;

          dN_dxi[0] = -0.5*(1.0 - val1);
          dN_dxi[1] = -val1;
          dN_dxi[2] =  0.5*(1.0 + val1);

      break;

      case 3:

          fact1 = 9.0/16.0;
          fact2 = 27.0/16.0;
          val1 = xi*xi;

          N[0] = -fact1 * (1 - xi)   * (1/9 - val1);
          N[1] =  fact2 * (1 - val1) * (1/3 - xi);
          N[2] =  fact2 * (1 - val1) * (1/3 + xi);
          N[3] = -fact1 * (1 + xi)   * ( 1/9 - val1);

          val2 = 3.0*val1;

          dN_dxi[0] = -fact1*(-1/9 - 2.0*xi   +  val2);
          dN_dxi[1] =  fact2*(-1   - 2.0/3*xi +  val2);
          dN_dxi[2] =  fact2*(1    - 2.0/3*xi -  val2);
          dN_dxi[3] = -fact1*(1/9  - 2.0*xi   -  val2);

      break;

      case 4:

          fact1 = 2.0/3.0;
          fact2 = 8.0/3.0;
          val1 = xi*xi;
          val2 = val1*xi;
          val3 = val2*xi;

          N[0] =  fact1 * (0.25*xi  - 0.25*val1  -  val2     + val3);
          N[1] = -fact2 * (0.5 *xi  - val1       -  0.5*val2 + val3);
          N[2] =    4.0 * (0.25     - 1.25*val1  -  0        + val3);
          N[3] =  fact2 * (0.5*xi   + val1       -  0.5*val2 - val3);
          N[4] = -fact1 * (0.25*xi  + 0.25*val1  -  val2     - val3);

          val4 = 4.0*val2;

          dN_dxi[0] =  fact1 * (0.25 - 0.5*xi  - 3.0*val1  + val4);
          dN_dxi[1] = -fact2 * (0.5  - 2.0*xi  - 1.5*val1  + val4);
          dN_dxi[2] =    4.0 * (  0  - 2.5*xi  -   0.0     + val4);
          dN_dxi[3] =  fact2 * (0.5  + 2.0*xi  - 1.5*val1  - val4);
          dN_dxi[4] = -fact1 * (0.25 + 0.5*xi  - 3.0*val1  - val4);

      break;

      default:

          printf("no basis functions defined for this degree = %5d \n", p);

      break;
  }

  return;
}



void Lagrange_BasisFuns1D(int p, double xi, double* N)
{
  double  fact1, fact2, val1, val2, val3, val4;

  switch(p)
  {
      case 0:

          N[0] = 1.0;

      break;

      case 1:

          N[0] = 0.5*(1.0 - xi);
          N[1] = 0.5*(1.0 + xi);

      break;

      case 2:

          val1 = xi*xi;

          N[0] = -0.5 * (xi - val1);
          N[1] = 1.0 - val1;
          N[2] = 0.5 *(xi + val1);

      break;

      case 3:

          fact1 = 9.0/16.0;
          fact2 = 27.0/16.0;
          val1 = xi*xi;

          N[0] = -fact1 * (1 - xi)   * (1/9 - val1);
          N[1] =  fact2 * (1 - val1) * (1/3 - xi);
          N[2] =  fact2 * (1 - val1) * (1/3 + xi);
          N[3] = -fact1 * (1 + xi)   * ( 1/9 - val1);
      break;

      case 4:

          fact1 = 2.0/3.0;
          fact2 = 8.0/3.0;
          val1 = xi*xi;
          val2 = val1*xi;
          val3 = val2*xi;

          N[0] =  fact1 * (0.25*xi  - 0.25*val1  -  val2     + val3);
          N[1] = -fact2 * (0.5 *xi  - val1       -  0.5*val2 + val3);
          N[2] =    4.0 * (0.25     - 1.25*val1  -  0        + val3);
          N[3] =  fact2 * (0.5*xi   + val1       -  0.5*val2 - val3);
          N[4] = -fact1 * (0.25*xi  + 0.25*val1  -  val2     - val3);

      break;

      default:

          printf("no basis functions defined for this degree = %5d \n", p);

      break;
  }

   return;
}



void LagrangeBasisFunsTria(int p, double xi1, double xi2, double* N)
{
  double  xi3 = 1.0 - xi1 - xi2;

  switch(p)
  {
      case 0:

          N[0] = 1.0;

      break;

      case 1:

          N[0] = xi3;
          N[1] = xi1;
          N[2] = xi2;

      break;

      case 2:

          N[0] = xi3*(2.0*xi3 - 1.0);
          N[1] = xi1*(2.0*xi1 - 1.0);
          N[2] = xi2*(2.0*xi2 - 1.0);
          N[3] = 4.0*xi1*xi3;
          N[4] = 4.0*xi1*xi2;
          N[5] = 4.0*xi2*xi3;

      break;

      default:

          printf("no basis functions defined for this degree = %5d \n", p);

      break;
  }
  return;
}


void LagrangeBasisFunsTria(int p, double xi1, double xi2, double* N, double* dN_dxi1, double* dN_dxi2)
{
  double  xi3 = 1.0 - xi1 - xi2;

  switch(p)
  {
      case 0:

          N[0] = 1.0;

          dN_dxi1[0] = 0.0;
          dN_dxi2[0] = 0.0;

      break;

      case 1:

          N[0] = xi3;
          N[1] = xi1;
          N[2] = xi2;

          dN_dxi1[0] = -1.0;
          dN_dxi1[1] =  1.0;
          dN_dxi1[2] =  0.0;

          dN_dxi2[0] = -1.0;
          dN_dxi2[1] =  0.0;
          dN_dxi2[2] =  1.0;

      break;

      case 2:

          N[0] = xi3*(2.0*xi3 - 1.0);
          N[1] = xi1*(2.0*xi1 - 1.0);
          N[2] = xi2*(2.0*xi2 - 1.0);
          N[3] = 4.0*xi1*xi3;
          N[4] = 4.0*xi1*xi2;
          N[5] = 4.0*xi2*xi3;

          dN_dxi1[0] = -4.0*xi3 + 1.0;
          dN_dxi1[1] =  4.0*xi1 - 1.0;
          dN_dxi1[2] =  0.0;
          dN_dxi1[3] =  4.0*(xi3 - xi1);
          dN_dxi1[4] =  4.0*xi2;
          dN_dxi1[5] = -4.0*xi2;

          dN_dxi2[0] = -4.0*xi3 + 1.0;
          dN_dxi2[1] =  0.0;
          dN_dxi2[2] =  4.0*xi2 - 1.0;
          dN_dxi2[3] = -4.0*xi1;
          dN_dxi2[4] =  4.0*xi1;
          dN_dxi2[5] =  4.0*(xi3 - xi2);

      break;

      default:

          printf("no basis functions defined for this degree = %5d \n", p);

      break;
  }
  return;
}




void BernsteinBasisFunsTria(int p, double xi1, double xi2, double* N)
{
  double  xi3 = 1.0 - xi1 - xi2;

  switch(p)
  {
      case 0:

          N[0] = 1.0;

      break;

      case 1:

          N[0] = xi3;
          N[1] = xi1;
          N[2] = xi2;

      break;

      case 2:

          N[0] = xi3*xi3;
          N[1] = xi1*xi1;
          N[2] = xi2*xi2;
          N[3] = 2.0*xi1*xi3;
          N[4] = 2.0*xi1*xi2;
          N[5] = 2.0*xi2*xi3;

      break;

      default:

          printf("no basis functions defined for this degree = %5d \n", p);

      break;
  }
  return;
}


void LagrangeBasisFunsQuad(int p, double xi1, double xi2, double* N)
{
  double  a1, a2, a3, b1, b2, b3;

  switch(p)
  {
      case 0:

          N[0] = 1.0;

      break;

      case 1:

          b1 = 0.5*(1.0 - xi2);
          b2 = 0.5*(1.0 + xi2);

          a1 = 0.5*(1.0 - xi1);
          a2 = 0.5*(1.0 + xi1);

          N[0] = b1*a1;
          N[1] = b1*a2;
          N[2] = b2*a2;
          N[3] = b2*a1;

      break;

      case 2:

          b1 = 0.5*(xi2*xi2 - xi2);
          b2 = 1.0 - xi2*xi2;
          b3 = 0.5*(xi2*xi2 + xi2);

          a1 = 0.5*(xi1*xi1 - xi1);
          a2 = 1.0 - xi1*xi1;
          a3 = 0.5*(xi1*xi1 + xi1);

          N[0] = b1*a1;
          N[1] = b1*a3;
          N[2] = b3*a3;
          N[3] = b3*a1;
          N[4] = b1*a2;
          N[5] = b2*a3;
          N[6] = b3*a2;
          N[7] = b2*a1;
          N[8] = b2*a2;

      break;

      default:

          printf("no basis functions defined for this degree = %5d \n", p);

      break;
  }

  return;
}



void LagrangeBasisFunsQuad(int p, double xi1, double xi2, double* N, double* dN_dxi1, double* dN_dxi2)
{
  double  a1, a2, a3, b1, b2, b3;
  double  da1, da2, da3, db1, db2, db3;

  switch(p)
  {
      case 0:

          N[0] = 1.0;

          dN_dxi1[0] = 0.0;
          dN_dxi2[0] = 0.0;

      break;

      case 1:

          b1 = 0.5*(1.0 - xi2);
          b2 = 0.5*(1.0 + xi2);

          db1 = -0.5;
          db2 =  0.5;

          a1 = 0.5*(1.0 - xi1);
          a2 = 0.5*(1.0 + xi1);

          da1 = -0.5;
          da2 =  0.5;

          N[0] = b1*a1;
          N[1] = b1*a2;
          N[2] = b2*a2;
          N[3] = b2*a1;

          dN_dxi1[0] = b1*da1;
          dN_dxi1[1] = b1*da2;
          dN_dxi1[2] = b2*da2;
          dN_dxi1[3] = b2*da1;

          dN_dxi2[0] = db1*a1;
          dN_dxi2[1] = db1*a2;
          dN_dxi2[2] = db2*a2;
          dN_dxi2[3] = db2*a1;

      break;

      case 2:

          b1 = 0.5*(xi2*xi2 - xi2);
          b2 = 1.0 - xi2*xi2;
          b3 = 0.5*(xi2*xi2 + xi2);

          db1 = 0.5*(2.0*xi2 - 1.0);
          db2 = -2.0*xi2;
          db3 = 0.5*(2.0*xi2 + 1.0);

          a1 = 0.5*(xi1*xi1 - xi1);
          a2 = 1.0 - xi1*xi1;
          a3 = 0.5*(xi1*xi1 + xi1);

          da1 = 0.5*(2.0*xi1 - 1.0);
          da2 = -2.0*xi1;
          da3 = 0.5*(2.0*xi1 + 1.0);

/*
          N[0] = b1*a1;
          N[1] = b1*a3;
          N[2] = b3*a3;
          N[3] = b3*a1;
          N[4] = b1*a2;
          N[5] = b2*a3;
          N[6] = b3*a2;
          N[7] = b2*a1;

          dN_dxi1[0] = b1*da1;
          dN_dxi1[1] = b1*da3;
          dN_dxi1[2] = b3*da3;
          dN_dxi1[3] = b3*da1;
          dN_dxi1[4] = b1*da2;
          dN_dxi1[5] = b2*da3;
          dN_dxi1[6] = b3*da2;
          dN_dxi1[7] = b2*da1;

          dN_dxi2[0] = db1*a1;
          dN_dxi2[1] = db1*a3;
          dN_dxi2[2] = db3*a3;
          dN_dxi2[3] = db3*a1;
          dN_dxi2[4] = db1*a2;
          dN_dxi2[5] = db2*a3;
          dN_dxi2[6] = db3*a2;
          dN_dxi2[7] = db2*a1;
*/
//
          N[0] = b1*a1;
          N[1] = b1*a3;
          N[2] = b3*a3;
          N[3] = b3*a1;
          N[4] = b1*a2;
          N[5] = b2*a3;
          N[6] = b3*a2;
          N[7] = b2*a1;
          N[8] = b2*a2;

          dN_dxi1[0] = b1*da1;
          dN_dxi1[1] = b1*da3;
          dN_dxi1[2] = b3*da3;
          dN_dxi1[3] = b3*da1;
          dN_dxi1[4] = b1*da2;
          dN_dxi1[5] = b2*da3;
          dN_dxi1[6] = b3*da2;
          dN_dxi1[7] = b2*da1;
          dN_dxi1[8] = b2*da2;

          dN_dxi2[0] = db1*a1;
          dN_dxi2[1] = db1*a3;
          dN_dxi2[2] = db3*a3;
          dN_dxi2[3] = db3*a1;
          dN_dxi2[4] = db1*a2;
          dN_dxi2[5] = db2*a3;
          dN_dxi2[6] = db3*a2;
          dN_dxi2[7] = db2*a1;
          dN_dxi2[8] = db2*a2;
//
      break;

      default:

          printf("no basis functions defined for this degree = %5d \n", p);

      break;
  }

  return;
}



void LagrangeBasisFunsTetra(int p, double xi1, double xi2, double xi3, double* N)
{
  double  xi4 = 1.0 - xi1 - xi2 - xi3;

  switch(p)
  {
      case 0:

          N[0] = 1.0;

      break;

      case 1:

          N[0] = xi4;
          N[1] = xi1;
          N[2] = xi2;
          N[3] = xi3;

      break;

      case 2:

          //xi4 = 1.0 - xi1 - xi2 - xi3;

          N[0] = xi4*(2.0*xi4 - 1.0);
          N[1] = xi1*(2.0*xi1 - 1.0);
          N[2] = xi2*(2.0*xi2 - 1.0);
          N[3] = xi3*(2.0*xi3 - 1.0);
          N[4] = 4.0*xi1*xi4;
          N[5] = 4.0*xi1*xi2;
          N[6] = 4.0*xi2*xi4;
          N[7] = 4.0*xi4*xi3;
          N[8] = 4.0*xi1*xi3;
          N[9] = 4.0*xi2*xi3;

      break;

      default:

          printf("no basis functions defined for this degree = %5d \n", p);

      break;
  }

}


void LagrangeBasisFunsTetra(int p, double xi1, double xi2, double xi3, double* N, double* dN_dxi1, double* dN_dxi2, double* dN_dxi3)
{
  double  xi4 = 1.0 - xi1 - xi2 - xi3;

  switch(p)
  {
      case 0:

          N[0] = 1.0;

          dN_dxi1[0] = 0.0;
          dN_dxi2[0] = 0.0;
          dN_dxi3[0] = 0.0;

      break;

      case 1:

          N[0] = xi4;
          N[1] = xi1;
          N[2] = xi2;
          N[3] = xi3;

          dN_dxi1[0] = -1.0;
          dN_dxi1[1] =  1.0;
          dN_dxi1[2] =  0.0;
          dN_dxi1[3] =  0.0;

          dN_dxi2[0] = -1.0;
          dN_dxi2[1] =  0.0;
          dN_dxi2[2] =  1.0;
          dN_dxi2[3] =  0.0;

          dN_dxi3[0] = -1.0;
          dN_dxi3[1] =  0.0;
          dN_dxi3[2] =  0.0;
          dN_dxi3[3] =  1.0;

      break;

      case 2:

          //xi4 = 1.0 - xi1 - xi2 - xi3;

          N[0] = xi4*(2.0*xi4 - 1.0);
          N[1] = xi1*(2.0*xi1 - 1.0);
          N[2] = xi2*(2.0*xi2 - 1.0);
          N[3] = xi3*(2.0*xi3 - 1.0);
          N[4] = 4.0*xi1*xi4;
          N[5] = 4.0*xi1*xi2;
          N[6] = 4.0*xi2*xi4;
          N[7] = 4.0*xi4*xi3;
          N[8] = 4.0*xi1*xi3;
          N[9] = 4.0*xi2*xi3;

          dN_dxi1[0] = -4.0*xi4 + 1.0;
          dN_dxi1[1] =  4.0*xi1 - 1.0;
          dN_dxi1[2] =  0.0;
          dN_dxi1[3] =  0.0;
          dN_dxi1[4] =  4.0*(xi4-xi1);
          dN_dxi1[5] =  4.0*xi2;
          dN_dxi1[6] = -4.0*xi2;
          dN_dxi1[7] = -4.0*xi3;
          dN_dxi1[8] =  4.0*xi3;
          dN_dxi1[9] =  0.0;

          dN_dxi2[0] = -4.0*xi4 + 1.0;
          dN_dxi2[1] =  0.0;
          dN_dxi2[2] =  4.0*xi2 - 1.0;
          dN_dxi2[3] =  0.0;
          dN_dxi2[4] = -4.0*xi1;
          dN_dxi2[5] =  4.0*xi1;
          dN_dxi2[6] =  4.0*(xi4-xi2);
          dN_dxi2[7] = -4.0*xi3;
          dN_dxi2[8] =  0.0;
          dN_dxi2[9] =  4.0*xi3;

          dN_dxi3[0] = -4.0*xi4 + 1.0;
          dN_dxi3[1] =  0.0;
          dN_dxi3[2] =  0.0;
          dN_dxi3[3] =  4.0*xi3 - 1.0;
          dN_dxi3[4] = -4.0*xi1;
          dN_dxi3[5] =  0.0;
          dN_dxi3[6] = -4.0*xi2;
          dN_dxi3[7] =  4.0*(xi4-xi3);
          dN_dxi3[8] =  4.0*xi1;
          dN_dxi3[9] =  4.0*xi2;

      break;

      default:

          printf("no basis functions defined for this degree = %5d \n", p);

      break;
  }

}



void LagrangeBasisFunsHexa(int p, double xi1, double xi2, double xi3, double* N)
{
  double  v11, v12, v21, v22, v31, v32;

  switch(p)
  {
      case 0:

          N[0] = 1.0;

      break;

      case 1:

          v11 = 1.0 - xi1;
          v12 = 1.0 + xi1;
          v21 = 1.0 - xi2;
          v22 = 1.0 + xi2;
          v31 = 1.0 - xi3;
          v32 = 1.0 + xi3;

          N[0] = 0.125*v11*v21*v31;
          N[1] = 0.125*v12*v21*v31;
          N[2] = 0.125*v12*v22*v31;
          N[3] = 0.125*v11*v22*v31;
          N[4] = 0.125*v11*v21*v32;
          N[5] = 0.125*v12*v21*v32;
          N[6] = 0.125*v12*v22*v32;
          N[7] = 0.125*v11*v22*v32;

      break;

      default:

          printf("no basis functions defined for this degree = %5d in 'LagrangeBasisFunsPenta' \n", p);

      break;
  }

  return;
}




void LagrangeBasisFunsHexa(int p, double xi1, double xi2, double xi3, double* N, double* dN_dxi1, double* dN_dxi2, double* dN_dxi3)
{
  double  v11, v12, v21, v22, v31, v32;

  switch(p)
  {
      case 0:

          N[0] = 1.0;

          dN_dxi1[0] = 0.0;
          dN_dxi2[0] = 0.0;
          dN_dxi3[0] = 0.0;

      break;

      case 1:

          v11 = 1.0 - xi1;
          v12 = 1.0 + xi1;
          v21 = 1.0 - xi2;
          v22 = 1.0 + xi2;
          v31 = 1.0 - xi3;
          v32 = 1.0 + xi3;

          N[0] = 0.125*v11*v21*v31;
          N[1] = 0.125*v12*v21*v31;
          N[2] = 0.125*v12*v22*v31;
          N[3] = 0.125*v11*v22*v31;
          N[4] = 0.125*v11*v21*v32;
          N[5] = 0.125*v12*v21*v32;
          N[6] = 0.125*v12*v22*v32;
          N[7] = 0.125*v11*v22*v32;

          dN_dxi1[0] = -0.125*v21*v31;
          dN_dxi1[1] =  0.125*v21*v31;
          dN_dxi1[2] =  0.125*v22*v31;
          dN_dxi1[3] = -0.125*v22*v31;
          dN_dxi1[4] = -0.125*v21*v32;
          dN_dxi1[5] =  0.125*v21*v32;
          dN_dxi1[6] =  0.125*v22*v32;
          dN_dxi1[7] = -0.125*v22*v32;

          dN_dxi2[0] = -0.125*v11*v31;
          dN_dxi2[1] = -0.125*v12*v31;
          dN_dxi2[2] =  0.125*v11*v31;
          dN_dxi2[3] =  0.125*v12*v31;
          dN_dxi2[4] = -0.125*v11*v32;
          dN_dxi2[5] = -0.125*v12*v32;
          dN_dxi2[6] =  0.125*v11*v32;
          dN_dxi2[7] =  0.125*v12*v32;

          dN_dxi3[0] = -0.125*v11*v21;
          dN_dxi3[1] = -0.125*v12*v21;
          dN_dxi3[2] = -0.125*v11*v22;
          dN_dxi3[3] = -0.125*v12*v22;
          dN_dxi3[4] =  0.125*v11*v21;
          dN_dxi3[5] =  0.125*v12*v21;
          dN_dxi3[6] =  0.125*v11*v22;
          dN_dxi3[7] =  0.125*v12*v22;

      break;

      default:

          printf("no basis functions defined for this degree = %5d in 'LagrangeBasisFunsPenta' \n", p);

      break;
  }

  return;
}


void LagrangeBasisFunsPrism(int p, double xi1, double xi2, double xi4, double* N, double* dN_dxi1, double* dN_dxi2, double* dN_dxi4)
{
  // xi1, xi2 and xi3 are the parametric coordinates for the triangle, and
  // xi4 is the parametric coordinate in the direction normal to the triangle
  double  xi3 = 1.0- xi1 - xi2;

  switch(p)
  {
      case 1:

          N[0] = xi3*(0.5*(1.0-xi4));
          N[1] = xi1*(0.5*(1.0-xi4));
          N[2] = xi2*(0.5*(1.0-xi4));
          N[3] = xi3*(0.5*(1.0+xi4));
          N[4] = xi1*(0.5*(1.0+xi4));
          N[5] = xi2*(0.5*(1.0+xi4));


          dN_dxi1[0]   = (-1.0)*(0.5*(1.0-xi4));
          dN_dxi1[1]   = ( 1.0)*(0.5*(1.0-xi4));
          dN_dxi1[2]   = ( 0.0)*(0.5*(1.0-xi4));
          dN_dxi1[3]   = (-1.0)*(0.5*(1.0+xi4));
          dN_dxi1[4]   = ( 1.0)*(0.5*(1.0+xi4));
          dN_dxi1[5]   = ( 0.0)*(0.5*(1.0+xi4));

          dN_dxi2[0]  = (-1.0)*(0.5*(1.0-xi4));
          dN_dxi2[1]  = ( 0.0)*(0.5*(1.0-xi4));
          dN_dxi2[2]  = ( 1.0)*(0.5*(1.0-xi4));
          dN_dxi2[3]  = (-1.0)*(0.5*(1.0+xi4));
          dN_dxi2[4]  = ( 0.0)*(0.5*(1.0+xi4));
          dN_dxi2[5]  = ( 1.0)*(0.5*(1.0+xi4));

          dN_dxi4[0]  = xi3*(-0.5);
          dN_dxi4[1]  = xi1*(-0.5);
          dN_dxi4[2]  = xi2*(-0.5);
          dN_dxi4[3]  = xi3*( 0.5);
          dN_dxi4[4]  = xi1*( 0.5);
          dN_dxi4[5]  = xi2*( 0.5);

      break;

      default:

          printf("no basis functions defined for this degree = %5d in 'LagrangeBasisFunsPenta' \n", p);

      break;
  }
  return;
}




void LagrangeBasisFunsPyramid(int p, double xi1, double xi2, double xi3, double* N, double* dN_dxi1, double* dN_dxi2, double* dN_dxi3)
{
 double  v11, v12, v21, v22, v31, v32;

  switch(p)
  {
      case 1:

          v11 = 1.0 - xi1;
          v12 = 1.0 + xi1;
          v21 = 1.0 - xi2;
          v22 = 1.0 + xi2;
          v31 = 1.0 - xi3;
          v32 = 1.0 + xi3;

          N[0] = 0.125*v11*v21*v31;
          N[1] = 0.125*v12*v21*v31;
          N[2] = 0.125*v12*v22*v31;
          N[3] = 0.125*v11*v22*v31;
          N[4] = 0.5*v32;

          dN_dxi1[0] = -0.125*v21*v31;
          dN_dxi1[1] =  0.125*v21*v31;
          dN_dxi1[2] =  0.125*v22*v31;
          dN_dxi1[3] = -0.125*v22*v31;
          dN_dxi1[4] =  0.0;

          dN_dxi2[0] = -0.125*v11*v31;
          dN_dxi2[1] = -0.125*v12*v31;
          dN_dxi2[2] =  0.125*v11*v31;
          dN_dxi2[3] =  0.125*v12*v31;
          dN_dxi2[4] =  0.0;

          dN_dxi3[0] = -0.125*v11*v21;
          dN_dxi3[1] = -0.125*v12*v21;
          dN_dxi3[2] = -0.125*v11*v22;
          dN_dxi3[3] = -0.125*v12*v22;
          dN_dxi3[4] =  0.5;

      break;

      default:

          printf("no basis functions defined for this degree = %5d in 'LagrangeBasisFunsPyramid' \n", p);

      break;
  }
  return;
}





void LagrangeBasisFunsLine1D(int p, double uu, double *xx, double *N, double *dN_dx, double& Jac)
{
  int  ii, jj, count, nlbf = p+1;

  double dx_du ;
  vector<double>  dN1(nlbf);

  Lagrange_BasisFuns1D(p, uu, N, &dN1[0]);

  if(p == 0)
  {
    Jac = 1.0;
  }
  else
  {
    Jac = 0.0;
    for(ii=0; ii<nlbf; ii++)
      Jac +=  (xx[ii] * dN1[ii]);
  }
  
  dx_du = 1.0/Jac ;

  // Compute derivatives of basis functions w.r.t physical coordinates
  for(ii=0; ii<nlbf; ii++)
    dN_dx[ii] = dN1[ii] * dx_du ;
  
  return;
}




void LagrangeBasisFunsLine3D(int p, double uu, double *xx, double* yy, double* zz, double *N, double *dN_dx, double& Jac)
{
   int  ii, jj, count, nlbf = p+1;

   double  du_dx, dx, dy, dz ;
   vector<double>  dN1(nlbf);

   Lagrange_BasisFuns1D(p, uu, N, &dN1[0]);

   if(p == 0)
   {
     Jac = 1.0;
   }
   else
   {
     dx = dy = dz = 0.0;
     for(ii=0; ii<nlbf; ii++)
     {
       dx +=  (xx[ii] * dN1[ii]);
       dy +=  (yy[ii] * dN1[ii]);
       dz +=  (zz[ii] * dN1[ii]);
     }
     Jac = sqrt(dx*dx+dy*dy+dz*dz);
   }

   du_dx = 1.0/Jac ;

  // Compute derivatives of basis functions w.r.t physical coordinates
  for(ii=0; ii<nlbf; ii++)
    dN_dx[ii] = dN1[ii] * du_dx ;

  return;
}





void LagrangeBasisFunsEdge2D(int p, double* param, double *xNode, double* yNode, double *N, double *normal, double& Jac)
{
   int  ii, jj, count, nlbf = p+1;

   double  du_dx, dx, dy ;
   vector<double>  dN1(nlbf);

   Lagrange_BasisFuns1D(p, param[0], N, &dN1[0]);

   if(p == 0)
   {
     Jac = 1.0;
   }
   else
   {
     dx = dy = 0.0;
     for(ii=0; ii<nlbf; ii++)
     {
       dx +=  (xNode[ii] * dN1[ii]);
       dy +=  (yNode[ii] * dN1[ii]);
     }
     Jac = sqrt(dx*dx+dy*dy);
   }

   // Compute the normal
   normal[0] =  dy/Jac;
   normal[1] = -dx/Jac;

   return;
}




void LagrangeBasisFunsEdge3D(int p, double* param, double *xNode, double* yNode, double *zNode, double *N, double *normal, double& Jac)
{
  return;
}




void LagrangeBasisFunsFaceTria(int degree, double* param, double *xNode, double* yNode, double *zNode, double *N, double *normal, double& Jac)
{
   int  ii, jj, count, nlbf;
   if(degree == 1)
   {
     nlbf = 3;
   }
   else if(degree == 2)
   {
     nlbf = 6;
   }

   Vector3d  du, dv;
   vector<double>  dN1(nlbf), dN2(nlbf);

   LagrangeBasisFunsTria(degree, param[0], param[1], N, &dN1[0], &dN2[0]);

   du.setZero();
   dv.setZero();
   for(ii=0; ii<nlbf; ii++)
   {
     du(0) += (xNode[ii] * dN1[ii]);
     du(1) += (yNode[ii] * dN1[ii]);
     du(2) += (zNode[ii] * dN1[ii]);

     dv(0) += (xNode[ii] * dN2[ii]);
     dv(1) += (yNode[ii] * dN2[ii]);
     dv(2) += (zNode[ii] * dN2[ii]);
   }

   Vector3d normal1 = du.cross(dv);

   Jac = normal1.norm();
   normal1 /= Jac;

   // Compute the normal
   normal[0] = normal1[0];
   normal[1] = normal1[1];
   normal[2] = normal1[2];

   return;
}


void LagrangeBasisFunsFaceQuad(int p, double* param, double *xNode, double* yNode, double *zNode, double *N, double *normal, double& Jac)
{
  return;
}





