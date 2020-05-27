
#include "LagrangeElem3DNavierStokesTetra4Node.h"
#include "SolutionData.h"
#include "UtilitiesQuadrature.h"
#include "stabilisationRoutines.h"
#include "KimMoinFlow.h"
#include "BasisFunctionsLagrange.h"


using namespace std;


LagrangeElem3DNavierStokesTetra4Node::LagrangeElem3DNavierStokesTetra4Node()
{
  degree = 1;
  npElem = 4;
  nlbf   = 4;
  ndof   = 4;
  nsize  = nlbf*ndof;

  nGP = 1;
}


LagrangeElem3DNavierStokesTetra4Node::~LagrangeElem3DNavierStokesTetra4Node()
{
}



void LagrangeElem3DNavierStokesTetra4Node::prepareElemData(vector<vector<double> >& node_coords)
{
    // compute Volume and basis function derivatives

    double  dvol, Jac, param[2];

    int   ii, gp;

    double xNode[npElem], yNode[npElem], zNode[npElem], xx, yy, zz;
    for(ii=0;ii<npElem;ii++)
    {
      //xNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][0];
      //yNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][1];

      xNode[ii] = node_coords[nodeNums[ii]][0];
      yNode[ii] = node_coords[nodeNums[ii]][1];
      zNode[ii] = node_coords[nodeNums[ii]][2];
    }

    // Gauss point coordinates and weights
    vector<double>  gpts1(nGP), gpts2(nGP), gpts3(nGP), gwts(nGP);
    getGaussPointsTetra(nGP, gpts1, gpts2, gpts3, gwts);

    Nv.resize(nGP);
    dNvdx.resize(nGP);
    dNvdy.resize(nGP);
    dNvdz.resize(nGP);

    for(int ii=0; ii<nGP; ii++)
    {
      Nv[ii].resize(nlbf);
      dNvdx[ii].resize(nlbf);
      dNvdy[ii].resize(nlbf);
      dNvdz[ii].resize(nlbf);

      Nv[ii].setZero();
      dNvdx[ii].setZero();
      dNvdy[ii].setZero();
      dNvdz[ii].setZero();
    }

    //printVector(nodeNums);

    elemVolGP.resize(nGP);

    elemVol=0.0;
    for(gp=0; gp<nGP; gp++)
    {
      param[0] = gpts1[gp];
      param[1] = gpts2[gp];
      param[2] = gpts3[gp];

      computeBasisFunctions3D(false, ELEM_SHAPE_TETRA, degree, param, xNode, yNode, zNode, &Nv[gp][0], &dNvdx[gp][0], &dNvdy[gp][0], &dNvdz[gp][0], Jac);

      //cout << " Jac = " << Jac << endl;
      /*
      if(Jac < 0.0)
      {
        cout << " Negative Jacobian in 'LagrangeElem2DNavierStokesQuad4Node::prepareElemData' " << endl;
        cout << " Jac = " << Jac << endl;
        exit[1];
      }
      */

      dvol = gwts[gp]*Jac;
      elemVolGP[gp] = dvol;

      elemVol += dvol;
    }//gp

    //double  volume = 0.5*( (xNode[0]-xNode[2])*(yNode[1]-yNode[3]) + (xNode[1]-xNode[3])*(yNode[2]-yNode[0]) );
    //cout << "volume = " << elemVol <<  '\t' <<  volume << endl;
    charlen = sqrt(4.0*abs(elemVol)/PI);

    globalDOFnums.resize(nsize);

    for(ii=0; ii<npElem; ii++)
    {
      globalDOFnums[ii*ndof]   = nodeNums[ii]*ndof;
      globalDOFnums[ii*ndof+1] = nodeNums[ii]*ndof + 1;
      globalDOFnums[ii*ndof+2] = nodeNums[ii]*ndof + 2;
      globalDOFnums[ii*ndof+3] = nodeNums[ii]*ndof + 3;
    }

    return;
}


//
int LagrangeElem3DNavierStokesTetra4Node::calcStiffnessAndResidual(vector<vector<double> >& node_coords, double* elemData, MatrixXdRM& Klocal, VectorXd& Flocal, double timeCur)
{
/*
  // Fully-implicit formulation

  double  dvol, b1, b2, b3, b4, b5, b6, b7, b8, Da, Db, c1, CI=4.0;
  double  Jac, totvol, pres, fact, fact2;
  double  param[3], bforce[3], tau[3], beta[6];

  VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf), dN_dz(nlbf);
  VectorXd  res(4), res2(3), dp(3), Du(3), vel(3), force(3), gradTvel(3), rStab(3);
  VectorXd  velPrev(3), velExt(3), velTemp(3), velDot(3);
  MatrixXd  grad(3,3), gradN(3,3), stress(3,3),  Dj(3,4), forVolume(4,4);
  Dj.setZero();

  int   index, ii, jj, gp, TI, TIp1, TIp2, TIp3, TJ, TJp1, TJp2, TJp3;

  elmDat = &(SolnData->ElemProp[elmType].data[0]);

  double  rho = elmDat[3];
  double  mu  = elmDat[4];

  bforce[0]  = 0.0;
  bforce[1]  = 0.0;
  bforce[2]  = 0.0;

  //bforce[0]  = elmDat[6]*timeFunction[0].prop ;
  //bforce[1]  = elmDat[7]*timeFunction[0].prop ;
  //bforce[2]  = elmDat[7]*timeFunction[0].prop ;

  double  dt = mpapTime.dt;
  double  af = SolnData->td(2);
  double  am = SolnData->td(1);
  double  acceFact = am*SolnData->td(9);
  double  muTaf = mu*af;

  double  tCur = mpapTime.cur - (1.0-af)*dt;

  double xNode[4], yNode[4], zNode[4], xx, yy, zz;
  for(ii=0;ii<npElem;ii++)
  {
    xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
    yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
  }

  forVolume(0,0) = 1.0;      forVolume(0,1) = 1.0;      forVolume(0,2) = 1.0;      forVolume(0,3) = 1.0;
  forVolume(1,0) = xNode[0]; forVolume(1,1) = xNode[1]; forVolume(1,2) = xNode[2]; forVolume(1,3) = xNode[3];
  forVolume(2,0) = yNode[0]; forVolume(2,1) = yNode[1]; forVolume(2,2) = yNode[2]; forVolume(2,3) = yNode[3];
  forVolume(3,0) = zNode[0]; forVolume(3,1) = zNode[1]; forVolume(3,2) = zNode[2]; forVolume(3,3) = zNode[3];

  totvol = forVolume.determinant()/6.0;
  //if(elenum < 20)
    //cout << " totvol = " << totvol << endl;

  double  h  = pow(6.0*totvol/PI, 1.0/3.0);

  double stabParam = h*h/(4.0*mu);
  tau[0] = elmDat[8]*stabParam;  // SUPG
  tau[1] = elmDat[9]*stabParam;  // PSPG
  tau[2] = elmDat[10]*stabParam; // LSIC

  //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", xNode[0], xNode[1], xNode[2], xNode[3]);
  //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n\n\n ", yNode[0], yNode[1], yNode[2], yNode[3]);

  //nGP = (int) elmDat[0] ;
  nGP = 1;

  vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;

  getGaussPointsTetra(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);

  if(Klocal.rows() != nsize)
  {
    Klocal.resize(nsize, nsize);
    Flocal.resize(nsize);
  }
  Klocal.setZero();
  Flocal.setZero();

  //cout << " AAAAAAAAAA " << endl;

  for(gp=0;gp<nGP;gp++)
  {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];
        param[2] = gausspoints3[gp];

        GeomData->computeBasisFunctions3D(0, 1, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        dvol = Jac * gaussweights[gp];

          vel(0) = computeValueCur(0, N);
          vel(1) = computeValueCur(1, N);
          vel(2) = computeValueCur(2, N);

          velPrev(0) = computeValuePrev(0, N);
          velPrev(1) = computeValuePrev(1, N);
          velPrev(2) = computeValuePrev(2, N);

          grad(0,0) = computeValueCur(0, dN_dx);
          grad(0,1) = computeValueCur(0, dN_dy);
          grad(0,2) = computeValueCur(0, dN_dz);

          grad(1,0) = computeValueCur(1, dN_dx);
          grad(1,1) = computeValueCur(1, dN_dy);
          grad(1,2) = computeValueCur(1, dN_dz);

          grad(2,0) = computeValueCur(2, dN_dx);
          grad(2,1) = computeValueCur(2, dN_dy);
          grad(2,2) = computeValueCur(2, dN_dz);

          Du.setZero();

          pres   = computeValue(3, N);
          dp(0)  = computeValue(3, dN_dx);
          dp(1)  = computeValue(3, dN_dy);
          dp(2)  = computeValue(3, dN_dz);

          velDot(0) = computeValueDotCur(0, N);
          velDot(1) = computeValueDotCur(1, N);
          velDot(2) = computeValueDotCur(2, N);

          // this is pseudo-stress
          stress = mu*grad;
          stress(0,0) -= pres;
          stress(1,1) -= pres;
          stress(2,2) -= pres;

          //cout << xx << '\t' << yy << endl;

          force.setZero();
          //force(0) = analy.computeForce(0, xx, yy, tCur);
          //force(1) = analy.computeForce(1, xx, yy, tCur);
          //cout << force(0) << '\t' << force(1) << endl;

          gradTvel = grad*vel;

          res2(0) = rho*(velDot(0) + gradTvel(0) - force(0)) ;
          res2(1) = rho*(velDot(1) + gradTvel(1) - force(1)) ;
          res2(2) = rho*(velDot(2) + gradTvel(2) - force(2)) ;

          rStab(0) = res2(0) - mu*Du(0) + dp(0) ;
          rStab(1) = res2(1) - mu*Du(1) + dp(1) ;
          rStab(2) = res2(2) - mu*Du(2) + dp(2) ;

          velTemp = velPrev;
          //velTemp = vel;

          //evaluateStabParams_algo1(&velTemp(0), h, rho, mu, dt,  beta, tau);

          //evaluateStabParams_algo2(&velTemp(0), h, rho, mu, dt,  beta, tau);

          //evaluateStabParams_algo3(velTemp, matG, dt, rho, mu, CI, tau);

          //tau[0] *= elmDat[8];  // SUPG
          //tau[1] *= elmDat[9];  // PSPG
          //tau[2] *= elmDat[10]; // LSIC

          for(ii=0;ii<nlbf;ii++)
          {
            TI   = 4*ii;
            TIp1 = TI+1;
            TIp2 = TI+2;
            TIp3 = TI+3;

            b1 = dN_dx[ii]*dvol;
            b2 = dN_dy[ii]*dvol;
            b3 = dN_dz[ii]*dvol;
            b4 = N[ii]*dvol;

            b5 = muTaf*b1;
            b6 = muTaf*b2;
            b7 = muTaf*b3;
            b8 = af*b4;

            Da = (vel(0)*b1 + vel(1)*b2 + vel(2)*b3)*tau[0]*rho;

            for(jj=0;jj<nlbf;jj++)
            {
              TJ   = 4*jj;
              TJp1 = TJ+1;
              TJp2 = TJ+2;
              TJp3 = TJ+3;

              fact2 = rho*acceFact*N(jj);

              // time acceleration term
              fact = b4*fact2 ;

              // diffusion term
              fact += (b5*dN_dx(jj)+b6*dN_dy(jj)+b7*dN_dz(jj));

              Klocal(TI,   TJ)   += fact;
              Klocal(TIp1, TJp1) += fact;
              Klocal(TIp2, TJp2) += fact;

              // convection term

              gradN = grad*(rho*N(jj));

              Db = rho*(vel(0)*dN_dx[jj] + vel(1)*dN_dy[jj] + vel(2)*dN_dz[jj] );

              gradN(0,0) += Db;
              gradN(1,1) += Db;
              gradN(2,2) += Db;
              gradN *= af;

              Klocal(TI,   TJ)   += (b4*gradN(0,0));
              Klocal(TI,   TJp1) += (b4*gradN(0,1));
              Klocal(TI,   TJp2) += (b4*gradN(0,2));

              Klocal(TIp1, TJ)   += (b4*gradN(1,0));
              Klocal(TIp1, TJp1) += (b4*gradN(1,1));
              Klocal(TIp1, TJp2) += (b4*gradN(1,2));

              Klocal(TIp2, TJ)   += (b4*gradN(2,0));
              Klocal(TIp2, TJp1) += (b4*gradN(2,1));
              Klocal(TIp2, TJp2) += (b4*gradN(2,2));

              // pressure term
              Klocal(TI,   TJp3) -= (b1*N(jj));
              Klocal(TIp1, TJp3) -= (b2*N(jj));
              Klocal(TIp2, TJp3) -= (b3*N(jj));

              // continuity equation
              Klocal(TIp3, TJ)   += (b8*dN_dx(jj));
              Klocal(TIp3, TJp1) += (b8*dN_dy(jj));
              Klocal(TIp3, TJp2) += (b8*dN_dz(jj));

              // SUPG and PSPG stabilisation terms
              //fact2 -= mu*d2N(jj);

              Dj(0,0) = gradN(0,0) + fact2;
              Dj(0,1) = gradN(0,1);
              Dj(0,2) = gradN(0,2);
              Dj(0,3) = dN_dx(jj);

              Dj(1,0) = gradN(1,0);
              Dj(1,1) = gradN(1,1) + fact2;
              Dj(1,2) = gradN(1,2);
              Dj(1,3) = dN_dy(jj);

              Dj(2,0) = gradN(2,0);
              Dj(2,1) = gradN(2,1);
              Dj(2,2) = gradN(2,2) + fact2;
              Dj(2,3) = dN_dz(jj);

              // SUPG
              Klocal(TI, TJ)     += Da*Dj(0,0);
              Klocal(TI, TJp1)   += Da*Dj(0,1);
              Klocal(TI, TJp2)   += Da*Dj(0,2);
              Klocal(TI, TJp3)   += Da*Dj(0,3);

              Klocal(TIp1, TJ)   += Da*Dj(1,0);
              Klocal(TIp1, TJp1) += Da*Dj(1,1);
              Klocal(TIp1, TJp2) += Da*Dj(1,2);
              Klocal(TIp1, TJp3) += Da*Dj(1,3);

              Klocal(TIp2, TJ)   += Da*Dj(2,0);
              Klocal(TIp2, TJp1) += Da*Dj(2,1);
              Klocal(TIp2, TJp2) += Da*Dj(2,2);
              Klocal(TIp2, TJp3) += Da*Dj(2,3);

              c1 = tau[0] * af * rho * N[jj];

              Klocal(TI,   TJ)   += ( c1 * b1 * rStab(0) );  
              Klocal(TI,   TJp1) += ( c1 * b2 * rStab(0) );
              Klocal(TI,   TJp2) += ( c1 * b3 * rStab(0) );

              Klocal(TIp1, TJ)   += ( c1 * b1 * rStab(1) );
              Klocal(TIp1, TJp1) += ( c1 * b2 * rStab(1) );
              Klocal(TIp1, TJp2) += ( c1 * b3 * rStab(1) );

              Klocal(TIp2, TJ)   += ( c1 * b1 * rStab(2) );
              Klocal(TIp2, TJp1) += ( c1 * b2 * rStab(2) );
              Klocal(TIp2, TJp2) += ( c1 * b3 * rStab(2) );

              // PSPG
              Klocal(TIp3, TJ)   += (b1*Dj(0,0) + b2*Dj(1,0) + b3*Dj(2,0))*tau[1];
              Klocal(TIp3, TJp1) += (b1*Dj(0,1) + b2*Dj(1,1) + b3*Dj(2,1))*tau[1];
              Klocal(TIp3, TJp2) += (b1*Dj(0,2) + b2*Dj(1,2) + b3*Dj(2,2))*tau[1];
              Klocal(TIp3, TJp3) += (b1*Dj(0,3) + b2*Dj(1,3) + b3*Dj(2,3))*tau[1];

              // LSIC stabilisation

              fact2 = rho*af*tau[2];

              Klocal(TI,   TJ)   += (b1*dN_dx(jj))*fact2;
              Klocal(TI,   TJp1) += (b1*dN_dy(jj))*fact2;
              Klocal(TI,   TJp2) += (b1*dN_dz(jj))*fact2;

              Klocal(TIp1, TJ)   += (b2*dN_dx(jj))*fact2;
              Klocal(TIp1, TJp1) += (b2*dN_dy(jj))*fact2;
              Klocal(TIp1, TJp2) += (b2*dN_dz(jj))*fact2;

              Klocal(TIp2, TJ)   += (b3*dN_dx(jj))*fact2;
              Klocal(TIp2, TJp1) += (b3*dN_dy(jj))*fact2;
              Klocal(TIp2, TJp2) += (b3*dN_dz(jj))*fact2;
            }

            Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) + b3*stress(0,2) );
            Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) + b3*stress(1,2) );
            Flocal(TIp2) -= (b4*res2(2) + b1*stress(2,0) + b2*stress(2,1) + b3*stress(2,2) );
            Flocal(TIp3) -= (b4*grad.trace());

            // SUPG stabilisation terms
            Flocal(TI)   -= Da*rStab(0);
            Flocal(TIp1) -= Da*rStab(1);
            Flocal(TIp2) -= Da*rStab(2);

            // PSPG stabilisation terms
            Flocal(TIp3) -= (tau[1]*(b1*rStab(0)+b2*rStab(1)+b3*rStab(2)));

            // LSIC stabilisation terms
            fact2 = tau[2]*rho*grad.trace();

            Flocal(TI)   -= b1*fact2;
            Flocal(TIp1) -= b2*fact2;
            Flocal(TIp2) -= b3*fact2;
          }
  }//gp
*/
  return 0;
}
//


/*
int LagrangeElem3DNavierStokesTetra4Node::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
{
  // Semi-implicit formulation - type A
  
  double  dvol, b1, b2, b3, b4, b5, b6, b7, b8, Da, Db, CI=4.0;
  double  Jac, totvol, pres, fact, fact2;
  double  param[3], bforce[3], tau[3];

  VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf), dN_dz(nlbf);
  VectorXd  res(4), res2(3), dp(3), Du(3), vel(3), force(3), gradTvel(3), rStab(3);
  VectorXd  velPrev(3), velExt(3), velTemp(3), velDot(3);
  MatrixXd  grad(3,3), gradN(3,3), stress(3,3),  Dj(3,4), forVolume(4,4), gradPrev(3,3);
  Dj.setZero();

  int   index, ii, jj, gp, TI, TIp1, TIp2, TIp3, TJ, TJp1, TJp2, TJp3;

  elmDat = &(SolnData->ElemProp[elmType].data[0]);

  double  rho = elmDat[3];
  double  mu  = elmDat[4];

  bforce[0]  = 0.0;
  bforce[1]  = 0.0;
  bforce[2]  = 0.0;

  //bforce[0]  = elmDat[6]*timeFunction[0].prop ;
  //bforce[1]  = elmDat[7]*timeFunction[0].prop ;
  //bforce[2]  = elmDat[7]*timeFunction[0].prop ;

  double  dt = mpapTime.dt;
  double  af = SolnData->td(2);
  double  am = SolnData->td(1);
  double  acceFact = am*SolnData->td(9);
  double  muTaf = mu*af;

  double  tCur = mpapTime.cur - (1.0-af)*dt;
  

  double xNode[4], yNode[4], zNode[4], xx, yy, zz;

  for(ii=0;ii<npElem;ii++)
  {
    xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
    yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
  }

  //forVolume(0,0) = 1.0; forVolume(0,1) = xNode[0]; forVolume(0,2) = yNode[0]; forVolume(0,3) = zNode[0];
  //forVolume(1,0) = 1.0; forVolume(1,1) = xNode[1]; forVolume(1,2) = yNode[1]; forVolume(1,3) = zNode[1];
  //forVolume(2,0) = 1.0; forVolume(2,1) = xNode[2]; forVolume(2,2) = yNode[2]; forVolume(2,3) = zNode[2];
  //forVolume(3,0) = 1.0; forVolume(3,1) = xNode[3]; forVolume(3,2) = yNode[3]; forVolume(3,3) = zNode[3];

  forVolume(0,0) = 1.0;      forVolume(0,1) = 1.0;      forVolume(0,2) = 1.0;      forVolume(0,3) = 1.0;
  forVolume(1,0) = xNode[0]; forVolume(1,1) = xNode[1]; forVolume(1,2) = xNode[2]; forVolume(1,3) = xNode[3];
  forVolume(2,0) = yNode[0]; forVolume(2,1) = yNode[1]; forVolume(2,2) = yNode[2]; forVolume(2,3) = yNode[3];
  forVolume(3,0) = zNode[0]; forVolume(3,1) = zNode[1]; forVolume(3,2) = zNode[2]; forVolume(3,3) = zNode[3];

  totvol = forVolume.determinant()/6.0;
  //cout << " totvol = " << totvol << endl;

  double h = pow(6.0*totvol/PI, 1.0/3.0);
  
  double stabParam = h*h/(4.0*mu);
  tau[0] = elmDat[8]*stabParam;  // SUPG
  tau[1] = elmDat[9]*stabParam;  // PSPG
  tau[2] = elmDat[10]*stabParam; // LSIC

  //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", xNode[0], xNode[1], xNode[2], xNode[3]);
  //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n\n\n ", yNode[0], yNode[1], yNode[2], yNode[3]);

  //nGP = (int) elmDat[0] ;
  nGP = 1 ;

  vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;

  getGaussPointsTetra(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);

  if(Klocal.rows() != nsize)
  {
    Klocal.resize(nsize, nsize);
    Flocal.resize(nsize);
  }
  Klocal.setZero();
  Flocal.setZero();

  //cout << " AAAAAAAAAA " << endl;

  for(gp=0;gp<nGP;gp++)
  {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];
        param[2] = gausspoints3[gp];

        GeomData->computeBasisFunctions3D(0, 1, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        dvol = Jac * gaussweights[gp];

          vel(0) = computeValueCur(0, N);
          vel(1) = computeValueCur(1, N);
          vel(2) = computeValueCur(2, N);

          velPrev(0) = computeValuePrev(0, N);
          velPrev(1) = computeValuePrev(1, N);
          velPrev(2) = computeValuePrev(2, N);

          velExt(0) = computeValueExtrap(0, N);
          velExt(1) = computeValueExtrap(1, N);
          velExt(2) = computeValueExtrap(2, N);

          grad(0,0) = computeValueCur(0, dN_dx);
          grad(0,1) = computeValueCur(0, dN_dy);
          grad(0,2) = computeValueCur(0, dN_dz);

          grad(1,0) = computeValueCur(1, dN_dx);
          grad(1,1) = computeValueCur(1, dN_dy);
          grad(1,2) = computeValueCur(1, dN_dz);

          grad(2,0) = computeValueCur(2, dN_dx);
          grad(2,1) = computeValueCur(2, dN_dy);
          grad(2,2) = computeValueCur(2, dN_dz);

          gradPrev(0,0) = computeValuePrev(0, dN_dx);
          gradPrev(0,1) = computeValuePrev(0, dN_dy);
          gradPrev(0,2) = computeValuePrev(0, dN_dz);

          gradPrev(1,0) = computeValuePrev(1, dN_dx);
          gradPrev(1,1) = computeValuePrev(1, dN_dy);
          gradPrev(1,2) = computeValuePrev(1, dN_dz);

          gradPrev(2,0) = computeValuePrev(2, dN_dx);
          gradPrev(2,1) = computeValuePrev(2, dN_dy);
          gradPrev(2,2) = computeValuePrev(2, dN_dz);

          Du.setZero();

          pres   = computeValue(3, N);
          dp(0)  = computeValue(3, dN_dx);
          dp(1)  = computeValue(3, dN_dy);
          dp(2)  = computeValue(3, dN_dz);

          velDot(0) = computeValueDotCur(0, N);
          velDot(1) = computeValueDotCur(1, N);
          velDot(2) = computeValueDotCur(2, N);

          // this is pseudo-stress
          stress = mu*grad;
          stress(0,0) -= pres;
          stress(1,1) -= pres;
          stress(2,2) -= pres;

          //cout << xx << '\t' << yy << endl;

          force.setZero();
          //force(0) = analy.computeForce(0, xx, yy, tCur);
          //force(1) = analy.computeForce(1, xx, yy, tCur);
          //cout << force(0) << '\t' << force(1) << endl;

          velExt = af*velExt + (1.0-af)*velPrev;

          gradTvel = grad*velExt;

          res2(0) = rho*(velDot(0) + gradTvel(0) - force(0)) ;
          res2(1) = rho*(velDot(1) + gradTvel(1) - force(1)) ;
          res2(2) = rho*(velDot(2) + gradTvel(2) - force(2)) ;

          rStab(0) = res2(0) - mu*Du(0) + dp(0) ;
          rStab(1) = res2(1) - mu*Du(1) + dp(1) ;
          rStab(2) = res2(2) - mu*Du(2) + dp(2) ;

          velTemp(0) = velPrev(0);
          velTemp(1) = velPrev(1);
          velTemp(2) = velPrev(2);

          //evaluateStabParams_algo2(&velTemp(0), h, rho, mu, dt,  beta, tau);

          //evaluateStabParams_algo3(velTemp, matG, dt, rho, mu, CI, tau);

          //tau[0] *= elmDat[8];  // SUPG
          //tau[1] *= elmDat[9];  // PSPG
          //tau[2] *= elmDat[10]; // LSIC


          for(ii=0;ii<nlbf;ii++)
          {
            TI   = 4*ii;
            TIp1 = TI+1;
            TIp2 = TI+2;
            TIp3 = TI+3;

            b1 = dN_dx[ii]*dvol;
            b2 = dN_dy[ii]*dvol;
            b3 = dN_dz[ii]*dvol;
            b4 = N[ii]*dvol;

            b5 = muTaf*b1;
            b6 = muTaf*b2;
            b7 = muTaf*b3;
            b8 = af*b4;

            Da = (velExt(0)*b1 + velExt(1)*b2 + velExt(2)*b3)*tau[0];

            for(jj=0;jj<nlbf;jj++)
            {
              TJ   = 4*jj;
              TJp1 = TJ+1;
              TJp2 = TJ+2;
              TJp3 = TJ+3;

              fact2 = rho*acceFact*N(jj);

              // time acceleration term
              fact = b4*fact2 ;

              // diffusion term
              fact += (b5*dN_dx(jj)+b6*dN_dy(jj)+b7*dN_dz(jj));

              Klocal(TI,   TJ)   += fact;
              Klocal(TIp1, TJp1) += fact;
              Klocal(TIp2, TJp2) += fact;

              // convection term

              Db = rho*(velExt(0)*dN_dx(jj) + velExt(1)*dN_dy(jj) + velExt(2)*dN_dz(jj));

              Klocal(TI,   TJ)   += (b8*Db);
              Klocal(TIp1, TJp1) += (b8*Db);
              Klocal(TIp2, TJp2) += (b8*Db);

              // pressure term
              Klocal(TI,   TJp3) -= (af*b1*N(jj));
              Klocal(TIp1, TJp3) -= (af*b2*N(jj));
              Klocal(TIp2, TJp3) -= (af*b3*N(jj));

              // continuity equation
              Klocal(TIp3, TJ)   += (b8*dN_dx(jj));
              Klocal(TIp3, TJp1) += (b8*dN_dy(jj));
              Klocal(TIp3, TJp2) += (b8*dN_dz(jj));

              // SUPG and PSPG stabilisation terms
              //fact2 -= mu*d2N(jj);
              //fact2 -= 0.0;

              Db *= af;

              Dj(0,0) = Db + fact2;
              Dj(0,1) = 0.0;
              Dj(0,2) = 0.0;
              Dj(0,3) = af*dN_dx(jj);

              Dj(1,0) = 0.0;
              Dj(1,1) = Db + fact2;
              Dj(1,2) = 0.0;
              Dj(1,3) = af*dN_dy(jj);

              Dj(2,0) = 0.0;
              Dj(2,1) = 0.0;
              Dj(2,2) = Db + fact2;
              Dj(2,3) = af*dN_dz(jj);

              // SUPG
              Klocal(TI, TJ)     += Da*Dj(0,0);
              Klocal(TI, TJp1)   += Da*Dj(0,1);
              Klocal(TI, TJp2)   += Da*Dj(0,2);
              Klocal(TI, TJp3)   += Da*Dj(0,3);

              Klocal(TIp1, TJ)   += Da*Dj(1,0);
              Klocal(TIp1, TJp1) += Da*Dj(1,1);
              Klocal(TIp1, TJp2) += Da*Dj(1,2);
              Klocal(TIp1, TJp3) += Da*Dj(1,3);

              Klocal(TIp2, TJ)   += Da*Dj(2,0);
              Klocal(TIp2, TJp1) += Da*Dj(2,1);
              Klocal(TIp2, TJp2) += Da*Dj(2,2);
              Klocal(TIp2, TJp3) += Da*Dj(2,3);

              // PSPG
              Klocal(TIp3, TJ)   += (b1*Dj(0,0) + b2*Dj(1,0) + b3*Dj(2,0))*tau[1];
              Klocal(TIp3, TJp1) += (b1*Dj(0,1) + b2*Dj(1,1) + b3*Dj(2,1))*tau[1];
              Klocal(TIp3, TJp2) += (b1*Dj(0,2) + b2*Dj(1,2) + b3*Dj(2,2))*tau[1];
              Klocal(TIp3, TJp3) += (b1*Dj(0,3) + b2*Dj(1,3) + b3*Dj(2,3))*tau[1];

              // LSIC stabilisation

              fact2 = rho*af*tau[2];
              Klocal(TI,   TJ)   += (b1*dN_dx(jj))*fact2;
              Klocal(TI,   TJp1) += (b1*dN_dy(jj))*fact2;
              Klocal(TI,   TJp2) += (b1*dN_dz(jj))*fact2;

              Klocal(TIp1, TJ)   += (b2*dN_dx(jj))*fact2;
              Klocal(TIp1, TJp1) += (b2*dN_dy(jj))*fact2;
              Klocal(TIp1, TJp2) += (b2*dN_dz(jj))*fact2;

              Klocal(TIp2, TJ)   += (b3*dN_dx(jj))*fact2;
              Klocal(TIp2, TJp1) += (b3*dN_dy(jj))*fact2;
              Klocal(TIp2, TJp2) += (b3*dN_dz(jj))*fact2;
            }

            Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) + b3*stress(0,2) );
            Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) + b3*stress(1,2) );
            Flocal(TIp2) -= (b4*res2(2) + b1*stress(2,0) + b2*stress(2,1) + b3*stress(2,2) );
            Flocal(TIp3) -= (b4*grad.trace());

            // SUPG stabilisation terms
            Flocal(TI)   -= Da*rStab(0);
            Flocal(TIp1) -= Da*rStab(1);
            Flocal(TIp2) -= Da*rStab(2);

            // PSPG stabilisation terms
            Flocal(TIp3) -= (tau[1]*(b1*rStab(0)+b2*rStab(1)+b3*rStab(2)));

            // LSIC stabilisation terms

            fact2 = tau[2]*rho*grad.trace();

            Flocal(TI)   -= b1*fact2;
            Flocal(TIp1) -= b2*fact2;
            Flocal(TIp2) -= b3*fact2;
          }
  }//gp

  return 0;
}
*/




/*
int LagrangeElem3DNavierStokesTetra4Node::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
{
  // Semi-implicit formulation - type B
  
  double  dvol, b1, b2, b3, b4, b5, b6, b7, b8, Da, Db, CI=4.0;
  double  Jac, totvol, pres, fact, fact2;
  double  param[3], bforce[3], tau[3];

  VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf), dN_dz(nlbf);
  VectorXd  res(4), res2(3), dp(3), Du(3), vel(3), force(3), gradTvel(3), rStab(3);
  VectorXd  velPrev(3), velExt(3), velTemp(3), velDot(3);
  MatrixXd  grad(3,3), gradN(3,3), stress(3,3),  Dj(3,4), forVolume(4,4), gradPrev(3,3);
  Dj.setZero();


  int   index, ii, jj, gp, TI, TIp1, TIp2, TIp3, TJ, TJp1, TJp2, TJp3;

  elmDat = &(SolnData->ElemProp[elmType].data[0]);

  double  rho = elmDat[3];
  double  mu  = elmDat[4];

  bforce[0]  = 0.0;
  bforce[1]  = 0.0;
  bforce[2]  = 0.0;

  //bforce[0]  = elmDat[6]*timeFunction[0].prop ;
  //bforce[1]  = elmDat[7]*timeFunction[0].prop ;
  //bforce[2]  = elmDat[7]*timeFunction[0].prop ;


  double  dt = mpapTime.dt;
  double  af = SolnData->td(2);
  double  am = SolnData->td(1);
  double  acceFact = am*SolnData->td(9);
  double  muTaf = mu*af;

  double  tCur = mpapTime.cur - (1.0-af)*dt;
  

  double xNode[4], yNode[4], zNode[4], xx, yy, zz;

  for(ii=0;ii<npElem;ii++)
  {
    xNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][0];
    yNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][1];
    zNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][2];
  }


  //forVolume(0,0) = 1.0; forVolume(0,1) = xNode[0]; forVolume(0,2) = yNode[0]; forVolume(0,3) = zNode[0];
  //forVolume(1,0) = 1.0; forVolume(1,1) = xNode[1]; forVolume(1,2) = yNode[1]; forVolume(1,3) = zNode[1];
  //forVolume(2,0) = 1.0; forVolume(2,1) = xNode[2]; forVolume(2,2) = yNode[2]; forVolume(2,3) = zNode[2];
  //forVolume(3,0) = 1.0; forVolume(3,1) = xNode[3]; forVolume(3,2) = yNode[3]; forVolume(3,3) = zNode[3];

  forVolume(0,0) = 1.0;      forVolume(0,1) = 1.0;      forVolume(0,2) = 1.0;      forVolume(0,3) = 1.0;
  forVolume(1,0) = xNode[0]; forVolume(1,1) = xNode[1]; forVolume(1,2) = xNode[2]; forVolume(1,3) = xNode[3];
  forVolume(2,0) = yNode[0]; forVolume(2,1) = yNode[1]; forVolume(2,2) = yNode[2]; forVolume(2,3) = yNode[3];
  forVolume(3,0) = zNode[0]; forVolume(3,1) = zNode[1]; forVolume(3,2) = zNode[2]; forVolume(3,3) = zNode[3];

  totvol = forVolume.determinant()/6.0;
  //cout << " totvol = " << totvol << endl;

  //double  h2  = totvol*4.0/PI;
  double h = pow(3.0*totvol/PI/4.0, 1.0/3.0);
  //double h = 1.0/6.0;
  //h = sqrt(3.0*h*h);
  
  //double  alpha = 1.0/12.0, stabParam;

  //tau = elmDat[8]*alpha*h*h/sqrt(mu);
  //tau = elmDat[8]*alpha*h*h/mu/dt;

    //cout << " volume = " << volume << endl;
  //double  h2 = 4.0*totvol/PI;

  double stabParam = h*h/(12.0*mu);
  tau[0] = elmDat[8]*stabParam;  // SUPG
  tau[1] = elmDat[9]*stabParam;  // PSPG
  tau[2] = elmDat[10]*stabParam; // LSIC

  //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", xNode[0], xNode[1], xNode[2], xNode[3]);
  //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n\n\n ", yNode[0], yNode[1], yNode[2], yNode[3]);



  nGP = (int) elmDat[0] ;

  vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;

  getGaussPointsTetra(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);


    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }
    Klocal.setZero();
    Flocal.setZero();

  //cout << " AAAAAAAAAA " << endl;
  //cout << nGP1 << '\t' << nGP2 << endl;

  for(gp=0;gp<nGP;gp++)
  {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];
        param[2] = gausspoints3[gp];

        GeomData->computeBasisFunctions3D(0, 1, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        dvol = Jac * gaussweights[gp];

          vel(0) = computeValueCur(0, N);
          vel(1) = computeValueCur(1, N);
          vel(2) = computeValueCur(2, N);

          velPrev(0) = computeValuePrev(0, N);
          velPrev(1) = computeValuePrev(1, N);
          velPrev(2) = computeValuePrev(2, N);

          velExt(0) = computeValueExtrap(0, N);
          velExt(1) = computeValueExtrap(1, N);
          velExt(2) = computeValueExtrap(2, N);

          grad(0,0) = computeValueCur(0, dN_dx);
          grad(0,1) = computeValueCur(0, dN_dy);
          grad(0,2) = computeValueCur(0, dN_dz);

          grad(1,0) = computeValueCur(1, dN_dx);
          grad(1,1) = computeValueCur(1, dN_dy);
          grad(1,2) = computeValueCur(1, dN_dz);

          grad(2,0) = computeValueCur(2, dN_dx);
          grad(2,1) = computeValueCur(2, dN_dy);
          grad(2,2) = computeValueCur(2, dN_dz);

          gradPrev(0,0) = computeValuePrev(0, dN_dx);
          gradPrev(0,1) = computeValuePrev(0, dN_dy);
          gradPrev(0,2) = computeValuePrev(0, dN_dz);

          gradPrev(1,0) = computeValuePrev(1, dN_dx);
          gradPrev(1,1) = computeValuePrev(1, dN_dy);
          gradPrev(1,2) = computeValuePrev(1, dN_dz);

          gradPrev(2,0) = computeValuePrev(2, dN_dx);
          gradPrev(2,1) = computeValuePrev(2, dN_dy);
          gradPrev(2,2) = computeValuePrev(2, dN_dz);

          Du.setZero();

          pres   = computeValue(3, N);
          dp(0)  = computeValue(3, dN_dx);
          dp(1)  = computeValue(3, dN_dy);
          dp(2)  = computeValue(3, dN_dz);

          velDot(0) = computeValueDotCur(0, N);
          velDot(1) = computeValueDotCur(1, N);
          velDot(2) = computeValueDotCur(2, N);

          // this is pseudo-stress
          stress = mu*grad;
          stress(0,0) -= pres;
          stress(1,1) -= pres;
          stress(2,2) -= pres;

          //cout << xx << '\t' << yy << endl;

          force.setZero();
          //force(0) = analy.computeForce(0, xx, yy, tCur);
          //force(1) = analy.computeForce(1, xx, yy, tCur);
          //cout << force(0) << '\t' << force(1) << endl;

          gradTvel = gradPrev*vel + grad*velPrev - gradPrev*velPrev;

          res2(0) = rho*(velDot(0) + gradTvel(0) - force(0)) ;
          res2(1) = rho*(velDot(1) + gradTvel(1) - force(1)) ;
          res2(2) = rho*(velDot(2) + gradTvel(2) - force(2)) ;

          rStab(0) = res2(0) - mu*Du(0) + dp(0) ;
          rStab(1) = res2(1) - mu*Du(1) + dp(1) ;
          rStab(2) = res2(2) - mu*Du(2) + dp(2) ;

          velTemp(0) = velPrev(0);
          velTemp(1) = velPrev(1);
          velTemp(2) = velPrev(2);

          //evaluateStabParams_algo2(&velTemp(0), h, rho, mu, dt,  beta, tau);

          //evaluateStabParams_algo3(velTemp, matG, dt, rho, mu, CI, tau);

          //tau[0] *= elmDat[8];  // SUPG
          //tau[1] *= elmDat[9];  // PSPG
          //tau[2] *= elmDat[10]; // LSIC


          for(ii=0;ii<nlbf;ii++)
          {
            TI   = 4*ii;
            TIp1 = TI+1;
            TIp2 = TI+2;
            TIp3 = TI+3;

            b1 = dN_dx[ii]*dvol;
            b2 = dN_dy[ii]*dvol;
            b3 = dN_dz[ii]*dvol;
            b4 = N[ii]*dvol;

            b5 = muTaf*b1;
            b6 = muTaf*b2;
            b7 = muTaf*b3;
            b8 = af*b4;

            Da = rho*(velPrev(0)*b1 + velPrev(1)*b2 + velPrev(2)*b3)*tau[0];

            for(jj=0;jj<nlbf;jj++)
            {
              TJ   = 4*jj;
              TJp1 = TJ+1;
              TJp2 = TJ+2;
              TJp3 = TJ+3;

              fact2 = rho*acceFact*N(jj);

              // time acceleration term
              fact = b4*fact2 ;

              // diffusion term
              fact += (b5*dN_dx(jj)+b6*dN_dy(jj)+b7*dN_dz(jj));

              Klocal(TI,   TJ)   += fact;
              Klocal(TIp1, TJp1) += fact;
              Klocal(TIp2, TJp2) += fact;

              // convection term - semi-implicit type B

              gradN = gradPrev*(rho*N(jj));

              Db = rho*(velPrev(0)*dN_dx(jj) + velPrev(1)*dN_dy(jj) + velPrev(2)*dN_dz(jj));

              gradN(0,0) += Db;
              gradN(1,1) += Db;
              gradN(2,2) += Db;

              Klocal(TI,   TJ)   += (b8*gradN(0,0));
              Klocal(TI,   TJp1) += (b8*gradN(0,1));
              Klocal(TI,   TJp2) += (b8*gradN(0,2));

              Klocal(TIp1, TJ)   += (b8*gradN(1,0));
              Klocal(TIp1, TJp1) += (b8*gradN(1,1));
              Klocal(TIp1, TJp2) += (b8*gradN(1,2));

              Klocal(TIp2, TJ)   += (b8*gradN(2,0));
              Klocal(TIp2, TJp1) += (b8*gradN(2,1));
              Klocal(TIp2, TJp2) += (b8*gradN(2,2));

              // pressure term
              Klocal(TI,   TJp3) -= (b1*N(jj));
              Klocal(TIp1, TJp3) -= (b2*N(jj));
              Klocal(TIp2, TJp3) -= (b3*N(jj));

              // continuity equation
              Klocal(TIp3, TJ)   += (b8*dN_dx(jj));
              Klocal(TIp3, TJp1) += (b8*dN_dy(jj));
              Klocal(TIp3, TJp2) += (b8*dN_dz(jj));

              // SUPG and PSPG stabilisation terms
              //fact2 -= mu*d2N(jj);
              //fact2 -= 0.0;

              gradN *= af;

              Dj(0,0) = gradN(0,0) + fact2;
              Dj(0,1) = gradN(0,1);
              Dj(0,2) = gradN(0,2);
              Dj(0,3) = af*dN_dx(jj);

              Dj(1,0) = gradN(1,0);
              Dj(1,1) = gradN(1,1) + fact2;
              Dj(1,2) = gradN(1,2);
              Dj(1,3) = af*dN_dy(jj);

              Dj(2,0) = gradN(2,0);
              Dj(2,1) = gradN(2,1);
              Dj(2,2) = gradN(2,2) + fact2;
              Dj(2,3) = af*dN_dz(jj);

              // SUPG
              Klocal(TI, TJ)     += Da*Dj(0,0);
              Klocal(TI, TJp1)   += Da*Dj(0,1);
              Klocal(TI, TJp2)   += Da*Dj(0,2);
              Klocal(TI, TJp3)   += Da*Dj(0,3);

              Klocal(TIp1, TJ)   += Da*Dj(1,0);
              Klocal(TIp1, TJp1) += Da*Dj(1,1);
              Klocal(TIp1, TJp2) += Da*Dj(1,2);
              Klocal(TIp1, TJp3) += Da*Dj(1,3);

              Klocal(TIp2, TJ)   += Da*Dj(2,0);
              Klocal(TIp2, TJp1) += Da*Dj(2,1);
              Klocal(TIp2, TJp2) += Da*Dj(2,2);
              Klocal(TIp2, TJp3) += Da*Dj(2,3);

              // PSPG
              Klocal(TIp3, TJ)   += (b1*Dj(0,0) + b2*Dj(1,0) + b3*Dj(2,0))*tau[1];
              Klocal(TIp3, TJp1) += (b1*Dj(0,1) + b2*Dj(1,1) + b3*Dj(2,1))*tau[1];
              Klocal(TIp3, TJp2) += (b1*Dj(0,2) + b2*Dj(1,2) + b3*Dj(2,2))*tau[1];
              Klocal(TIp3, TJp3) += (b1*Dj(0,3) + b2*Dj(1,3) + b3*Dj(2,3))*tau[1];

              // LSIC stabilisation

              fact2 = rho*af*tau[2];

              Klocal(TI,   TJ)   += (b1*dN_dx(jj))*fact2;
              Klocal(TI,   TJp1) += (b1*dN_dy(jj))*fact2;
              Klocal(TI,   TJp2) += (b1*dN_dz(jj))*fact2;

              Klocal(TIp1, TJ)   += (b2*dN_dx(jj))*fact2;
              Klocal(TIp1, TJp1) += (b2*dN_dy(jj))*fact2;
              Klocal(TIp1, TJp2) += (b2*dN_dz(jj))*fact2;

              Klocal(TIp2, TJ)   += (b3*dN_dx(jj))*fact2;
              Klocal(TIp2, TJp1) += (b3*dN_dy(jj))*fact2;
              Klocal(TIp2, TJp2) += (b3*dN_dz(jj))*fact2;
            }

            Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) + b3*stress(0,2) );
            Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) + b3*stress(1,2) );
            Flocal(TIp2) -= (b4*res2(2) + b1*stress(2,0) + b2*stress(2,1) + b3*stress(2,2) );
            Flocal(TIp3) -= (b4*grad.trace());

            // SUPG stabilisation terms
            Flocal(TI)   -= Da*rStab(0);
            Flocal(TIp1) -= Da*rStab(1);
            Flocal(TIp2) -= Da*rStab(2);

            // PSPG stabilisation terms
            Flocal(TIp3) -= (tau[1]*(b1*rStab(0)+b2*rStab(1)+b3*rStab(2)));

            // LSIC stabilisation terms
            fact2 = tau[2]*rho*grad.trace();

            Flocal(TI)   -= b1*fact2;
            Flocal(TIp1) -= b2*fact2;
            Flocal(TIp2) -= b3*fact2;
          }
  }//gp

  return 0;
}
*/




int  LagrangeElem3DNavierStokesTetra4Node::calcError(int index)
{
/*
    //cout << " Error ... " << endl;

    // computing error for Navier-Stokes
    ///////////////////////////////////////////////////////////

    //Stokes2DEx1  analy;
    //Stokes2DEx2  analy;

    //Stokes2DEx3  analy;
    //Kovasznay  analy;
    //analy.SetPressure(0.0);

    //PearsonVortex analy;

    double rho = elmDat[4];
    double mu  = elmDat[5];

    KimMoinFlowUnsteadyNavierStokes  analy(rho, mu, 1.0);
    //Kovasznay  analy;

    int      ii, jj, gp, TI;
    double   Jac, dvol, diff, fact, val;
    VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf);
    double  param[2], geom[2], v[2];
    MatrixXd  F(2,2);

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsQuad(nGP, gausspoints1, gausspoints2, gaussweights);

    elemError = 0.0;

    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        //GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_QUAD, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        dvol = gaussweights[gp] * Jac;

        //geom[0] = computeGeomOrig(0, N);
        //geom[1] = computeGeomOrig(1, N);

        if(index < 3) // L2 norm in x-velocity (index=0), y-velocity (index=1) and pressure (index=2)
        {
            //val = analy.computeValue(index, geom[0], geom[1], 0.0, mpapTime.cur);

            //val -= computeValue(index, N);

            //printf(" computed, exact and difference \t %12.8f \t %12.8f \t %12.8f  \t %12.8f   \t %12.8f \n", vx, vx2, vy, vy2, diff);

            elemError += ( (val*val) * dvol );
        }
        else
        {
            //v[0] = analy.computeValue(0, geom[0], geom[1], 0.0, mpapTime.cur);
            //v[1] = analy.computeValue(1, geom[0], geom[1], 0.0, mpapTime.cur);

            //v[0] -= computeValue(0, N);
            //v[1] -= computeValue(1, N);

            F.setZero();
            //analy.computeDerivatives(geom[0], geom[1], 0.0, mpapTime.cur, &(F(0,0)));

            //F(0,0) -= computeValue(0, dN_dx);
            //F(0,1) -= computeValue(0, dN_dy);
            //F(1,0) -= computeValue(1, dN_dx);
            //F(1,1) -= computeValue(1, dN_dy);

            val = v[0]*v[0]+v[1]*v[1] + F(0,0)*F(0,0)+F(0,1)*F(0,1)+F(1,0)*F(1,0)+F(1,1)*F(1,1);

            //printf(" computed, exact and difference \t %12.8f \t %12.8f \t %12.8f \n", dx, dy, diff);

            elemError += ( val * dvol );
        }
    }//gp
*/
    return 0;
}



