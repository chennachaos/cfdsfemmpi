/*
 * 3    2
 * 
 * 0    1
*/


#include "LagrangeElem2DNavierStokesQuad4Node.h"
#include "SolutionData.h"
#include "UtilitiesQuadrature.h"
#include "stabilisationRoutines.h"
#include "KimMoinFlow.h"
#include "BasisFunctionsLagrange.h"


using namespace std;


LagrangeElem2DNavierStokesQuad4Node::LagrangeElem2DNavierStokesQuad4Node()
{
  degree = 1;
  npElem = 4;
  nlbf   = 4;
  ndof   = 3;
  nsize  = nlbf*ndof;

  nGP = 4;
}


LagrangeElem2DNavierStokesQuad4Node::~LagrangeElem2DNavierStokesQuad4Node()
{
}


void LagrangeElem2DNavierStokesQuad4Node::prepareElemData(vector<vector<double> >& node_coords)
{
    // compute Volume and basis function derivatives

    double  dvol, Jac, param[2];

    int   ii, gp;

    double xNode[npElem], yNode[npElem], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = node_coords[SolnData->node_map_get_old[nodeNums[ii]]][0];
      yNode[ii] = node_coords[SolnData->node_map_get_old[nodeNums[ii]]][1];

      //xNode[ii] = node_coords[nodeNums[ii]][0];
      //yNode[ii] = node_coords[nodeNums[ii]][1];
    }

    // Gauss point coordinates and weights
    vector<double>  gpts1(nGP), gpts2(nGP), gwts(nGP);
    getGaussPointsQuad(nGP, gpts1, gpts2, gwts);

    Nv.resize(nGP);
    dNvdx.resize(nGP);
    dNvdy.resize(nGP);

    for(int ii=0; ii<nGP; ii++)
    {
      Nv[ii].resize(nlbf);
      dNvdx[ii].resize(nlbf);
      dNvdy[ii].resize(nlbf);

      Nv[ii].setZero();
      dNvdx[ii].setZero();
      dNvdy[ii].setZero();
    }

    //printVector(nodeNums);

    elemVolGP.resize(nGP);

    elemVol=0.0;
    for(gp=0; gp<nGP; gp++)
    {
      param[0] = gpts1[gp];
      param[1] = gpts2[gp];

      computeBasisFunctions2D(false, ELEM_SHAPE_QUAD, degree, param, xNode, yNode, &Nv[gp][0], &dNvdx[gp][0], &dNvdy[gp][0], Jac);

      //cout << " Jac = " << Jac << endl;
      if(Jac < 0.0)
      {
        cout << " Negative Jacobian in 'LagrangeElem2DNavierStokesQuad4Node::prepareElemData' " << endl;
        cout << " Jac = " << Jac << endl;
        exit(1);
      }

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
    }

    return;
}



/*
int LagrangeElem2DNavierStokesQuad4Node::calcStiffnessAndResidual(vector<vector<double> >& node_coords, double* elemData, double* timeData, VectorXd& solnPrev, VectorXd& solnPrev2, VectorXd& solnCur, VectorXd& solnDotCur, MatrixXd& Klocal, VectorXd& Flocal, double timeCur)
{
    // Fully-implicit formulation

    int ii, jj, gp, TI, TIp1, TIp2, count, TJ, TJp1, TJp2;

    double  dvol, b1, b2, b3, b4, b5, b6, b7, b8, xx, yy;
    double  pres, Da, Db, rad, urdr, urdr2, tau[3], CI=4.0;
    double  fact, fact1, fact2, param[2];

    double  xNode[4], yNode[4];
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = node_coords[nodeNums[ii]][0];
      yNode[ii] = node_coords[nodeNums[ii]][1];
    }

    VectorXd  res(3), res2(2), dp(2), Du(2), vel(2), velPrev(2), velExtp(2), velDot(2), force(2), gradTvel(2), rStab(3);
    MatrixXd  Dj(2, 3), grad(2,2), gradN(2,2), stress(2,2), gradPrev(2,2);
    Dj.setZero();
    VectorXd  velTemp(3);


    MatrixXd  matG(3,3);
    matG.setZero();
    fact = 4.0/charlen/charlen;
    matG(0,0) = fact;
    matG(1,1) = fact;

    double  rho = elemData[0];
    double  mu  = elemData[1];

    double  am = timeData[1];
    double  af = timeData[2];
    double  acceFact = timeData[8];
    double  muTaf = mu*af;
    double  dt = timeData[5]/af;

    //KimMoinFlowUnsteadyNavierStokes  analy(rho, mu, 1.0);

    bool axsy = false;
    //axsy = ((int)elemData[2] == 1);

    double  tCur  = timeCur;
    double  tPrev = tCur - dt;


    Klocal.setZero();
    Flocal.setZero();

    for(gp=0; gp<nGP; gp++)
    {
        // compute the gradient of velocity
        xx = 0.0; yy = 0.0;
        vel[0] = vel[1] = 0.0;
        velPrev[0] = velPrev[1] = 0.0;
        velExtp[0] = velExtp[1] = 0.0;
        velDot[0] = velDot[1] = 0.0;
        grad.setZero();
        Du.setZero();
        pres = 0.0;
        dp.setZero();

        for(ii=0; ii<npElem; ii++)
        {
            xx += xNode[ii]*Nv[gp][ii];
            yy += yNode[ii]*Nv[gp][ii];

            TI   = nodeNums[ii]*ndof;
            TIp1 = TI+1;
            TIp2 = TI+2;

            b1 = solnPrev[TI];
            b2 = solnPrev[TIp1];

            velPrev[0] += b1*Nv[gp][ii];
            velPrev[1] += b2*Nv[gp][ii];

            b1 = solnCur[TI];
            b2 = solnCur[TIp1];
            b3 = solnCur[TIp2];

            vel[0]     += b1*Nv[gp][ii];
            vel[1]     += b2*Nv[gp][ii];
            pres       += b3*Nv[gp][ii];

            grad(0,0)  += b1*dNvdx[gp][ii];
            grad(0,1)  += b1*dNvdy[gp][ii];
            grad(1,0)  += b2*dNvdx[gp][ii];
            grad(1,1)  += b2*dNvdy[gp][ii];

            velDot[0] += solnDotCur[TI]*Nv[gp][ii];
            velDot[1] += solnDotCur[TIp1]*Nv[gp][ii];

            dp(0)     += b3*dNvdx[gp][ii];
            dp(1)     += b3*dNvdy[gp][ii];
        }

        // this is pseudo-stress
        stress = mu*grad;
        stress(0,0) -= pres;
        stress(1,1) -= pres;

        force.setZero();
        //force(0) = af*analy.computeForce(0, geom[0], geom[1], 0.0, tCur)+(1.0-af)*analy.computeForce(0, geom[0], geom[1], 0.0, tPrev);
        //force(1) = af*analy.computeForce(1, geom[0], geom[1], 0.0, tCur)+(1.0-af)*analy.computeForce(1, geom[0], geom[1], 0.0, tPrev);
        //cout << force(0) << '\t' << force(1) << endl;


        dvol = elemVolGP[gp];
        if(axsy)
        {
          rad = xx;

          urdr  = vel(0)/rad;
          urdr2 = urdr/rad;
          dvol *= (2.0*PI*rad);
        }

        gradTvel = grad*vel;

        res2(0) = rho*(velDot(0) + gradTvel(0) - force(0)) ;
        res2(1) = rho*(velDot(1) + gradTvel(1) - force(1)) ;

        rStab(0) = res2(0) - mu*Du(0) + dp(0) ;
        rStab(1) = res2(1) - mu*Du(1) + dp(1) ;

        velTemp(0) = velPrev(0);
        velTemp(1) = velPrev(1);
        velTemp(2) = 0.0;

        //evaluateStabParams_algo2(&velTemp(0), h, rho, mu, dt,  beta, tau);
        evaluateStabParams_algo3(velTemp, matG, dt, rho, mu, CI, tau);

        tau[0] *= elemData[8];                             // SUPG
        tau[1] *= elemData[9];                             // PSPG
        tau[2] *= elemData[10];                            // LSIC

        //cout << tau[0] << '\t' << tau[1] << '\t' << tau[2] << endl;

        for(ii=0;ii<nlbf;ii++)
        {
            TI   = ndof*ii;
            TIp1 = TI+1;
            TIp2 = TI+2;

            b1 = dNvdx[gp][ii]*dvol;
            b2 = dNvdy[gp][ii]*dvol;
            b4 = Nv[gp][ii]*dvol;

            b5 = muTaf*b1;
            b6 = muTaf*b2;
            b8 = af*b4;

            Da = (vel(0)*b1 + vel(1)*b2)*tau[0];

            for(jj=0;jj<nlbf;jj++)
            {
                TJ   = ndof*jj;
                TJp1 = TJ+1;
                TJp2 = TJ+2;

                fact2 = rho*acceFact*Nv[gp][jj];

                // time acceleration term
                fact = b4*fact2 ;

                // diffusion term
                fact += b5*dNvdx[gp][jj]+b6*dNvdy[gp][jj];

                Klocal(TI,   TJ)   += fact;
                Klocal(TIp1, TJp1) += fact;

                // convection term

                gradN = grad*(rho*Nv[gp][jj]);

                Db = rho*(vel(0)*dNvdx[gp][jj] + vel(1)*dNvdy[gp][jj]);

                gradN(0,0) += Db;
                gradN(1,1) += Db;

                Klocal(TI,   TJ)   += (b8*gradN(0,0));
                Klocal(TI,   TJp1) += (b8*gradN(0,1));
                Klocal(TIp1, TJ)   += (b8*gradN(1,0));
                Klocal(TIp1, TJp1) += (b8*gradN(1,1));

                // pressure term
                Klocal(TI,   TJp2) -= (b1*af*Nv[gp][jj]);
                Klocal(TIp1, TJp2) -= (b2*af*Nv[gp][jj]);

                // continuity equation
                Klocal(TIp2, TJ)   += (b8*dNvdx[gp][jj]);
                Klocal(TIp2, TJp1) += (b8*dNvdy[gp][jj]);
                Klocal(TIp2, TJp2) += 0.0;

                // SUPG and PSPG stabilisation terms
                //fact2 -= mu*d2Nv[gp][jj];
                fact2 -= 0.0;

                gradN *= af;

                Dj(0,0) = gradN(0,0) + fact2;
                Dj(0,1) = gradN(0,1);
                Dj(0,2) = af*dNvdx[gp][jj];
                Dj(1,0) = gradN(1,0);
                Dj(1,1) = gradN(1,1) + fact2;
                Dj(1,2) = af*dNvdy[gp][jj];

                // SUPG
                Klocal(TI, TJ)     += Da*Dj(0,0);
                Klocal(TI, TJp1)   += Da*Dj(0,1);
                Klocal(TI, TJp2)   += Da*Dj(0,2);

                Klocal(TIp1, TJ)   += Da*Dj(1,0);
                Klocal(TIp1, TJp1) += Da*Dj(1,1);
                Klocal(TIp1, TJp2) += Da*Dj(1,2);

                Klocal(TI,   TJ)   += ( (tau[0]*af) * b1 * rStab(0) * Nv[gp][jj] );
                Klocal(TI,   TJp1) += ( (tau[0]*af) * b2 * rStab(0) * Nv[gp][jj] );
                Klocal(TIp1, TJ)   += ( (tau[0]*af) * b1 * rStab(1) * Nv[gp][jj] );
                Klocal(TIp1, TJp1) += ( (tau[0]*af) * b2 * rStab(1) * Nv[gp][jj] );

                // PSPG
                Klocal(TIp2, TJ)   += (b1*Dj(0,0) + b2*Dj(1,0))*tau[1];
                Klocal(TIp2, TJp1) += (b1*Dj(0,1) + b2*Dj(1,1))*tau[1];
                Klocal(TIp2, TJp2) += (b1*Dj(0,2) + b2*Dj(1,2))*tau[1];

                // LSIC stabilisation

                //fact2 = rho*af*tau[2];
                //Klocal(TI,   TJ)   += (b1*dN_dx[jj])*fact2;
                //Klocal(TI,   TJp1) += (b1*dN_dy[jj])*fact2;

                //Klocal(TIp1, TJ)   += (b2*dN_dx[jj])*fact2;
                //Klocal(TIp1, TJp1) += (b2*dN_dy[jj])*fact2;

                if(axsy)
                {
                  // diffusion term
                  Klocal(TI, TJ)     += (mu*b4*Nv[gp][jj]/rad/rad);
                  Klocal(TI, TJp2)   -= (b4*Nv[gp][jj]/rad);

                  // continuity equation
                  Klocal(TIp2, TJ)   -= (b4*Nv[gp][jj]/rad);
                }
            }

            Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) );
            Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) );
            Flocal(TIp2) -= (b4*grad.trace());

            // SUPG stabilisation terms
            Flocal(TI)   -= Da*rStab(0);
            Flocal(TIp1) -= Da*rStab(1);

            // PSPG stabilisation terms
            Flocal(TIp2) -= (tau[1]*(b1*rStab(0)+b2*rStab(1)));

            // LSIC stabilisation terms
            //fact2 = tau[2]*rho*grad.trace();

            //Flocal(TI)   -= b1*fact2;
            //Flocal(TIp1) -= b2*fact2;

            if(axsy)
            {
                Flocal(TI)   += (-b4*(mu*vel(0)/rad/rad));
                Flocal(TI)   += (b4*pres/rad);
                Flocal(TIp2) += (b4*vel(0)/rad);
            }
        }
    } //gp

    //printMatrix(Klocal);

    return 0;
}
*/




/*
int LagrangeElem2DNavierStokesQuad4Node::calcStiffnessAndResidual(vector<vector<double> >& node_coords, double* elemData, double* timeData, VectorXd& solnPrev, VectorXd& solnPrev2, VectorXd& solnCur, VectorXd& solnDotCur, MatrixXd& Klocal, VectorXd& Flocal, double timeCur)
{
    // stabilised formulation
    // semi-implicit scheme - type A

    int ii, jj, gp, TI, TIp1, TIp2, count, TJ, TJp1, TJp2;

    double  dvol, b1, b2, b3, b4, b5, b6, b7, b8, xx, yy;
    double  pres, Da, Db, rad, urdr, urdr2, tau[3], CI=4.0;
    double  fact, fact1, fact2, param[2];

    double  xNode[4], yNode[4];
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = node_coords[nodeNums[ii]][0];
      yNode[ii] = node_coords[nodeNums[ii]][1];
    }

    VectorXd  res(3), res2(2), dp(2), Du(2), vel(2), velPrev(2), velExtp(2), velDot(2), force(2), gradTvel(2), rStab(3);
    MatrixXd  Dj(2, 3), grad(2,2), gradN(2,2), stress(2,2), gradPrev(2,2);
    Dj.setZero();
    VectorXd  velTemp(3);


    MatrixXd  matG(3,3);
    matG.setZero();
    fact = 4.0/charlen/charlen;
    matG(0,0) = fact;
    matG(1,1) = fact;

    double  rho = elemData[0];
    double  mu  = elemData[1];

    double  am = timeData[1];
    double  af = timeData[2];
    double  acceFact = timeData[8];
    double  muTaf = mu*af;
    double  dt = timeData[5]/af;

    //KimMoinFlowUnsteadyNavierStokes  analy(rho, mu, 1.0);

    bool axsy = false;
    //axsy = ((int)elemData[2] == 1);

    double  tCur  = timeCur;
    double  tPrev = tCur - dt;


    Klocal.setZero();
    Flocal.setZero();

    for(gp=0; gp<nGP; gp++)
    {
        // compute the gradient of velocity
        xx = 0.0; yy = 0.0;
        vel[0] = vel[1] = 0.0;
        velPrev[0] = velPrev[1] = 0.0;
        velExtp[0] = velExtp[1] = 0.0;
        velDot[0] = velDot[1] = 0.0;
        grad.setZero();
        gradPrev.setZero();
        Du.setZero();
        pres = 0.0;
        dp.setZero();

        for(ii=0; ii<npElem; ii++)
        {
            xx += xNode[ii]*Nv[gp][ii];
            yy += yNode[ii]*Nv[gp][ii];

            TI   = nodeNums[ii]*ndof;
            TIp1 = TI+1;
            TIp2 = TI+2;

            b1 = solnPrev[TI];
            b2 = solnPrev[TIp1];

            velPrev[0] += b1*Nv[gp][ii];
            velPrev[1] += b2*Nv[gp][ii];

            gradPrev(0,0) += b1*dNvdx[gp][ii];
            gradPrev(0,1) += b1*dNvdy[gp][ii];
            gradPrev(1,0) += b2*dNvdx[gp][ii];
            gradPrev(1,1) += b2*dNvdy[gp][ii];

            velExtp[0] += (2.0*solnPrev[TI]-solnPrev2[TI])*Nv[gp][ii];
            velExtp[1] += (2.0*solnPrev[TIp1]-solnPrev2[TIp1])*Nv[gp][ii];;

            b1 = solnCur[TI];
            b2 = solnCur[TIp1];
            b3 = solnCur[TIp2];

            vel[0]     += b1*Nv[gp][ii];
            vel[1]     += b2*Nv[gp][ii];
            pres       += b3*Nv[gp][ii];

            grad(0,0)  += b1*dNvdx[gp][ii];
            grad(0,1)  += b1*dNvdy[gp][ii];
            grad(1,0)  += b2*dNvdx[gp][ii];
            grad(1,1)  += b2*dNvdy[gp][ii];

            velDot[0] += solnDotCur[TI]*Nv[gp][ii];
            velDot[1] += solnDotCur[TIp1]*Nv[gp][ii];

            dp(0)     += b3*dNvdx[gp][ii];
            dp(1)     += b3*dNvdy[gp][ii];
        }

        // this is pseudo-stress
        stress = mu*grad;
        stress(0,0) -= pres;
        stress(1,1) -= pres;

        //cout << xx << '\t' << yy << endl;
        force.setZero();
        //force(0) = af*analy.computeForce(0, geom[0], geom[1], 0.0, tCur)+(1.0-af)*analy.computeForce(0, geom[0], geom[1], 0.0, tPrev);
        //force(1) = af*analy.computeForce(1, geom[0], geom[1], 0.0, tCur)+(1.0-af)*analy.computeForce(1, geom[0], geom[1], 0.0, tPrev);
        //cout << force(0) << '\t' << force(1) << endl;

        velExtp = af*velExtp + (1.0-af)*velPrev;

        gradTvel = grad*velExtp;

        res2(0) = rho*(velDot(0) + gradTvel(0) - force(0)) ;
        res2(1) = rho*(velDot(1) + gradTvel(1) - force(1)) ;

        rStab(0) = res2(0) - mu*Du(0) + dp(0) ;
        rStab(1) = res2(1) - mu*Du(1) + dp(1) ;

        dvol = elemVolGP[gp];
        if(axsy)
        {
            rad = xx;

            urdr  = vel(0)/rad;
            urdr2 = urdr/rad;
            dvol *= (2.0*PI*rad);

            rStab(0) -= mu*(grad(0,0)/rad - urdr2 );
            rStab(1) -= mu*(grad(1,0)/rad );
        }

        velTemp(0) = velPrev(0);
        velTemp(1) = velPrev(1);
        velTemp(2) = 0.0;

        //evaluateStabParams_algo1(&velTemp(0), h, rho, mu, dt,  beta, tau);
        //evaluateStabParams_algo2(&velTemp(0), h, rho, mu, dt,  beta, tau);
        evaluateStabParams_algo3(velTemp, matG, dt, rho, mu, CI, tau);

        tau[0] *= elemData[8];                             // SUPG
        tau[1] *= elemData[9];                             // PSPG
        tau[2] *= elemData[10];                            // LSIC


        for(ii=0;ii<nlbf;ii++)
        {
            TI   = ndof*ii;
            TIp1 = TI+1;
            TIp2 = TI+2;

            b1 = dNvdx[gp][ii]*dvol;
            b2 = dNvdy[gp][ii]*dvol;
            b4 = Nv[gp][ii]*dvol;

            b5 = muTaf*b1;
            b6 = muTaf*b2;
            b8 = af*b4;

            Da = rho*(velPrev(0)*b1 + velPrev(1)*b2)*tau[0];

            for(jj=0;jj<nlbf;jj++)
            {
                TJ   = ndof*jj;
                TJp1 = TJ+1;
                TJp2 = TJ+2;

                fact2 = rho*acceFact*Nv[gp][jj];

                // time acceleration term
                fact = b4*fact2 ;

                // diffusion term
                fact += b5*dNvdx[gp][jj]+b6*dNvdy[gp][jj];

                Klocal(TI,   TJ)   += fact;
                Klocal(TIp1, TJp1) += fact;

                // convection term

                Db = rho*(velExtp(0)*dNvdx[gp][jj] + velExtp(1)*dNvdy[gp][jj]);

                Klocal(TI,   TJ)   += (b8*Db);
                Klocal(TIp1, TJp1) += (b8*Db);

                // pressure term
                Klocal(TI,   TJp2) -= (af*b1*Nv[gp][jj]);
                Klocal(TIp1, TJp2) -= (af*b2*Nv[gp][jj]);

                // continuity equation
                Klocal(TIp2, TJ)   += (b8*dNvdx[gp][jj]);
                Klocal(TIp2, TJp1) += (b8*dNvdy[gp][jj]);

                // SUPG and PSPG stabilisation terms
                //fact2 -= mu*d2Nv[gp][jj];
                fact2 -= 0.0;

                Db *= af;

                Dj(0,0) = Db + fact2;
                Dj(0,1) = 0.0;
                Dj(0,2) = af*dNvdx[gp][jj];
                Dj(1,0) = 0.0;
                Dj(1,1) = Db + fact2;
                Dj(1,2) = af*dNvdy[gp][jj];

                // SUPG
                Klocal(TI, TJ)     += Da*Dj(0,0);
                Klocal(TI, TJp1)   += Da*Dj(0,1);
                Klocal(TI, TJp2)   += Da*Dj(0,2);

                Klocal(TIp1, TJ)   += Da*Dj(1,0);
                Klocal(TIp1, TJp1) += Da*Dj(1,1);
                Klocal(TIp1, TJp2) += Da*Dj(1,2);

                // PSPG
                Klocal(TIp2, TJ)   += (b1*Dj(0,0) + b2*Dj(1,0))*tau[1];
                Klocal(TIp2, TJp1) += (b1*Dj(0,1) + b2*Dj(1,1))*tau[1];
                Klocal(TIp2, TJp2) += (b1*Dj(0,2) + b2*Dj(1,2))*tau[1];

                // LSIC stabilisation
                //fact2 = rho*af*tau[2];
                //Klocal(TI,   TJ)   += (b1*dN_dx[jj])*fact2;
                //Klocal(TI,   TJp1) += (b1*dN_dy[jj])*fact2;
                //Klocal(TIp1, TJ)   += (b2*dN_dx[jj])*fact2;
                //Klocal(TIp1, TJp1) += (b2*dN_dy[jj])*fact2;

                if(axsy)
                {
                  // diffusion term
                  Klocal(TI, TJ)     += (mu*b4*Nv[gp][jj]/rad/rad);
                  Klocal(TI, TJp2)   -= (b4*Nv[gp][jj]/rad);

                  // continuity equation
                  Klocal(TIp2, TJ)   -= (b4*Nv[gp][jj]/rad);
                }
            }

            Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) );
            Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) );
            Flocal(TIp2) -= (b4*grad.trace());

            // SUPG stabilisation terms
            Flocal(TI)   -= Da*rStab(0);
            Flocal(TIp1) -= Da*rStab(1);

            // PSPG stabilisation terms
            Flocal(TIp2) -= (tau[1]*(b1*rStab(0)+b2*rStab(1)));

            // LSIC stabilisation terms
            //fact2 = tau[2]*rho*grad.trace();
            //Flocal(TI)   -= b1*fact2;
            //Flocal(TIp1) -= b2*fact2;

            if(axsy)
            {
                Flocal(TI)   += (-b4*(mu*vel(0)/rad/rad));
                Flocal(TI)   += (b4*pres/rad);
                Flocal(TIp2) += (b4*vel(0)/rad);
            }
        }
    } //gp

    return 0;
}
*/




//
int LagrangeElem2DNavierStokesQuad4Node::calcStiffnessAndResidual(vector<vector<double> >& node_coords, double* elemData, MatrixXdRM& Klocal, VectorXd& Flocal, double timeCur)
{
    // stabilised formulation
    // semi-implicit scheme - type B

    int ii, jj, gp, TI, TIp1, TIp2, count, TJ, TJp1, TJp2;

    double  dvol, b1, b2, b3, b4, b5, b6, b7, b8, xx, yy;
    double  pres, Da, Db, rad, urdr, urdr2, tau[3], CI=4.0;
    double  fact, fact1, fact2, param[2];

    double  xNode[4], yNode[4];
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = node_coords[nodeNums[ii]][0];
      yNode[ii] = node_coords[nodeNums[ii]][1];
    }

    VectorXd  res(3), res2(2), dp(2), Du(2), vel(2), velPrev(2), velExtp(2), velDot(2), force(2), gradTvel(2), rStab(3);
    MatrixXd  Dj(2, 3), grad(2,2), gradN(2,2), stress(2,2), gradPrev(2,2);
    Dj.setZero();
    VectorXd  velTemp(3);


    MatrixXd  matG(3,3);
    matG.setZero();
    fact = 4.0/charlen/charlen;
    matG(0,0) = fact;
    matG(1,1) = fact;

    double  rho = elemData[0];
    double  mu  = elemData[1];

    double  am = SolnData->td[1];
    double  af = SolnData->td[2];
    double  acceFact = SolnData->td[8];
    double  muTaf = mu*af;
    double  dt = SolnData->td[5]/af;

    //KimMoinFlowUnsteadyNavierStokes  analy(rho, mu, 1.0);

    bool axsy = false;
    //axsy = ((int)elemData[2] == 1);

    double  tCur  = timeCur;
    double  tPrev = tCur - dt;


    Klocal.setZero();
    Flocal.setZero();

    for(gp=0; gp<nGP; gp++)
    {
        // compute the gradient of velocity
        xx = 0.0; yy = 0.0;
        vel[0] = vel[1] = 0.0;
        velPrev[0] = velPrev[1] = 0.0;
        velExtp[0] = velExtp[1] = 0.0;
        velDot[0] = velDot[1] = 0.0;
        grad.setZero();
        gradPrev.setZero();
        Du.setZero();
        pres = 0.0;
        dp.setZero();

        for(ii=0; ii<npElem; ii++)
        {
            xx += xNode[ii]*Nv[gp][ii];
            yy += yNode[ii]*Nv[gp][ii];

            TI   = nodeNums[ii]*ndof;
            TIp1 = TI+1;
            TIp2 = TI+2;

            b1 = SolnData->solnPrev[TI];
            b2 = SolnData->solnPrev[TIp1];

            velPrev[0] += b1*Nv[gp][ii];
            velPrev[1] += b2*Nv[gp][ii];

            gradPrev(0,0) += b1*dNvdx[gp][ii];
            gradPrev(0,1) += b1*dNvdy[gp][ii];
            gradPrev(1,0) += b2*dNvdx[gp][ii];
            gradPrev(1,1) += b2*dNvdy[gp][ii];

            b1 = SolnData->solnCur[TI];
            b2 = SolnData->solnCur[TIp1];
            b3 = SolnData->solnCur[TIp2];

            vel[0]     += b1*Nv[gp][ii];
            vel[1]     += b2*Nv[gp][ii];
            pres       += b3*Nv[gp][ii];

            grad(0,0)  += b1*dNvdx[gp][ii];
            grad(0,1)  += b1*dNvdy[gp][ii];
            grad(1,0)  += b2*dNvdx[gp][ii];
            grad(1,1)  += b2*dNvdy[gp][ii];

            velDot[0] += SolnData->solnDotCur[TI]*Nv[gp][ii];
            velDot[1] += SolnData->solnDotCur[TIp1]*Nv[gp][ii];

            dp(0)     += b3*dNvdx[gp][ii];
            dp(1)     += b3*dNvdy[gp][ii];
        }

        // this is pseudo-stress
        stress = mu*grad;
        stress(0,0) -= pres;
        stress(1,1) -= pres;

        //cout << xx << '\t' << yy << endl;
        force.setZero();
        //force(0) = af*analy.computeForce(0, geom[0], geom[1], 0.0, tCur)+(1.0-af)*analy.computeForce(0, geom[0], geom[1], 0.0, tPrev);
        //force(1) = af*analy.computeForce(1, geom[0], geom[1], 0.0, tCur)+(1.0-af)*analy.computeForce(1, geom[0], geom[1], 0.0, tPrev);
        //cout << force(0) << '\t' << force(1) << endl;

        gradTvel = gradPrev*vel + grad*velPrev - gradPrev*velPrev;

        res2(0) = rho*(velDot(0) + gradTvel(0) - force(0)) ;
        res2(1) = rho*(velDot(1) + gradTvel(1) - force(1)) ;

        rStab(0) = res2(0) - mu*Du(0) + dp(0) ;
        rStab(1) = res2(1) - mu*Du(1) + dp(1) ;

        dvol = elemVolGP[gp];
        if(axsy)
        {
            rad = xx;

            urdr  = vel(0)/rad;
            urdr2 = urdr/rad;
            dvol *= (2.0*PI*rad);

            rStab(0) -= mu*(grad(0,0)/rad - urdr2 );
            rStab(1) -= mu*(grad(1,0)/rad );
        }

        velTemp(0) = velPrev(0);
        velTemp(1) = velPrev(1);
        velTemp(2) = 0.0;

        //evaluateStabParams_algo1(&velTemp(0), h, rho, mu, dt,  beta, tau);
        //evaluateStabParams_algo2(&velTemp(0), h, rho, mu, dt,  beta, tau);
        evaluateStabParams_algo3(velTemp, matG, dt, rho, mu, CI, tau);

        tau[0] *= elemData[8];                             // SUPG
        tau[1] *= elemData[9];                             // PSPG
        tau[2] *= elemData[10];                            // LSIC

        for(ii=0;ii<nlbf;ii++)
        {
            TI   = ndof*ii;
            TIp1 = TI+1;
            TIp2 = TI+2;

            b1 = dNvdx[gp][ii]*dvol;
            b2 = dNvdy[gp][ii]*dvol;
            b4 = Nv[gp][ii]*dvol;

            b5 = muTaf*b1;
            b6 = muTaf*b2;
            b8 = af*b4;

            Da = rho*(velPrev(0)*b1 + velPrev(1)*b2)*tau[0];

            for(jj=0;jj<nlbf;jj++)
            {
                TJ   = ndof*jj;
                TJp1 = TJ+1;
                TJp2 = TJ+2;

                fact2 = rho*acceFact*Nv[gp][jj];

                // time acceleration term
                fact = b4*fact2 ;

                // diffusion term
                fact += ( b5*dNvdx[gp][jj]+b6*dNvdy[gp][jj] );

                Klocal(TI,   TJ)   += fact;
                Klocal(TIp1, TJp1) += fact;

                // convection term

                gradN = gradPrev*(rho*Nv[gp][jj]);

                Db = rho*(velPrev(0)*dNvdx[gp][jj] + velPrev(1)*dNvdy[gp][jj]);

                gradN(0,0) += Db;
                gradN(1,1) += Db;

                Klocal(TI,   TJ)   += (b8*gradN(0,0));
                Klocal(TI,   TJp1) += (b8*gradN(0,1));
                Klocal(TIp1, TJ)   += (b8*gradN(1,0));
                Klocal(TIp1, TJp1) += (b8*gradN(1,1));

                // pressure term
                Klocal(TI,   TJp2) -= (b1*af*Nv[gp][jj]);
                Klocal(TIp1, TJp2) -= (b2*af*Nv[gp][jj]);

                // continuity equation
                Klocal(TIp2, TJ)   += (b8*dNvdx[gp][jj]);
                Klocal(TIp2, TJp1) += (b8*dNvdy[gp][jj]);

                // SUPG and PSPG stabilisation terms

                gradN *= af;

                Dj(0,0) = gradN(0,0) + fact2;
                Dj(0,1) = gradN(0,1);
                Dj(0,2) = af*dNvdx[gp][jj];
                Dj(1,0) = gradN(1,0);
                Dj(1,1) = gradN(1,1) + fact2;
                Dj(1,2) = af*dNvdy[gp][jj];

                if(axsy)
                {
                  Dj(0,0) -= muTaf*(dNvdx[gp][jj]/rad - Nv[gp][jj]/rad/rad);
                  Dj(1,1) -= muTaf*(dNvdx[gp][jj]/rad);
                }

                // SUPG
                Klocal(TI, TJ)     += Da*Dj(0,0);
                Klocal(TI, TJp1)   += Da*Dj(0,1);
                Klocal(TI, TJp2)   += Da*Dj(0,2);

                Klocal(TIp1, TJ)   += Da*Dj(1,0);
                Klocal(TIp1, TJp1) += Da*Dj(1,1);
                Klocal(TIp1, TJp2) += Da*Dj(1,2);

                // PSPG stabilisation
                Klocal(TIp2, TJ)   += (b1*Dj(0,0) + b2*Dj(1,0))*tau[1];
                Klocal(TIp2, TJp1) += (b1*Dj(0,1) + b2*Dj(1,1))*tau[1];
                Klocal(TIp2, TJp2) += (b1*Dj(0,2) + b2*Dj(1,2))*tau[1];

                // LSIC stabilisation

                //fact = af*rho*tau[2];

                //Klocal(TI,   TJ)   += (b1*fact*dN_dx[jj]);
                //Klocal(TI,   TJp1) += (b1*fact*dN_dy[jj]);

                //Klocal(TIp1, TJ)   += (b2*fact*dN_dx[jj]);
                //Klocal(TIp1, TJp1) += (b2*fact*dN_dy[jj]);

                if(axsy)
                {
                  // diffusion term
                  Klocal(TI, TJ)     += (b4 * (mu/rad/rad) * (af*Nv[gp][jj]) );
                  Klocal(TI, TJp2)   -= (b4 * Nv[gp][jj]/rad);

                  // continuity equation
                  Klocal(TIp2, TJ)   += (b4 * af*Nv[gp][jj]/rad);
                }
            }

            Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) );
            Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) );
            Flocal(TIp2) -= (b4*grad.trace());

            // SUPG stabilisation terms
            Flocal(TI)   -= Da*rStab(0);
            Flocal(TIp1) -= Da*rStab(1);

            // PSPG stabilisation terms
            Flocal(TIp2) -= (tau[1]*(b1*rStab(0)+b2*rStab(1)));

            // LSIC stabilisation terms
            //fact2 = tau[2]*rho*grad.trace();
            //Flocal(TI)   -= b1*fact2;
            //Flocal(TIp1) -= b2*fact2;

            if(axsy)
            {
                Flocal(TI)   -= (b4 * (mu/rad/rad) * vel(0) );
                Flocal(TI)   += (b4 * pres/rad);
                Flocal(TIp2) -= (b4 * vel(0)/rad);
            }
        }
    } //gp

    //printMatrix(Klocal); printf("\n\n\n");
    //printVector(Flocal); printf("\n\n\n");

    return 0;
}
//



int  LagrangeElem2DNavierStokesQuad4Node::calcError(int index)
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

