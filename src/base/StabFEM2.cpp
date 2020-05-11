
#include "StabFEM.h"
#include "UtilitiesGeneral.h"
#include "ElementBase.h"
#include <chrono>
#include "KimMoinFlow.h"



int StabFEM::setSolver(int slv, int *parm, bool cIO)
{
    Eigen::setNbThreads(0);

    if(solverPetsc != NULL)
      delete solverPetsc;
    solverPetsc = NULL;


    solverPetsc = (SolverPetsc*) new SolverPetsc;

    prepareMatrixPattern();

    //solverPetsc->setSolverAndParameters();

    if(solverPetsc != NULL)
      solverPetsc->checkIO = cIO;

    return;
}




int StabFEM::prepareMatrixPattern()
{
    cout <<  "\n     StabFEM::prepareMatrixPattern()  .... STARTED ...\n" <<  endl;

    int  size1, size2, row, col;
    int  tempDOF, domTemp, npElem, ind;
    int  r, c, r1, c1, count=0, count1=0, count2=0, ii, jj, ee, dd, ind1, ind2;
    int  *tt1, *tt2, val1, val2, nnz, nnz_max_row, n1, n2, kk, e1, e2, ll, aa, bb, nn, dof;
    int  side, start1, start2, nr1, nr2;

    PetscInt  *colTemp;
    PetscScalar  *arrayTemp;
    double  tstart, tend;

    /////////////////////////////////////////////////////////////
    //
    // prepare the matrix pattern
    /////////////////////////////////////////////////////////////

    vector<vector<int> >  ID;
    vector<vector<bool> >  NodeType;

    NodeType.resize(nNode);
    ID.resize(nNode);

    for(ii=0;ii<nNode;++ii)
    {
      NodeType[ii].resize(ndof);
      ID[ii].resize(ndof);

      for(jj=0;jj<ndof;++jj)
      {
        NodeType[ii][jj] = false;
        ID[ii][jj] = -1;
      }
    }

    // fix the specified Dirichlet BCs
    for(ii=0; ii<nDBC; ++ii)
    {
      NodeType[DirichletBCs[ii][0]][DirichletBCs[ii][1]] = true;
    }

    totalDOF = 0;
    for(ii=0;ii<nNode;++ii)
    {
      for(jj=0;jj<ndof;++jj)
      {
        if(!NodeType[ii][jj])
        {
          ID[ii][jj] = totalDOF++;
          assyForSoln.push_back(ii*ndof+jj);
        }
      }
    }

    cout << " Mesh statistics .....\n" << endl;
    cout << " nElem          = " << '\t' << nElem << endl;
    cout << " nNode          = " << '\t' << nNode  << endl;
    cout << " npElem         = " << '\t' << npElem << endl;
    cout << " ndof           = " << '\t' << ndof << endl;
    cout << " Total DOF      = " << '\t' << totalDOF << endl;


    for(ee=0;ee<nElem;++ee)
    {
      npElem = elems[ee]->nodeNums.size();

      ind = ndof*npElem;

      elems[ee]->forAssyVec.resize(ind);

      for(ii=0; ii<npElem; ++ii)
      {
        ind = ndof*ii;

        kk = elems[ee]->nodeNums[ii];

        for(jj=0;jj<ndof;++jj)
        {
          elems[ee]->forAssyVec[ind+jj] = ID[kk][jj];
        }
      }

      //printVector(elems[ee]->forAssyVec);
    }


    bool pp=false;
    //pp=true;
    if(pp)
    {
       printf("   ID array \n\n");
       for(ii=0;ii<nNode;++ii)
       {
          for(jj=0;jj<ndof;++jj)
            cout << '\t' << ID[ii][jj];
          cout << endl;
       }
       printf("\n\n\n");

       printf("  assyForSoln array \n\n");
       for(ii=0;ii<totalDOF;++ii)
       {
          cout << assyForSoln[ii] << endl;
       }
       printf("\n\n\n");
    }

    printf("\n element DOF values initialised \n\n");

    // remove data objects

    //ID.clear();
    //forAssyMat.clear();

    printf("\n     StabFEM::prepareMatrixPattern()  .... FINISHED ...\n\n");

    cout <<  " StabFEM::setSolverDataForFullyImplicit() ... STARTED " << endl;

    cout << " Total DOF   = " << '\t' << totalDOF << endl;


    ndofs_local = totalDOF;

    VectorXi nnzVec(totalDOF);

    ind = 50;
    if(totalDOF < 50) ind = totalDOF;

    for(ii=0; ii<totalDOF; ii++)
      nnzVec(ii) = ind;

    cout <<  " Initialising petsc solver " << endl;

    // Initialize the petsc solver
    solverPetsc->initialise(totalDOF, totalDOF, &nnzVec[0], &nnzVec[0]);


    //Create parallel matrix, specifying only its global dimensions.
    //When using MatCreate(), the matrix format can be specified at
    //runtime. Also, the parallel partitioning of the matrix is
    //determined by Petsc at runtime.
    //Performance tuning note: For problems of substantial size,
    //preallocation of matrix memory is crucial for attaining good
    //performance. See the matrix chapter of the users manual for details.

    ind = npElem*ndof;
    ind = ind*ind;
    PetscScalar  Klocal[ind];
    for(ii=0; ii<ind; ii++)  Klocal[ii] = 0.0;

    vector<int>  vecIntTemp;
    for(ee=0; ee<nElem; ee++)
    {
        size1 = elems[ee]->forAssyVec.size();
        vecIntTemp = elems[ee]->forAssyVec;

        errpetsc = MatSetValues(solverPetsc->mtx, size1, &vecIntTemp[0], size1, &vecIntTemp[0], Klocal, INSERT_VALUES);

        //printVector(elems[ee]->forAssyVec);
    } //for(ee=0;)


    solverPetsc->currentStatus = PATTERN_OK;

    PetscPrintf(MPI_COMM_WORLD, " HBSplineCutFEM::prepareMatrixPattern()  .... FINISHED. Took %f  milliseconds \n", (tend-tstart)*1000);

    cout <<  " StabFEM::setSolverDataForFullyImplicit() ... ENDED " << endl;

    return;
}


int  StabFEM::solveFullyImplicit()
{
    cout << " Solving with the Fully-Implicit Scheme " << endl;

    int parm[3];

    setSolver(1, parm, false);

    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////

    soln.setZero();      solnPrev.setZero();      solnPrev2.setZero();      solnCur.setZero();
    solnDot.setZero();   solnDotPrev.setZero();   solnDotCur.setZero();

    int  stepsCompleted=0;
    int  aa, bb, ee, ii, jj, kk, count, row, col, ind, n1, n2, size1, size2;

    double  norm_rhs, fact, fact1, fact2;
    double  timeNow=dt, timeFact=0.0;

    VectorXd  td(100), reacVec(nNode*ndof);
    VectorXd  TotalForce(3);
    ind = npElem*ndof;
    VectorXd  Flocal(ind);
    MatrixXd  Klocal(ind, ind);
    PetscScalar *arrayTempSoln;

    SetTimeParametersFluid(tis, rhoInf, dt, td);

    //KimMoinFlowUnsteadyNavierStokes  analy(elemData[0], elemData[1], 1.0);


    //setInitialConditions();
    //timeFact = 1.0;
    //Time loop
    for(int tstep=0; tstep<stepsMax; tstep++)
    {
        cout << " Time = " << timeNow << endl;

        //
        if(timeNow <= 1.0)
        {
          timeFact = timeNow;
          //timeFact = stepsCompleted/5000.0;
        }
        else
        {
          timeFact = 1.0;
        }
        //
        //timeFact = 1.0;

        assignBoundaryConditions(timeNow, dt, timeFact);

        cout << " Iteration loop " << endl;
        for(int iter=0; iter<10; iter++)
        {
            solnDot    = td[9]*soln + td[10]*solnPrev + td[15]*solnDotPrev;

            solnCur    = td[2]*soln    + (1.0-td[2])*solnPrev;
            solnDotCur = td[1]*solnDot + (1.0-td[1])*solnDotPrev;

            //cout << " Zeroing the solver " << endl;
            solverPetsc->zeroMtx();

            reacVec.setZero();

            //cout << " Element loop " << endl;

            // loop over elements and compute matrix and residual
            for(ee=0; ee<nElem; ee++)
            {
                elems[ee]->calcStiffnessAndResidual(node_coords, elemData, &td[0], solnPrev, solnPrev2, solnCur, solnDotCur, Klocal, Flocal, timeNow);

                size1 = elems[ee]->forAssyVec.size();

                //cout << "Applying boundary conditions " << ee << endl;
                if(iter == 0)
                {
                    // apply boundary conditions
                    for(ii=0; ii<size1; ii++)
                    {
                        aa = elems[ee]->forAssyVec[ii];

                        if(aa == -1) // this DOF has a prescribed value
                        {
                            fact = solnApplied[elems[ee]->globalDOFnums[ii]];

                            for(jj=0; jj<size1; jj++)
                            {
                                if( elems[ee]->forAssyVec[jj] != -1 )
                                {
                                    Flocal(jj) -= Klocal(jj, ii) * fact;
                                }
                            }
                        }
                    }
                } // if(iter == 0)

                //cout << "Assembling matrices and vectors " << endl;

                // assemble matrices and vectors
                for(ii=0; ii<size1; ii++)
                {
                  row = elems[ee]->forAssyVec[ii];

                  reacVec[elems[ee]->globalDOFnums[ii]] += Flocal[ii];

                  if(row != -1)
                  {
                    VecSetValue(solverPetsc->rhsVec, row, Flocal(ii), ADD_VALUES);

                    for(jj=0; jj<size1; jj++)
                    {
                      col = elems[ee]->forAssyVec[jj];

                      if(col != -1)
                      {
                        MatSetValue(solverPetsc->mtx, row, col, Klocal(ii,jj), ADD_VALUES);
                        //solverPetsc->assembleMatrixAndVectorCutFEM(start, start, elems[ee]->forAssyVec, dof_map_old_to_new, Klocal, Flocal);
                      }
                    }
                  }//if(row != -1)
                } //for(ii=0;)
            } //Element Loop

            //cout << "Adding boundary conditions " << endl;

            // add boundary conditions
            if(iter == 0)
            {
              for(ii=0; ii<nDBC; ++ii)
              {
                n1 = DirichletBCs[ii][0];
                n2 = DirichletBCs[ii][1];

                jj = n1*ndof+n2;

                soln[jj]  += solnApplied[jj];
              }
            }

            VecAssemblyBegin(solverPetsc->rhsVec);
            VecAssemblyEnd(solverPetsc->rhsVec);

            VecNorm(solverPetsc->rhsVec, NORM_2, &norm_rhs);
            solverPetsc->currentStatus = ASSEMBLY_OK;

            printf(" RHS norm = %E \n", norm_rhs);

            if(norm_rhs < 1.0e-8)
            {
                cout << " Solution converged below the specified tolerance " << endl;
                break;
            }
            else
            {
                cout << "Solving the matrix system " << endl;

                if( solverPetsc->factoriseAndSolve() )
                {
                  cerr << " PETSc solver not converged " << endl;
                  return -1;
                }

                cout << "Adding the solution increment " << endl;

                VecGetArray(solverPetsc->solnVec, &arrayTempSoln);

                //printVector(rhsVec);
                //printVector(solnVec);

                for(ii=0; ii<totalDOF; ii++)
                {
                  soln[assyForSoln[ii]]   +=  arrayTempSoln[ii];
                }

                VecRestoreArray(solverPetsc->solnVec, &arrayTempSoln);
            }
        } //Iteration Loop


        cout << "Postprocessing " << endl;
        postProcess();

        double TotalForce[2] = {0.0, 0.0};
        for(ii=0; ii<nOutputFaceLoads; ++ii)
        {
          TotalForce[0] +=  reacVec[outputEdges[ii][0]*ndof];
          TotalForce[1] +=  reacVec[outputEdges[ii][0]*ndof+1];
        }

        fout_convdata <<  timeNow << '\t' << TotalForce[0] << '\t' << TotalForce[1] << endl;
        cout << endl; cout << endl;

        solnPrev2    =  solnPrev;
        solnPrev     =  soln;
        solnDotPrev  =  solnDot;


        timeNow = timeNow + dt;

        if(timeNow > timeFinal)
          break;
 
    } //Time loop


    /*
    cout << " Computing errors \n " << endl;
    double totalError = 0.0;
    //cout << " index = " << index << endl;
    for(int index=0; index<4; index++)
    {
      totalError = 0.0;
      for(ee=0; ee<nElem; ee++)
      {
        //Compute the element force vector, including residual force
        totalError += elems[ee]->CalculateError(node_coords, elemData, timeData, velo, veloDot, pres, 5.0, index);
      }

      totalError = sqrt(totalError);

      if(index == 0)
        printf(" \n\n \t L2 Error in X-velocity = %12.6E \n\n " , totalError);
      else if(index == 1)
        printf(" \n\n \t L2 Error in Y-velocity = %12.6E \n\n " , totalError);
      else if(index == 2)
        printf(" \n\n \t L2 Error in pressure   = %12.6E \n\n " , totalError);
      else
        printf(" \n\n \t H1 Error in velocity   = %12.6E \n\n " , totalError);
    }
    */


    return 0;
}




void StabFEM::addExternalForces(double loadFact)
{
    int  nn, dof, ii, ind;
    double specVal=0.0;

    VectorXd  vecTemp, Flocal;
    vecTemp.setZero();

    // specified nodal forces
    for(ii=0;ii<nodeForcesData.size();++ii)
    {
      nn  = (int) (nodeForcesData[ii][0] - 1);
      dof = (int) (nodeForcesData[ii][1] - 1);
      specVal = nodeForcesData[ii][2];

      ind = nn*ndof+dof;

      vecTemp[ind] += specVal;
    }
    //printVector(vecTemp);

    return;
}

void StabFEM::computeElementErrors(int ind)
{
    cout << " Computing errors \n " << endl;

    double totalError = 0.0, timeNow;
    for(int index=0; index<4; index++)
    {
      totalError = 0.0;
      for(int ee=0; ee<nElem; ++ee)
      {
        //totalError += elems[ee]->CalculateError(node_coords, elemData, timeData, solnVTK, veloDot, pres, timeNow, index);
      }

      totalError = sqrt(totalError);

      if(index == 0)
        printf(" \n\n \t L2 Error in X-velocity = %12.6E \n\n " , totalError);
      else if(index == 1)
        printf(" \n\n \t L2 Error in Y-velocity = %12.6E \n\n " , totalError);
      else if(index == 2)
        printf(" \n\n \t L2 Error in pressure   = %12.6E \n\n " , totalError);
      else
        printf(" \n\n \t H1 Error in velocity   = %12.6E \n\n " , totalError);
    }

    return;
}



