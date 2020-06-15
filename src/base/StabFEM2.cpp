
#include "StabFEM.h"
#include "UtilitiesGeneral.h"
#include "ElementBase.h"
#include <chrono>
#include "KimMoinFlow.h"
#include "metis.h"


int StabFEM::setSolver(int slv, int *parm, bool cIO)
{
    Eigen::setNbThreads(0);

    if(solverPetsc != NULL)
      delete solverPetsc;
    solverPetsc = NULL;


    solverPetsc = (SolverPetsc*) new SolverPetsc;

    computerTimePattern = MPI_Wtime();
    prepareMatrixPattern();
    computerTimePattern = MPI_Wtime() - computerTimePattern;

    if(solverPetsc != NULL)
      solverPetsc->checkIO = cIO;

    return 0;
}




int StabFEM::prepareMatrixPattern()
{
    PetscPrintf(MPI_COMM_WORLD, "\n\n     StabFEM::prepareMatrixPattern()  .... STARTED ...\n\n");

    int  size1, size2, row, col;
    int  tempDOF, domTemp, ind;
    int  r, c, r1, c1, count=0, count1=0, count2=0, ii, jj, ee, dd, ind1, ind2, nsize;
    int  *tt1, *tt2, val1, val2, nnz, nnz_max_row, n1, n2, kk, e1, e2, ll, aa, bb, nn, dof;
    int  side, start1, start2, nr1, nr2;

    PetscInt  *colTemp;
    PetscScalar  *arrayTemp;
    double  tstart, tend;

    int ndomains = n_mpi_procs, subdomain=0;


    /////////////////////////////////////////////////////////////
    //
    // prepare the matrix pattern
    /////////////////////////////////////////////////////////////

    // set sizes of some data arrays
    vector<bool>  vecBoolTempFalse(ndof, false);
	  NodeTypeOld.resize(nNode_global, vecBoolTempFalse);
    NodeTypeNew.resize(nNode_global, vecBoolTempFalse);

    vector<int>  vecIntTempM1(ndof, -1);
    NodeDofArrayOld.resize(nNode_global, vecIntTempM1);
    NodeDofArrayNew.resize(nNode_global, vecIntTempM1);

    // fix the specified Dirichlet BCs
    for(ii=0; ii<nDBC; ++ii)
    {
      NodeTypeOld[DirichletBCs[ii][0]][DirichletBCs[ii][1]] = true;
    }

    ntotdofs_global = 0;
    for(ii=0;ii<nNode_global;++ii)
    {
      for(jj=0;jj<ndof;++jj)
      {
        if(!NodeTypeOld[ii][jj])
        {
          NodeDofArrayOld[ii][jj] = ntotdofs_global++;
        }
      }
    }

    if(this_mpi_proc == 0)
    {
      cout << " Mesh statistics .....\n" << endl;
      cout << " nElem_global     = " << '\t' << nElem_global << endl;
      cout << " nNode_global     = " << '\t' << nNode_global  << endl;
      cout << " npElem           = " << '\t' << npElem << endl;
      cout << " ndof             = " << '\t' << ndof << endl;
      cout << " ntotdofs_local   = " << '\t' << ntotdofs_local  << endl;
      cout << " ntotdofs_global  = " << '\t' << ntotdofs_global << endl;
      cout << " n_mpi_procs      = " << '\t' << n_mpi_procs << endl;
    }

    node_map_get_old.resize(nNode_global, 0);
    node_map_get_new.resize(nNode_global, 0);

    dof_map_get_old.resize(ntotdofs_global, 0);
    dof_map_get_new.resize(ntotdofs_global, 0);

    elem_proc_id.resize(nElem_global);
    node_proc_id.resize(nNode_global);

    fill(elem_proc_id.begin(), elem_proc_id.end(), 0);
    fill(node_proc_id.begin(), node_proc_id.end(), 0);

    if(n_mpi_procs == 1)
    {
        elem_start = 0;
        elem_end   = nElem_global-1;

        nElem_local = nElem_global;
        nNode_local = nNode_global;
        ntotdofs_local  = ntotdofs_global;

        row_start = 0;
        row_end   = ntotdofs_global-1;

        for(ii=0; ii<nNode_global; ii++)
        {
          node_map_get_old[ii] = ii;
          node_map_get_new[ii] = ii;
        }
        for(ii=0; ii<ntotdofs_global; ii++)
        {
          dof_map_get_old[ii] = ii;
          dof_map_get_new[ii] = ii;
        }

        for(ii=0;ii<nNode_global;++ii)
        {
          NodeTypeNew[ii]     = NodeTypeOld[ii];
          NodeDofArrayNew[ii] = NodeDofArrayOld[ii];
        }
    }
    else
    {
        PetscPrintf(MPI_COMM_WORLD, "\n\n Before partitionMesh ... \n\n");
        partitionMesh();
        PetscPrintf(MPI_COMM_WORLD, "\n\n After partitionMesh ... \n\n"); 

        errpetsc = MPI_Barrier(MPI_COMM_WORLD);

        PetscPrintf(MPI_COMM_WORLD, "\n\n Before prepareDataForParallel ... \n\n");
        prepareDataForParallel();
        PetscPrintf(MPI_COMM_WORLD, "\n\n After prepareDataForParallel ... \n\n");

        errpetsc = MPI_Barrier(MPI_COMM_WORLD);
    }

    SolnData.node_map_get_old = node_map_get_old;
    SolnData.node_map_get_new = node_map_get_new;

    errpetsc = MPI_Barrier(MPI_COMM_WORLD);

    for(ee=0; ee<nElem_global; ++ee)
    {
      elems[ee]->prepareElemData(node_coords);

      npElem = elems[ee]->nodeNums.size();

      nsize = ndof*npElem;

      elems[ee]->forAssyVec.resize(nsize);

      for(ii=0; ii<npElem; ++ii)
      {
        ind = ndof*ii;

        kk = elems[ee]->nodeNums[ii];

        for(jj=0;jj<ndof;++jj)
        {
          elems[ee]->forAssyVec[ind+jj] = NodeDofArrayNew[kk][jj];
        }
      }
    }

    assyForSoln.resize(ntotdofs_global);
    ind = 0;
    for(ii=0;ii<nNode_global;++ii)
    {
      for(jj=0;jj<ndof;++jj)
      {
        if(NodeDofArrayNew[ii][jj] != -1)
        {
          assyForSoln[ind++] = ii*ndof+jj;
        }
      }
    }

    PetscPrintf(MPI_COMM_WORLD, "\n\n Preparing matrix pattern DONE \n\n");
    PetscPrintf(MPI_COMM_WORLD, "\n\n Element DOF values initialised \n\n");

    cout << " Total DOF   = " << '\t' << ntotdofs_local << '\t' << ntotdofs_global << endl;

    PetscInt  *diag_nnz, *offdiag_nnz;

    errpetsc  = PetscMalloc1(ntotdofs_local,  &diag_nnz);CHKERRQ(errpetsc);
    errpetsc  = PetscMalloc1(ntotdofs_local,  &offdiag_nnz);CHKERRQ(errpetsc);

    n1 = 500; n2 = 250;
    if( (n1 > ntotdofs_local) || (n2 > ntotdofs_local) )
    {
      n1 = ntotdofs_local;
      n2 = n1;
    }

    for(ii=0; ii<ntotdofs_local; ii++)
    {
      diag_nnz[ii]    = n1;
      offdiag_nnz[ii] = n2;
    }

    PetscPrintf(MPI_COMM_WORLD, "\n\n Initialising petsc solver \n\n");

    // Initialize the petsc solver
    solverPetsc->initialise(ntotdofs_local, ntotdofs_global, diag_nnz, offdiag_nnz);


    //Create parallel matrix, specifying only its global dimensions.
    //When using MatCreate(), the matrix format can be specified at
    //runtime. Also, the parallel partitioning of the matrix is
    //determined by Petsc at runtime.
    //Performance tuning note: For problems of substantial size,
    //preallocation of matrix memory is crucial for attaining good
    //performance. See the matrix chapter of the users manual for details.

    PetscPrintf(MPI_COMM_WORLD, " Initialise the Matrix pattern \n", errpetsc);

    ind = npElem*ndof;
    ind = ind*ind;
    PetscScalar  Klocal[ind];
    for(ii=0; ii<ind; ii++)  Klocal[ii] = 0.0;

    vector<int>  vecIntTemp;
    for(ee=0; ee<nElem_global; ee++)
    {
      if(elem_proc_id[ee] == this_mpi_proc)
      {
        size1 = elems[ee]->forAssyVec.size();
        vecIntTemp = elems[ee]->forAssyVec;

        errpetsc = MatSetValues(solverPetsc->mtx, size1, &vecIntTemp[0], size1, &vecIntTemp[0], Klocal, INSERT_VALUES);
      }
    } //for(ee=0;)

    solverPetsc->currentStatus = PATTERN_OK;

    PetscPrintf(MPI_COMM_WORLD, " StabFEM::prepareMatrixPattern()  .... FINISHED. Took %f  milliseconds \n", (tend-tstart)*1000);

    return 0;
}




int StabFEM::partitionMesh()
{
    PetscPrintf(MPI_COMM_WORLD, "\n     StabFEM::partitionMesh()  .... STARTED ...\n");

    if(this_mpi_proc == 0)
    {
        int  ee, ii, jj, kk, n2;
        int  nparts = n_mpi_procs, subdomain=0;

        /////////////////////////////////////////////////////////////////////////////
        //
        // Partition the mesh. Here, METIS is used.
        // 
        /////////////////////////////////////////////////////////////////////////////

        PetscInt  *eptr, *eind;

        errpetsc  = PetscMalloc1(nElem_global+1,  &eptr);CHKERRQ(errpetsc);
        errpetsc  = PetscMalloc1(nElem_global*npElem,  &eind);CHKERRQ(errpetsc);

        vector<int>  vecTemp2;

        eptr[0] = 0;
        kk = 0;
        for(ee=0; ee<nElem_global; ee++)
        {
            eptr[ee+1] = (ee+1)*npElem;

            //vecTemp2 = elems[ee]->nodeNums ;
            vecTemp2 = elemConn[ee];
            //printVector(vecTemp2);

            for(ii=0; ii<npElem; ii++)
              eind[kk+ii] = vecTemp2[ii] ;

            kk += npElem;
        }

        int  ncommon_nodes;
        if(ndim == 2)
          ncommon_nodes = 2;    // 3-noded tria or 4-noded quad
        else
        {
          if(npElem == 4)       // 4-noded tetra element
            ncommon_nodes = 3;
          else
            ncommon_nodes = 4;  // 8-noded hexa element
        }

        idx_t objval;
        idx_t options[METIS_NOPTIONS];

        METIS_SetDefaultOptions(options);

        // Specifies the partitioning method.
        //options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB;          // Multilevel recursive bisectioning.
        options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;        // Multilevel k-way partitioning.

        //options[METIS_OPTION_NSEPS] = 10;

        //options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;     // Edge-cut minimization
        options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;     // Total communication volume minimization

        options[METIS_OPTION_NUMBERING] = 0;  // C-style numbering is assumed that starts from 0.

        cout << " Executing METIS subroutine " << endl;

        // METIS partition routine
        //int ret = METIS_PartMeshNodal(&nElem_global, &nNode_global, eptr, eind, NULL, NULL, &nparts, NULL, options, &objval, elem_proc_id, node_proc_id);
        int ret = METIS_PartMeshDual(&nElem_global, &nNode_global, eptr, eind, NULL, NULL, &ncommon_nodes, &nparts, NULL, options, &objval, &elem_proc_id[0], &node_proc_id[0]);

        if(ret == METIS_OK)
          cout << " METIS partition routine successful "  << endl;
        else
          cout << " METIS partition routine FAILED "  << endl;

        errpetsc = PetscFree(eptr); CHKERRQ(errpetsc);
        errpetsc = PetscFree(eind); CHKERRQ(errpetsc);

        if( 1 < 0)
        {
          for(ee=0; ee<nNode_global; ee++)
            cout << ee << '\t' << node_proc_id[ee] << endl;
          cout << endl;  cout << endl;  cout << endl;

          for(ee=0; ee<nElem_global; ee++)
            cout << ee << '\t' << elem_proc_id[ee] << endl;
          cout << endl;  cout << endl;  cout << endl;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    errpetsc = MPI_Bcast(&elem_proc_id[0], nElem_global, MPI_INT, 0, MPI_COMM_WORLD); CHKERRQ(errpetsc);
    MPI_Barrier(MPI_COMM_WORLD);

    errpetsc = MPI_Bcast(&node_proc_id[0], nNode_global, MPI_INT, 0, MPI_COMM_WORLD); CHKERRQ(errpetsc);
    MPI_Barrier(MPI_COMM_WORLD);

    for(int ee=0; ee<nElem_global; ee++)
    {
      elems[ee]->setSubdomainId(elem_proc_id[ee]);
    }

    PetscPrintf(MPI_COMM_WORLD, "\n     StabFEM::partitionMesh()  .... FINISHED ...\n");

    return 0;
}


int StabFEM::prepareDataForParallel()
{
      int n_subdomains = n_mpi_procs, subdomain=0, ee, ii, jj, kk, n1, n2, ind;

      nElem_local = std::count(elem_proc_id.begin(), elem_proc_id.end(), this_mpi_proc);
      cout << " nElem_local =  " << nElem_local << '\t' << this_mpi_proc << '\t' << n_mpi_procs << endl;

      nNode_owned = std::count(node_proc_id.begin(), node_proc_id.end(), this_mpi_proc);
      cout << " nNode_owned =  " << nNode_owned << '\t' << this_mpi_proc << '\t' << n_mpi_procs << endl;

      MPI_Barrier(MPI_COMM_WORLD);

      // find nodes local to each processor generate the list of locally owned nodes
      std::vector<int>  nodelist_owned(nNode_owned);

      kk=0;
      for(ii=0; ii<nNode_global; ii++)
      {
        if( node_proc_id[ii] == this_mpi_proc )
        {
          nodelist_owned[kk++] = ii;
        }
      }
      cout << " Locally owned nodes " << '\t' << this_mpi_proc << endl;
      //printVector(nodelist_owned);

      MPI_Barrier(MPI_COMM_WORLD);

      // create the vector (of size n_mpi_procs)
      // consisting of nNode_owned from all the processors in the communication
      vector<int>  nNode_owned_vector(n_mpi_procs), nNode_owned_sum(n_mpi_procs);

      MPI_Allgather(&nNode_owned, 1, MPI_INT, &nNode_owned_vector[0], 1, MPI_INT, MPI_COMM_WORLD);
      //printVector(nNode_owned_vector);

      // compute the numbers of first and last nodes in the local processor
      nNode_owned_sum = nNode_owned_vector;
      for(ii=1; ii<n_mpi_procs; ii++)
      {
        nNode_owned_sum[ii] += nNode_owned_sum[ii-1];
      }
      //printVector(nNode_owned_sum);
      node_start = 0;
      if(this_mpi_proc > 0)
        node_start = nNode_owned_sum[this_mpi_proc-1];
      node_end   = nNode_owned_sum[this_mpi_proc]-1;

      cout << " node_start =  " << node_start << '\t' << node_end << '\t' << this_mpi_proc << endl;

      MPI_Barrier(MPI_COMM_WORLD);

      std::vector<int>  displs(n_mpi_procs);

      displs[0] = 0;
      for(ii=0; ii<n_mpi_procs-1; ii++)
        displs[ii+1] = displs[ii] + nNode_owned_vector[ii];

      // create a global list of nodelist_owned
      // which will serve as a mapping from NEW node numbers to OLD node numbers
      errpetsc = MPI_Allgatherv(&nodelist_owned[0], nNode_owned, MPI_INT, &node_map_get_old[0], &nNode_owned_vector[0], &displs[0], MPI_INT, MPI_COMM_WORLD);

      // create an array for mapping from OLD node numbers to NEW node numbers
      // Also, generate NodeTypeNew array for computing the local and global DOF size
      // as well as creating the element-wise array for element matrix/vector assembly
      for(ii=0; ii<nNode_global; ii++)
      {
        n1 = node_map_get_old[ii];
        node_map_get_new[n1] = ii;

        //for(jj=0; jj<ndof; jj++)
        //{
          //NodeTypeNew[ii][jj] = NodeTypeOld[n1][jj];
        //}
      }

      // update elem<->node connectivity with new node numbers
      for(ee=0; ee<nElem_global; ee++)
      {
          for(ii=0; ii<npElem; ii++)
            elemConn[ee][ii] = node_map_get_new[elemConn[ee][ii]];

          elems[ee]->nodeNums = elemConn[ee];
      }

      // update Dirichlet BC information with new node numbers
      for(ii=0; ii<nDBC; ii++)
      {
          n1 = node_map_get_new[DirichletBCs[ii][0]];
          DirichletBCs[ii][0] = n1;

          NodeTypeNew[n1][DirichletBCs[ii][1]] = true;
      }

      MPI_Barrier(MPI_COMM_WORLD);

      // compute NodeDofArrayNew
      ind = 0;
      for(ii=0; ii<nNode_global; ii++)
      {
        for(jj=0; jj<ndof; jj++)
        {
          if(NodeTypeNew[ii][jj] == false)
          {
            NodeDofArrayNew[ii][jj] = ind++;
          }
        }
      }

      if(ind != ntotdofs_global)
      {
        cerr << "Something wrong with NodeDofArrayNew " << endl;
        cout << "ind = " << ind << '\t' << ntotdofs_global << '\t' << this_mpi_proc << endl;
      }

      MPI_Barrier(MPI_COMM_WORLD);

      // compute first and last row indices of the rows owned by the local processor
      row_start  =  1e9;
      row_end    = -1e9;
      ntotdofs_local = 0;
      for(ii=node_start; ii<=node_end; ii++)
      {
        for(jj=0; jj<ndof; jj++)
        {
          if(NodeTypeNew[ii][jj] == false)
          {
            ind = NodeDofArrayNew[ii][jj];
            row_start  = min(row_start, ind);
            row_end    = max(row_end,   ind);
            ntotdofs_local++;
          }
        }
      }

      cout << "ntotdofs_local = " << ntotdofs_local << '\t' << ntotdofs_global << '\t' << this_mpi_proc << endl;

      cout << "row_start  = " << row_start  << '\t' << row_end     << '\t' << this_mpi_proc << endl;

      // check if the sum of local problem sizes is equal to that of global problem size
      ind=0;
      errpetsc = MPI_Allreduce(&ntotdofs_local, &ind, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      if(ind != ntotdofs_global)
      {
        cerr << "Sum of local problem sizes is not equal to global size" << endl;
        cout << "ind = " << ind << '\t' << ntotdofs_global << '\t' << this_mpi_proc << endl;
      }

  return 0;
}



int  StabFEM::solveFullyImplicit()
{
    PetscPrintf(MPI_COMM_WORLD, " Solving with the Fully-Implicit Scheme \n\n\n");

    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////


    int  stepsCompleted=0;
    int  aa, bb, ee, ii, jj, kk, count, row, col, ind, n1, n2, size1, size2;

    double  norm_rhs, fact, fact1, fact2, timerVal;
    double  timeNow=dt, timeFact=0.0;

    VectorXd  reacVec(nNode_global*ndof);
    VectorXd  TotalForce(3);
    ind = npElem*ndof;
    VectorXd  Flocal(ind);
    MatrixXd  Klocal(ind, ind);
    PetscScalar *arrayTempSoln;
    Vec            vec_SEQ;
    VecScatter     ctx;

    //KimMoinFlowUnsteadyNavierStokes  analy(elemData[0], elemData[1], 1.0);

    computerTimeTimeLoop = MPI_Wtime();

    //setInitialConditions();
    //timeFact = 1.0;
    //Time loop
    for(int tstep=0; tstep<stepsMax; tstep++)
    {
        PetscPrintf(MPI_COMM_WORLD, " ==================================================================== \n");
        PetscPrintf(MPI_COMM_WORLD, " Time = %f \n", timeNow);
        PetscPrintf(MPI_COMM_WORLD, " ==================================================================== \n");

        SolnData.setTimeParam();
        SolnData.timeUpdate();

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
        //if(this_mpi_proc == 0) printVector(SolnData.solnApplied);

        MPI_Barrier(MPI_COMM_WORLD);

        PetscPrintf(MPI_COMM_WORLD, " Iteration loop \n");
        for(int iter=0; iter<10; iter++)
        {
            SolnData.updateIterStep();

            solverPetsc->zeroMtx();

            reacVec.setZero();

            MPI_Barrier(MPI_COMM_WORLD);

            PetscPrintf(MPI_COMM_WORLD, "\n Element loop \n");

            timerVal = MPI_Wtime();

            // loop over elements and compute matrix and residual
            for(ee=0; ee<nElem_global; ee++)
            {
              if(elems[ee]->getSubdomainId() == this_mpi_proc)
              {
                elems[ee]->calcStiffnessAndResidual(node_coords, elemData, Klocal, Flocal, timeNow);

                size1 = elems[ee]->forAssyVec.size();

                //PetscPrintf(MPI_COMM_WORLD, "\n Applying boundary conditions %d \n ", ee);
                if(iter == 0)
                {
                    // apply boundary conditions
                    for(ii=0; ii<size1; ii++)
                    {
                        aa = elems[ee]->forAssyVec[ii];

                        if(aa == -1) // this DOF has a prescribed value
                        {
                            fact = SolnData.solnApplied[elems[ee]->globalDOFnums[ii]];

                            // check if fact is zero. We don't need to
                            // execute the for loop if fact is zero.
                            if( abs(fact) > 1.0e-10)
                            {
                              for(jj=0; jj<size1; jj++)
                              {
                                if( elems[ee]->forAssyVec[jj] != -1 )
                                {
                                    Flocal(jj) -= Klocal(jj, ii) * fact;
                                }
                              }
                            }
                        }
                    }
                } // if(iter == 0)

                //PetscPrintf(MPI_COMM_WORLD, "\n Assembling matrices and vectors \n");
                solverPetsc->assembleMatrixAndVectorSerial(elems[ee]->forAssyVec, Klocal, Flocal);
              } // if(elem_proc_id[ee] == this_proc_id)
            } //Element Loop

            timerVal = MPI_Wtime() - timerVal;
            computerTimeAssembly += timerVal;
            PetscPrintf(MPI_COMM_WORLD, "\n\n Elapsed time for matrix assembly = %f seconds \n\n", timerVal );

            MPI_Barrier(MPI_COMM_WORLD);

            //PetscPrintf(MPI_COMM_WORLD, "\n Adding boundary conditions \n");

            // add boundary conditions
            if(iter == 0)
            {
              for(ii=0; ii<nDBC; ++ii)
              {
                n1 = int(DirichletBCs[ii][0]);
                n2 = int(DirichletBCs[ii][1]);

                jj = n1*ndof+n2;

                SolnData.soln[jj]  += SolnData.solnApplied[jj];
              }
            }

            MPI_Barrier(MPI_COMM_WORLD);

            //for(ii=0; ii<10; ii++)
              //VecSetValue(solverPetsc->rhsVec, 0, 1.0, ADD_VALUES);

            VecAssemblyBegin(solverPetsc->rhsVec);
            VecAssemblyEnd(solverPetsc->rhsVec);

            VecNorm(solverPetsc->rhsVec, NORM_2, &norm_rhs);
            solverPetsc->currentStatus = ASSEMBLY_OK;

            MPI_Barrier(MPI_COMM_WORLD);

            //VecView(solverPetsc->rhsVec, PETSC_VIEWER_STDOUT_WORLD);

            PetscPrintf(MPI_COMM_WORLD, " Iteration = %d  \t RHS norm = %E \n", (iter+1), norm_rhs);

            if(norm_rhs < 1.0e-8)
            {
                PetscPrintf(MPI_COMM_WORLD, " Solution converged below the specified tolerance. \n\n");
                break;
            }
            else
            {
                PetscPrintf(MPI_COMM_WORLD, "Assembly done. Solving the matrix system. \n");

                timerVal = MPI_Wtime();
                if( solverPetsc->factoriseAndSolve() )
                {
                  PetscPrintf(MPI_COMM_WORLD, " PETSc solver not converged. \n\n");
                  return -1;
                }
                timerVal = MPI_Wtime() - timerVal;
                computerTimeSolver += timerVal;
                PetscPrintf(MPI_COMM_WORLD, "\n\n Elapsed time for PETSc solver = %f seconds \n\n", timerVal );

                /////////////////////////////////////////////////////////////////////////////
                // get the solution vector onto all the processors
                /////////////////////////////////////////////////////////////////////////////

                if(n_mpi_procs > 1)
                {
                  VecScatterCreateToAll(solverPetsc->solnVec, &ctx, &vec_SEQ);
                  VecScatterBegin(ctx, solverPetsc->solnVec, vec_SEQ, INSERT_VALUES, SCATTER_FORWARD);
                  VecScatterEnd(ctx, solverPetsc->solnVec, vec_SEQ, INSERT_VALUES, SCATTER_FORWARD);

                  VecGetArray(vec_SEQ, &arrayTempSoln);
                }
                else
                {
                  VecGetArray(solverPetsc->solnVec, &arrayTempSoln);
                }

                // update solution vector
                for(ii=0; ii<ntotdofs_global; ii++)
                {
                  SolnData.soln[assyForSoln[ii]]   +=  arrayTempSoln[ii];
                }

                if(n_mpi_procs > 1)
                {
                  VecRestoreArray(vec_SEQ, &arrayTempSoln);
                  VecScatterDestroy(&ctx);
                  VecDestroy(&vec_SEQ);
                }
                else
                {
                  VecRestoreArray(solverPetsc->solnVec, &arrayTempSoln);
                }
            }
        } //Iteration Loop

        PetscPrintf(MPI_COMM_WORLD, " Postprocessing... \n\n");
        if(this_mpi_proc == 0)
        {
            timerVal = MPI_Wtime();
            postProcess();
            computerTimePostprocess += (MPI_Wtime() - timerVal);

            double TotalForce[3] = {0.0, 0.0, 0.0};
            for(ii=0; ii<nOutputFaceLoads; ++ii)
            {
              TotalForce[0] +=  reacVec[outputEdges[ii][0]*ndof];
              TotalForce[1] +=  reacVec[outputEdges[ii][0]*ndof+1];
            }

            fout_convdata <<  timeNow << '\t' << TotalForce[0] << '\t' << TotalForce[1] << endl;
            cout << endl; cout << endl;
        }

        timeNow = timeNow + dt;

        if(timeNow > timeFinal)
          break;
 
    } //Time loop

    computerTimeTimeLoop = MPI_Wtime() - computerTimeTimeLoop;

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

    double totalError = 0.0, timeNow = 5.0;
    //cout << " index = " << index << endl;
    for(int index=0; index<4; index++)
    {
      totalError = 0.0;
      for(int ee=0; ee<nElem_global; ee++)
      {
        //Compute the element force vector, including residual force
        //totalError += elems[ee]->CalculateError(node_coords, elemData, timeData, soln, solnDot, timeNow, index);
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



