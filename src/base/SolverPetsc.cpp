
#include "SolverPetsc.h"
#include "petscmat.h"
#include "petscksp.h"



SolverPetsc::SolverPetsc()
{
  FREED = false;
}


SolverPetsc::~SolverPetsc()
{
  //PetscPrintf(MPI_COMM_WORLD, "SolverPetsc::~SolverPetsc() \n");
  cout << "SolverPetsc::~SolverPetsc() " << endl;
  //free();
/*
  errpetsc = VecDestroy(&soln);
  //cout << " errpetsc = " << errpetsc << endl;
  errpetsc = VecDestroy(&solnPrev);
  //cout << " errpetsc = " << errpetsc << endl;
  errpetsc = VecDestroy(&rhsVec);
  //cout << " errpetsc = " << errpetsc << endl;
  errpetsc = VecDestroy(&reac);
  //cout << " errpetsc = " << errpetsc << endl;
  errpetsc = MatDestroy(&mtx);
  //cout << " errpetsc = " << errpetsc << endl;
  //errpetsc = KSPGetPC(ksp,&pc);
  //cout << " errpetsc = " << errpetsc << endl;
  //errpetsc = PCDestroy(&pc);
  //cout << " errpetsc = " << errpetsc << endl;
  errpetsc = PCReset(pc);
  //cout << " errpetsc = " << errpetsc << endl;
  errpetsc = KSPDestroy(&ksp);
  //errpetsc = KSPReset(ksp);CHKERRQ(errpetsc);
  //cout << " errpetsc = " << errpetsc << endl;
*/
  //PetscPrintf(MPI_COMM_WORLD, "SolverPetsc::~SolverPetsc() \n");
  cout << "SolverPetsc::~SolverPetsc() " << endl;
}



// Initialise Petsc solver
// size_local  - number of local rows/columns in the matrix
// size_global - number of global rows/columns in the matrix
// diag_nnz - number of nonzeros in the diagonal matrix
// offdiag_nnz - number of nonzeros in the off-diagonal matrix
// 
int SolverPetsc::initialise(int size_local, int size_global, int* diag_nnz, int* offdiag_nnz)
{
    nRow = nCol = size_global;

    int dummy = 50;

    // Create PETSc vector
    errpetsc = VecCreate(PETSC_COMM_WORLD, &solnVec);
    CHKERRQ(errpetsc);

    errpetsc = VecSetSizes(solnVec, size_local, size_global);
    CHKERRQ(errpetsc);

    errpetsc = VecSetFromOptions(solnVec);
    CHKERRQ(errpetsc);

    errpetsc = VecDuplicate(solnVec, &rhsVec);
    CHKERRQ(errpetsc);

    errpetsc = VecSetOption(rhsVec, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

    errpetsc = VecDuplicate(solnVec, &reacVec);
    CHKERRQ(errpetsc);

    //PetscPrintf(PETSC_COMM_WORLD, " Creating PETSc matrices \n", errpetsc)

    // Create PETSc matrix
    errpetsc = MatCreate(PETSC_COMM_WORLD, &mtx);
    CHKERRQ(errpetsc);

    errpetsc = MatSetSizes(mtx, size_local, size_local, size_global, size_global);
    CHKERRQ(errpetsc);

    errpetsc = MatSetFromOptions(mtx);
    CHKERRQ(errpetsc);

    errpetsc = MatMPIAIJSetPreallocation(mtx, dummy, diag_nnz, dummy, offdiag_nnz);
    CHKERRQ(errpetsc);

    errpetsc = MatSeqAIJSetPreallocation(mtx, dummy, diag_nnz);
    CHKERRQ(errpetsc);


    errpetsc = MatSetOption(mtx, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    CHKERRQ(errpetsc);

    errpetsc = MatSetOption(mtx, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
    CHKERRQ(errpetsc);

    errpetsc = MatSetOption(mtx, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
    CHKERRQ(errpetsc);

    PetscPrintf(MPI_COMM_WORLD, "\n\n Creating KSP context ... \n\n");

    // Create the KSP context
    errpetsc = KSPCreate(PETSC_COMM_WORLD, &ksp);
    CHKERRQ(errpetsc);

    // Set the operators for the KSP context
    errpetsc = KSPSetOperators(ksp, mtx, mtx);
    CHKERRQ(errpetsc);

    //  Set whether to use non-zero initial guess or not
    //  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    //  KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);

    PetscPrintf(MPI_COMM_WORLD, "\n\n Setting KSP context from input file ... \n\n");

    // Set KSP options from the input file
    // This is convenient as it allows to choose different options
    // from the input files instead of recompiling the code
    errpetsc = KSPSetFromOptions(ksp);    CHKERRQ(errpetsc);

    PetscPrintf(MPI_COMM_WORLD, "\n\n Creating PC context ... \n\n");

    // Get the PC context
    errpetsc = KSPGetPC(ksp, &pc);    CHKERRQ(errpetsc);

    PetscPrintf(MPI_COMM_WORLD, "\n\n Setting PC context from input file ... \n\n");

    // Set PC options from the input file
    errpetsc = PCSetFromOptions(pc);    CHKERRQ(errpetsc);

    currentStatus = SOLVER_EMPTY;

    return 0;
}


int SolverPetsc::setSolverAndParameters()
{
    return 0;
}




int SolverPetsc::zeroMtx()
{
  errpetsc = MatAssemblyBegin(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(errpetsc);
  errpetsc = MatAssemblyEnd(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(errpetsc);

  MatZeroEntries(mtx);

  VecAssemblyBegin(rhsVec);
  VecAssemblyEnd(rhsVec);

  VecZeroEntries(rhsVec);

  VecAssemblyBegin(reacVec);
  VecAssemblyEnd(reacVec);

  VecZeroEntries(reacVec);

  return 0;
}



int SolverPetsc::free()
{
  //if(FREED) return 0;

  PetscPrintf(MPI_COMM_WORLD, "SolverPetsc::free() \n");

  errpetsc = VecDestroy(&solnVec);CHKERRQ(errpetsc);
  cout << " errpetsc = " << errpetsc << endl;
  //errpetsc = VecDestroy(&solnVecPrev);CHKERRQ(errpetsc);
  cout << " errpetsc = " << errpetsc << endl;
  errpetsc = VecDestroy(&rhsVec);CHKERRQ(errpetsc);
  cout << " errpetsc = " << errpetsc << endl;
  errpetsc = VecDestroy(&reacVec);CHKERRQ(errpetsc);
  cout << " errpetsc = " << errpetsc << endl;
  errpetsc = MatDestroy(&mtx);CHKERRQ(errpetsc);
  cout << " errpetsc = " << errpetsc << endl;

  //errpetsc = KSPGetPC(ksp,&pc);CHKERRQ(errpetsc);
  //cout << " errpetsc = " << errpetsc << endl;
  //errpetsc = PCDestroy(&pc);CHKERRQ(errpetsc);
  //cout << " errpetsc = " << errpetsc << endl;
  errpetsc = PCReset(pc);CHKERRQ(errpetsc);
  //cout << " errpetsc = " << errpetsc << endl;
  errpetsc = KSPDestroy(&ksp);CHKERRQ(errpetsc);
  cout << " errpetsc = " << errpetsc << endl;
  errpetsc = KSPReset(ksp);CHKERRQ(errpetsc);
  cout << " errpetsc = " << errpetsc << endl;

  FREED = true;

  PetscPrintf(MPI_COMM_WORLD, "SolverPetsc::free() \n");

  return 0;
}




int SolverPetsc::printMatrix(int dig, int dig2, bool gfrmt, int indent, bool interactive)
{
  errpetsc = MatAssemblyBegin(mtx, MAT_FINAL_ASSEMBLY);//CHKERRQ(errpetsc);
  errpetsc = MatAssemblyEnd(mtx, MAT_FINAL_ASSEMBLY);//CHKERRQ(errpetsc);
 
  //PetscViewerCreate(PETSC_COMM_WORLD, &viewer_matx);
  //PetscViewerDrawOpen();
  //PetscViewerSetFormat(viewer_matx, PETSC_VIEWER_ASCII_MATLAB);
  //MatView(mtx, PETSC_VIEWER_STDOUT_WORLD);
  //MatView(A,viewer);

  return 0;
}


double SolverPetsc::getMatrixCoefficient(int row, int col)
{ 
  //return mtx(row,col,true);
  return 0.0;
}



int SolverPetsc::factorise()
{
  if (currentStatus != ASSEMBLY_OK)
  {
    cerr << " SolverPetsc::factorise ... assemble matrix first!" << endl;
    return -1;
  }

  if (checkIO)
  {
    // search for "nan" entries in matrix coefficients

    //if (prgNAN(mtx.x.x,NE)) prgError(1,fct,"nan matrix coefficient!");
  }

  currentStatus = FACTORISE_OK;

  return 0;
}


int SolverPetsc::solve()
{
  //PetscPrintf(MPI_COMM_WORLD, "  SolverPetsc::solve() ...  \n\n");

  //time_t tstart, tend; 

  if (currentStatus != FACTORISE_OK)
  {
    cerr << " SolverPetsc::solve ... factorise matrix first!" << endl;
    return -1;
  }

  errpetsc = MatAssemblyBegin(mtx,MAT_FINAL_ASSEMBLY); CHKERRQ(errpetsc);
  errpetsc = MatAssemblyEnd(mtx,MAT_FINAL_ASSEMBLY); CHKERRQ(errpetsc);

  //PetscPrintf(MPI_COMM_WORLD, "  SolverPetsc::solve() ... Matrix Assembly ...  \n\n");

  //PetscViewerCreate(PETSC_COMM_WORLD, &viewer_matx);
  //PetscViewerDrawOpen();
  //PetscViewerSetFormat(viewer_matx, PETSC_VIEWER_ASCII_MATLAB);

  //MatView(mtx,PETSC_VIEWER_STDOUT_WORLD);

  errpetsc = VecAssemblyBegin(rhsVec); CHKERRQ(errpetsc);
  errpetsc = VecAssemblyEnd(rhsVec); CHKERRQ(errpetsc);

  //PetscPrintf(MPI_COMM_WORLD, "  SolverPetsc::solve() ... rhsVec Assembly ...  \n\n");

  errpetsc = VecAssemblyBegin(solnVec); CHKERRQ(errpetsc);
  errpetsc = VecAssemblyEnd(solnVec); CHKERRQ(errpetsc);

  //PetscPrintf(MPI_COMM_WORLD, "  SolverPetsc::solve() ... soln Assembly ...  \n\n");

  errpetsc = VecZeroEntries(solnVec); CHKERRQ(errpetsc);

  //PetscPrintf(MPI_COMM_WORLD, "  SolverPetsc::solve() ... vec zero ...  \n\n");

  //VecView(rhsVec, PETSC_VIEWER_STDOUT_WORLD);

  //tstart = time(0);

  errpetsc = KSPSolve(ksp, rhsVec, solnVec); CHKERRQ(errpetsc);

  //PetscPrintf(MPI_COMM_WORLD, "  SolverPetsc::solve() ... KSP solve ...  \n\n");

  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);

  PetscInt its;

  errpetsc = KSPGetIterationNumber(ksp, &its); CHKERRQ(errpetsc);

  if(reason<0)
  {
    PetscPrintf(MPI_COMM_WORLD, "\n Divergence... %d iterations. \n\n", its);
    cout <<  reason << endl;
    exit(1);
    return -1;
  }
  else
  {
    PetscPrintf(MPI_COMM_WORLD, "\n Convergence in %d iterations.\n\n", its);
  }

  //errpetsc = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(errpetsc);
  //VecView(soln, PETSC_VIEWER_STDOUT_WORLD);

  //tend = time(0); 

  //cout << "SolverPetsc::solve()  took "<< difftime(tend, tstart) <<" second(s)."<< endl;

  return 0;
}






int SolverPetsc::factoriseAndSolve()
{
  char fct[] = "SolverPetsc::factoriseAndSolve";

  if(currentStatus != ASSEMBLY_OK)
  {
    cerr << " SolverPetsc::factoriseAndSolve ... assemble matrix first!" << endl;
    return -1;
  }

  factorise();
  solve();

  return 0;
}





int SolverPetsc::assembleMatrixAndVectorSerial(vector<int>& forAssyElem, MatrixXd& Klocal, VectorXd& Flocal)
{
  int  size1 = forAssyElem.size();

  MatrixXdRM Klocal2 = Klocal;

  VecSetValues(rhsVec, size1, &forAssyElem[0], &Flocal[0], ADD_VALUES);
  MatSetValues(mtx,    size1, &forAssyElem[0], size1, &forAssyElem[0], &Klocal2(0,0), ADD_VALUES);

  return 0;
}






int SolverPetsc::assembleMatrixAndVectorParallel(vector<int>& forAssyElem, vector<int>& dof_map, MatrixXdRM& Klocal, VectorXd& Flocal)
{
  int ii, jj, r, kk=0;

  int  size1 = forAssyElem.size();
  int  size2 = size1*size1;

  PetscInt  row1[size1];
  PetscScalar  array[size2];

  for(ii=0;ii<size1;ii++)
  {
    r = dof_map[forAssyElem[ii]];

    VecSetValue(rhsVec, r, Flocal(ii), ADD_VALUES);

    row1[ii] = r;

    for(jj=0;jj<size1;jj++)
    {
      array[kk++] = Klocal(ii,jj);
    }
  }

  MatSetValues(mtx, size1, row1, size1, row1, array, ADD_VALUES);

  return 0;
}

