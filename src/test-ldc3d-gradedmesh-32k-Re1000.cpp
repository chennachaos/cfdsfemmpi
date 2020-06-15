/*
! Program for incompressible Navier-Stokes using Stabilised Finite Element Method
!
!
! Author: Dr. Chennakesava Kadapa
! Date  : 06-May-2020
! Place : Swansea, UK
!
!
*/


#include "headersVTK.h"
#include "headersBasic.h"
#include "headersEigen.h"
#include "StabFEM.h"

using namespace std;



int main(int argc, char* argv[])
{
    PetscInitialize(NULL, NULL, "../testdata/ldc3d-gradedmesh-32k-Re1000/petsc_options.dat", NULL);

    string  meshfile = "../testdata/ldc3d-gradedmesh-32k-Re1000/LDC3d-hex-graded-32k-mesh";
    string  controlfile = "../testdata/ldc3d-gradedmesh-32k-Re1000/control-parameters-Re1000-dt0p1.dat";

    StabFEM  stabfem;

    stabfem.readInputData(meshfile);

    stabfem.readControlParameters(controlfile);

    stabfem.prepareInputData();

    int parm[3];
    stabfem.setSolver(1, parm, false);

    stabfem.solveFullyImplicit();

    PetscPrintf(MPI_COMM_WORLD, " Simulation is completed. ... Checking the solution... \n ");

    string  resultfile = "../testdata/ldc3d-gradedmesh-32k-Re1000/LDC3d-gradedmesh-32k-Re1000-ref.dat";

    if(stabfem.checkResult(resultfile) != 0)
    {
      PetscPrintf(MPI_COMM_WORLD, "\n\n\n Result does not match \n\n\n");
      return 1;
    }
    PetscPrintf(MPI_COMM_WORLD, "\n\n\n Result matches \n\n\n");

    PetscPrintf(MPI_COMM_WORLD, "\n\n\n Test is successful \n\n\n ");

    PetscFinalize(); //CHKERRQ(ierr);

    return 0;
}

