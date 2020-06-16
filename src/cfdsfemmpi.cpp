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
    //PetscInitialize(NULL, NULL, petscfile.c_str(), NULL);
    PetscInitialize(NULL, NULL, "../input/petsc_options_mumps.dat", NULL);
    //PetscInitialize(NULL, NULL, "petsc_options_mumps.dat", NULL);
    PetscMemorySetGetMaximumUsage();

    //Set the input file name
    //The file name is specified from the command line
    if(argc == 0)
    {
        cerr << " Error in input data " << endl;
        cerr << " You must enter name of input file" << endl;
        cerr << " Aborting..." << endl;
    }

    string  meshfile    = argv[1];
    string  controlfile = argv[2];
    //string  petscfile   = argv[3];


    PetscPrintf(MPI_COMM_WORLD, "\n    Input mesh file         = %s \n ", meshfile.c_str());
    PetscPrintf(MPI_COMM_WORLD, "\n    Control parameters file = %s \n ", controlfile.c_str());
    //PetscPrintf(MPI_COMM_WORLD, "\n    PETSc options file      = %s \n ", petscfile.c_str());


    StabFEM  stabfem;

    stabfem.readInputData(meshfile);

    stabfem.readControlParameters(controlfile);

    stabfem.prepareInputData();

    double timerStart, timerEnd;
    int parm[3];

    timerStart = MPI_Wtime();
    stabfem.setSolver(1, parm, false);
    timerEnd = MPI_Wtime(); 
    PetscPrintf(MPI_COMM_WORLD, "\n\n Elapsed time = %f seconds \n\n", timerEnd - timerStart );

    //stabfem.postProcess();
    //stabfem.diffStiffTest();

    stabfem.solveFullyImplicit();

    //string  outputfile = "solution.dat";
    //stabfem.writeResult(outputfile);
    //stabfem.postProcess();

    //stabfem.printComputerTimes();
    MPI_Barrier(MPI_COMM_WORLD);

    PetscPrintf(MPI_COMM_WORLD, "\n\n\n Program is successful \n\n\n ");

    stabfem.deallocatePetscObjects();
    MPI_Barrier(MPI_COMM_WORLD);

    PetscFinalize(); //CHKERRQ(ierr);

    return 0;
}

