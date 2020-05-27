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
    //Set the input file name
    //The file name is specified from the command line
    if(argc == 0)
    {
        cerr << " Error in input data " << endl;
        cerr << " You must enter name of input file" << endl;
        cerr << " Aborting..." << endl;
    }

    PetscInitialize(NULL, NULL, "petsc_options.dat", NULL);

    string  meshfile = argv[1];
    string  controlfile = argv[2];


    StabFEM  stabfem;

    stabfem.readInputData(meshfile);

    stabfem.readControlParameters(controlfile);

    stabfem.prepareInputData();

    //stabfem.partitionMesh();

    int parm[3];
    stabfem.setSolver(1, parm, false);

    //stabfem.solveFullyImplicit();

    stabfem.postProcess();

    cout << " Program is successful \n " << endl;

    PetscFinalize(); //CHKERRQ(ierr);

    return 0;
}

