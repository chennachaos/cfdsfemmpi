
#ifndef incl_SolverPetsc_h
#define incl_SolverPetsc_h

#include "headersBasic.h"
#include "headersEigen.h"
#include "petscksp.h"
#include "petscmat.h"


using namespace std;

enum { SOLVER_EMPTY, PATTERN_OK, INIT_OK, ASSEMBLY_OK, FACTORISE_OK };

class SolverPetsc
{
  public:

    Vec  rhsVec, solnVec, solnVecPrev, reacVec;
    Mat  mtx; // linear system matrix
    KSP  ksp; // linear solver context
    PC   pc; // preconditioner context

    PetscInt nRow, nCol, nnz;

    int currentStatus;

    bool  checkIO;

    PetscReal norm; // norm of solution error

    PetscErrorCode errpetsc;

    ////////////////////////////
    //
    // member functions
    //
    ///////////////////////////

    SolverPetsc();

    virtual ~SolverPetsc();

    virtual int initialise(int size_local, int size_global, int* diag_nnz, int* offdiag_nnz);

    int setSolverAndParameters();

    virtual int zeroMtx();

    virtual int free();

    virtual int printMatrix(int dig=8, int dig2=4, bool gfrmt=true, int indent = 0, bool interactive = false);

    virtual double getMatrixCoefficient(int,int);

    virtual int assembleMatrixAndVectorSerial(vector<int>& forAssyElem, MatrixXd& Klocal, VectorXd& Flocal);

    virtual int assembleMatrixAndVectorParallel(vector<int>& forAssyElem, vector<int>& dof_map, MatrixXdRM& Klocal, VectorXd& Flocal);

    virtual int factorise();

    virtual int solve();

    virtual int factoriseAndSolve();
};



#endif

