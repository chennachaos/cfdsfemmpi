
#ifndef incl_StabFEM_CLASS_h
#define incl_StabFEM_CLASS_h
#define EIGEN_SUPERLU_SUPPORT

#include <string.h>
#include <vector>
#include <fstream>
#include "ElementBase.h"
#include "SolverPetsc.h"


using std::vector;
using std::cout;
using std::string;
using std::ofstream;

class ElementBase;


#define ELEMENT_TYPE_NAMES_2D {\
                            "LagrangeElem2DNavierStokesTria3Node", \
                            "LagrangeElem2DNavierStokesQuad4Node", NULL};

#define ELEMENT_TYPE_NAMES_3D {\
                            "LagrangeElem3DNavierStokesTetra4Node", \
                            "LagrangeElem3DNavierStokesHexa8Node", NULL};

#define ELEMENT_TYPE_NAMES_FACE {\
                                "BernsteinElem2DEdge3Node", \
                                "BernsteinElem3DFaceTria6Node", \
                                "BernsteinElem3DFaceQuad9Node", NULL};


class StabFEM
{
    //private:

    public:

        int  ndim, ndof, npElem;
        int  nDBC, nFBC, nOutputFaceLoads, fileCount;
        int  stepsMax, outputFreq, tis;
        int  AlgoType;

        PetscInt  nElem_global, nNode_global, nElem_local, nNode_local, nNode_owned;
        PetscInt  n_mpi_procs, this_mpi_proc;
        PetscInt  node_start, node_end, elem_start, elem_end;
        PetscInt  row_start, row_end, ntotdofs_local, ntotdofs_global;

        PetscErrorCode  errpetsc;

        double  conv_tol, rhoInf, timeFinal, dt;
        double  elemData[50], timeData[50];

        vector<vector<double> >  node_coords;               //!< coordinates of the nodes (or control points)
        vector<vector<int> >     outputEdges;               //!< data for computing drag/lift forces
        vector<vector<int> >     elemConn;                  //!< element-node connectivity array

        vector<int>  assyForSoln, OutputNodes;
        vector<int>  node_map_get_old, node_map_get_new;
        vector<int>  dof_map_get_old, dof_map_get_new;
        vector<int>  elem_proc_id, node_proc_id;

        vector<vector<int> >  NodeDofArrayOld, NodeDofArrayNew;
        vector<vector<bool> >  NodeTypeOld, NodeTypeNew;

        vector<vector<double> >  DirichletBCs;          //!< Dirichlet BCs
        vector<vector<double> >  NeumannBCs;                //!< Neumann BCs
        vector<vector<double> >  InitialConds;              //!< Initial conditions
        vector<vector<double> >  OutputData;                //!< data for output
        vector<vector<double> >  nodeForcesData;
        vector<vector<double> >  ElemFaceLoadData;

        SolutionData  SolnData;
        ElementBase  **elems;
        ElementBase  **elemsFaces;

        VectorXd  ForceVectorExternal;
        VectorXd  totalForce, totalMoment, centroid;

        string  infilename;
        ofstream  fout_convdata;

        SolverPetsc  *solverPetsc;

    public:

        StabFEM();

        ~StabFEM();

        ///////////////////////////////////////////////////////////
        //
        // PRE-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        void  assignBoundaryConditions(double timeCur, double dt, double timeFact);

        void  prepareInputData();

        void  readInputData(string& fname);

        void  readControlParameters(string& fname);

        void  writeNodalData();

        int  writeResult(string&);

        int  readResult(string&);

        int  checkResult(string&);

        ///////////////////////////////////////////////////////////
        //
        // SOLUTION PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        int   setSolver(int, int *parm = NULL, bool cIO = false);

        int   prepareMatrixPattern();

        int   partitionMesh();

        int   prepareDataForParallel();

        void  setTimeParam();

        void  timeUpdate();

        void  addExternalForces(double loadFact);

        void  computeElementErrors(int);

        void  setInitialConditions();

        int   solveFullyImplicit();

        int   diffStiffTest();

        ///////////////////////////////////////////////////////////
        //
        // POST-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        void  postProcess();
};






#endif






