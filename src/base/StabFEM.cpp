
#include <algorithm>
#include <chrono>
#include "StabFEM.h"
#include "KimMoinFlow.h"
#include "UtilitiesGeneral.h"
#include "LagrangeElem2DNavierStokesTria3Node.h"
#include "LagrangeElem2DNavierStokesQuad4Node.h"
#include "LagrangeElem3DNavierStokesTetra4Node.h"
#include "LagrangeElem3DNavierStokesHexa8Node.h"


using namespace std;


StabFEM::StabFEM()
{
    MPI_Comm_size(MPI_COMM_WORLD, &n_mpi_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_mpi_proc);

    ndof = 0; nElem_global = 0; nNode_global = 0; npElem = 0; fileCount = 0;
    ntotdofs_local = ntotdofs_global = 0;

    computerTimeAssembly = computerTimeSolver = computerTimePostprocess = computerTimePattern = computerTimeTimeLoop = 0.0;

    AlgoType = 2;

    elems = nullptr;

    solverPetsc = nullptr;
}


StabFEM::~StabFEM()
{
    //cout << " StabFEM::~StabFEM() " << this_mpi_proc << endl;
    if(elems != nullptr)
    {
      for(int ii=0;ii<nElem_global;++ii)
        delete elems[ii];

      delete [] elems;
      elems = nullptr;
    }

    if(elemsFaces != nullptr)
    {
      for(int ii=0;ii<ElemFaceLoadData.size();++ii)
        delete elemsFaces[ii];

      delete [] elemsFaces;
      elemsFaces = nullptr;
    }

    if(solverPetsc != nullptr)
    {
      delete solverPetsc;
      solverPetsc = nullptr;
    }

    //cout << " StabFEM::~StabFEM() " << this_mpi_proc << endl;
}




int  StabFEM::deallocatePetscObjects()
{
  if(solverPetsc != nullptr)
    solverPetsc->free();

  return 1;
}





void  StabFEM::readInputData(string&  fname)
{
    PetscPrintf(MPI_COMM_WORLD, "\n\n Reading input data \n\n");

    infilename = fname;

    std::ifstream  infile( string(infilename+".dat") );

    if(infile.fail())
    {
        PetscPrintf(MPI_COMM_WORLD, "\n\n Could not open the input nodes file \n\n");
        exit(1);
    }

    string  line, stringVal, stringVec[10];
    int  ii, arrayInt[100];
    double  tempDbl;

    // read the dimension
    infile >> stringVal >> ndim;
    ndof = ndim+1;

    // read number of nodes per element
    infile >> stringVal >> npElem;

    // read number of nodes
    infile >> stringVal >> nNode_global;

    // read number of elements
    infile >> stringVal >> nElem_global;

    // read number of Dirichlet BCs
    infile >> stringVal >> nDBC;

    // read number of Force BCs
    infile >> stringVal >> nFBC;

    // read number of Output nodes
    infile >> stringVal >> nOutputFaceLoads;

    if(this_mpi_proc == 0)
    {
      cout << " ndim              =  " << ndim << endl;
      cout << " ndof              =  " << ndof << endl;
      cout << " nNode_global      =  " << nNode_global << endl;
      cout << " nElem_global      =  " << nElem_global << endl;
      cout << " npElem            =  " << npElem << endl;
      cout << " nDBC              =  " << nDBC << endl;
      cout << " nFBC              =  " << nFBC << endl;
      cout << " nOutputFaceLoads  =  " << nOutputFaceLoads << endl;
    }

    // read nodal coordinates
    ////////////////////////////////////////////

    PetscPrintf(MPI_COMM_WORLD, "\n\n reading nodes \n\n");

    node_coords.resize(nNode_global);
    for(ii=0; ii<nNode_global; ++ii)
      node_coords[ii].resize(ndim);

    infile >> stringVal ;
    PetscPrintf(MPI_COMM_WORLD, "\n\n reading, %s \n\n", stringVal.c_str());
    if(ndim == 2)
    {
      for(ii=0; ii<nNode_global; ++ii)
      {
        infile >> tempDbl >> node_coords[ii][0] >> node_coords[ii][1] ;
      }
    }
    else
    {
      for(ii=0; ii<nNode_global; ++ii)
      {
        infile >> tempDbl >> node_coords[ii][0] >> node_coords[ii][1] >> node_coords[ii][2];
      }
    }

    // read elements
    ////////////////////////////////////////////
    infile >> stringVal ;
    PetscPrintf(MPI_COMM_WORLD, "\n\n reading %s \t %d \n\n", stringVal.c_str(), npElem);

    elemConn.resize(nElem_global);

    if(npElem == 3)
    {
      for(int ee=0; ee<nElem_global; ++ee)
      {
        elemConn[ee].resize(npElem);

        infile >> arrayInt[0] >> arrayInt[1] >> arrayInt[2] >> arrayInt[3] >> arrayInt[4] >> arrayInt[5] >> arrayInt[6];

        for(ii=0; ii<npElem; ++ii)
          elemConn[ee][ii] = arrayInt[4+ii]-1;
      }
    }
    else if(npElem == 4)
    {
      for(int ee=0; ee<nElem_global; ++ee)
      {
        elemConn[ee].resize(npElem);

        infile >> arrayInt[0] >> arrayInt[1] >> arrayInt[2] >> arrayInt[3] >> arrayInt[4] >> arrayInt[5] >> arrayInt[6] >> arrayInt[7];

        for(ii=0; ii<npElem; ++ii)
          elemConn[ee][ii] = arrayInt[4+ii]-1;
      }
    }
    else if(npElem == 8)
    {
      for(int ee=0; ee<nElem_global; ++ee)
      {
        elemConn[ee].resize(npElem);

        infile >> arrayInt[0] >> arrayInt[1] >> arrayInt[2] >> arrayInt[3] >> arrayInt[4] >> arrayInt[5] >> arrayInt[6] >> arrayInt[7] >> arrayInt[8] >> arrayInt[9] >> arrayInt[10] >> arrayInt[11] ;

        for(ii=0; ii<npElem; ++ii)
          elemConn[ee][ii] = arrayInt[4+ii]-1;
      }
    }
    else
    {
      cerr << " Invalid npElem " << npElem << endl;
      exit(-1);
    }

    //
    // Read Dirichlet BC data
    //
    ////////////////////////////////////////////

    vector<double>  vecDblTemp(3);

    infile >> stringVec[0] >> stringVec[1] >> stringVec[2] ;

    for(ii=0; ii<nDBC; ++ii)
    {
      infile >> arrayInt[0] >> arrayInt[1] >> tempDbl ;

      vecDblTemp[0] = arrayInt[0]-1;
      vecDblTemp[1] = arrayInt[1]-1;
      vecDblTemp[2] = tempDbl;

      DirichletBCs.push_back(vecDblTemp);
    }

    //
    // Read Output data
    //
    ////////////////////////////////////////////

    infile >> stringVal ;

    outputEdges.resize(nOutputFaceLoads);
    for(ii=0; ii<nOutputFaceLoads; ++ii)
    {
      outputEdges[ii].resize(1);

      //infile >> arrayInt[0] >> arrayInt[1];
      infile >> arrayInt[0];

      outputEdges[ii][0] = arrayInt[0]-1;
      //outputEdges[ii][1] = arrayInt[1]-1;
    }

    infile.close();


    fout_convdata.open(string(infilename+"-forces.dat"), ios::out | ios::trunc );

    if(fout_convdata.fail())
    {
       cout << " Could not open the output file for writing forces " << endl;
       exit(1);
    }

    fout_convdata.setf(ios::fixed);
    fout_convdata.setf(ios::showpoint);
    fout_convdata.precision(14);

    PetscPrintf(MPI_COMM_WORLD, "\n\n Input files have been read successfully \n\n");

    return;
}



void StabFEM::readControlParameters(string& fname)
{
    PetscPrintf(MPI_COMM_WORLD, "\n\n StabFEM::readControlParameters \n\n");

    ifstream  infile(fname);
    if(infile.fail())
    {
       cout << " Could not open 'control-parameters.dat' file " << endl;
       exit(-1);
    }

    //time integration parameters
    timeData[1] = 1.0;   timeData[2] = 0.0;

    string  stringVal;

    //density
    infile >> stringVal >> elemData[0];

    //viscosity
    infile >> stringVal >> elemData[1];

    //Body force in X-, Y- and Z- direction
    infile >> stringVal >> elemData[2] >> elemData[3] >> elemData[4];

    //Stabilisation: SUPG, PSPG, LSIC
    infile >> stringVal >> elemData[8] >> elemData[9] >> elemData[10];

    //tis
    infile >> stringVal >> tis;

    //rhoInf = 0.0;
    infile >> stringVal >> rhoInf;

    // timestep
    infile >> stringVal >> dt;

    // final time
    infile >> stringVal >> timeFinal;

    // Maximum number of time steps
    infile >> stringVal >> stepsMax;

    // Output file frequency
    infile >> stringVal >> outputFreq;

    // convergence tolerance
    infile >> stringVal >> conv_tol;

    infile.close();

    PetscPrintf(MPI_COMM_WORLD, "\n\n Control parameters are successfully read \n\n");

    return;
}





void StabFEM::prepareInputData()
{
    PetscPrintf(MPI_COMM_WORLD, "\n\n  StabFEM::prepareInputData()  .... STARTED ...\n\n");

    int ii, jj, kk, ee, nn, ind, n1, n2, dof;

    assert(ndim > 0 && ndim < 4);

    // ==================================================
    //
    // Check the  consistency of input data
    //
    // ==================================================

    //checkInputData();

    ///////////////////////////////////////////////////////////////////
    //
    ///////////////////////////////////////////////////////////////////

    // create elements and prepare element data
    elems = new ElementBase* [nElem_global];

    for(ee=0;ee<nElem_global;++ee)
    {
      if(ndim == 2)
      {
        if(npElem == 3)
          elems[ee] = new LagrangeElem2DNavierStokesTria3Node;
        else if(npElem == 4)
          elems[ee] = new LagrangeElem2DNavierStokesQuad4Node;
      }
      else
      {
        if(npElem == 4)
          elems[ee] = new LagrangeElem3DNavierStokesTetra4Node;
        else if(npElem == 8)
          elems[ee] = new LagrangeElem3DNavierStokesHexa8Node;
      }

      elems[ee]->SolnData = &(SolnData);

      elems[ee]->nodeNums = elemConn[ee];
      //elems[ee]->prepareElemData(node_coords);
    }


    ///////////////////////////////////////////////////////////////////
    //
    // set SolnData details
    //
    ///////////////////////////////////////////////////////////////////

    SolnData.initialise(nNode_global*ndof);

    PetscPrintf(MPI_COMM_WORLD, "\n\n  StabFEM::prepareInputData()  .... FINISHED ...\n\n");

    return;
}




void StabFEM::assignBoundaryConditions(double timeCur, double dt, double timeFact)
{
    int ii, n1, n2, ind;
    for(ii=0; ii<nDBC; ++ii)
    {
        n1 = int(DirichletBCs[ii][0]);
        n2 = int(DirichletBCs[ii][1]);

        ind = n1*ndof+n2;

        SolnData.solnApplied[ind] = DirichletBCs[ii][2] * timeFact - SolnData.soln[ind];
        //cout << ii << '\t' << n1 << '\t' << n2 << '\t' << ind << '\t' << solnApplied[ind] << endl;

        //solnApplied[ind] = analy.computeValue(n2, node_coords[n1][0], node_coords[n1][1], 0.0, timeNow) - soln[ind];
    }
    //printVector(solnApplied);

    return;
}




void StabFEM::setInitialConditions()
{
    double  xx=0.0, yy=0.0, zz=0.0, fact;

    for(int ii=0; ii<nNode_global; ++ii)
    {
        xx = node_coords[ii][0];
        yy = node_coords[ii][1];
        zz = node_coords[ii][2];

        //veloPrev(ii*2) = 2.0*yy*(3.0-yy)/3.0;
        SolnData.solnPrev(ii*ndim) = 1.0*yy;

        //solnPrev(ii*ndim) = 16.0*0.45*yy*zz*(0.41-yy)*(0.41-zz)/0.41/0.41/0.41/0.41;
    }
    SolnData.soln = SolnData.solnPrev;

    return;
}






void StabFEM::setTimeParam()
{
  //SolnData.setTimeParam();

  return;
}




void StabFEM::writeNodalData()
{
  return;
}



int  StabFEM::writeResult(string&  fname)
{
    ofstream fout(fname);

    if(fout.fail())
    {
       cout << " Could not open the file in 'StabFEM::writeResult' ... " << endl;
       exit(1);
    }

    fout.setf(ios::fixed);
    fout.setf(ios::showpoint);
    fout.precision(8);

    fout << nNode_global <<endl;

    int  ii, nn, ind;
    if(ndim == 2)
    {
      for(int ii=0; ii<nNode_global; ii++)
      {
        ind = node_map_get_new[ii]*ndof;
        fout << SolnData.soln[ind] << '\t' << SolnData.soln[ind+1] << '\t' << SolnData.soln[ind+2] << endl;
      }
    }
    else
    {
      for(int ii=0; ii<nNode_global; ii++)
      {
        ind = node_map_get_new[ii]*ndof;
        fout << SolnData.soln[ind] << '\t' << SolnData.soln[ind+1] << '\t' << SolnData.soln[ind+2] << '\t' << SolnData.soln[ind+3] << endl;
      }
    }

    fout.close();

    return  0;
}






int  StabFEM::checkResult(string&  fname)
{
    std::ifstream  datfile(fname);

    if(datfile.fail())
    {
        cout << " Could not open the file in 'StabFEM::checkResult' ... " << endl;
        exit(1);
    }

    string  line, stringVal, stringVec[10];
    double  tempDbl;

    int  nnodes, ii, ind;
    // read number of nodes
    datfile  >>  nnodes;
    if(nnodes != nNode_global)
      return  -1;

    double  val[ndof];
    bool  flag = true;
    if(ndim == 2)
    {
      for(ii=0; ii<nNode_global; ++ii)
      {
        datfile >> val[0] >> val[1] >> val[2];
        //cout << val[0] << '\t' << val[1] << '\t' << val[2] << endl;

        ind = node_map_get_new[ii]*ndof;
        val[0] -= SolnData.soln[ind];
        val[1] -= SolnData.soln[ind+1];
        val[2] -= SolnData.soln[ind+2];

        if( sqrt(val[0]*val[0] + val[1]*val[1] + val[2]*val[2]) > 1.0e-4)
        {
          flag = false;
          break;
        }
      }
    }
    else
    {
      for(ii=0; ii<nNode_global; ++ii)
      {
        datfile >> val[0] >> val[1] >> val[2] >> val[3];
      }
    }

    datfile.close();
    if(!flag)
      return -2;

    return  0;
}




int  StabFEM::readResult(string&  fname)
{
    ofstream fout(fname);

    if(fout.fail())
    {
       cout << " Could not open the file in 'StabFEM::readResult' ... " << endl;
       exit(1);
    }

    fout.setf(ios::fixed);
    fout.setf(ios::showpoint);
    fout.precision(8);


  return  0;
}





int StabFEM::diffStiffTest()
{
    SolnData.setTimeParam();

    int  elemnum = 0;
    cout << " Diff Test for element " << elemnum << endl;
    int  nlbf = elems[elemnum]->nlbf;
    int  nsize = elems[elemnum]->nsize;

    double ddd = 1.0;
    int    i, ii, jnd, jj, j, k, ind, bb, nn;

    double dd[6]  = {-3.*ddd, -2.*ddd, -ddd, +ddd, +2.*ddd, +3.*ddd };
    double  r[6*nsize], timeNow = 1.0;

    MatrixXd  Klocal1(nsize, nsize), Klocal2(nsize, nsize), Kdiff(nsize, nsize);
    VectorXd  Flocal(nsize);

    Klocal1.setZero();
    Klocal2.setZero();
    Flocal.setZero();

    // loop over columns of s

    for (jnd=0; jnd<nlbf; jnd++) // nodes
    {
      for (jj=0; jj<ndof; jj++)  // dof
      {
        j = jnd*ndof + jj;

        nn  = elems[elemnum]->nodeNums[jnd];
        ind = nn*ndof+jj;

        // loop over perturbations

        for (k=0; k<6; k++)
        {
          // apply pertubation
          SolnData.soln[ind] += dd[k];

          SolnData.solnCur[ind]  =  SolnData.soln[ind];

          // calculate residual
          elems[elemnum]->calcStiffnessAndResidual(node_coords, elemData, Klocal1, Flocal, timeNow);
          //printMatrix(Klocal1);
          //printVector(Flocal);

          // remove pertubation
          SolnData.soln[ind] -= dd[k];

          SolnData.solnCur[ind]  =  SolnData.soln[ind];

          // loop over rows of s, store residual

          for(i=0; i<nsize; i++)
            r[i*6+k] = Flocal[i];
        }

        // loop over rows of s
        for (i=0; i<nsize; i++)
        Klocal2(j,i) = (   +       r[i*6+0]
                           -  9. * r[i*6+1]
                           + 45. * r[i*6+2]
                           - 45. * r[i*6+3]
                           +  9. * r[i*6+4]
                           -       r[i*6+5] ) / (60. * ddd);
    }
  }

  cout << " performing Difference calculations " << endl;

  // calculate stiffness
  elems[elemnum]->calcStiffnessAndResidual(node_coords, elemData, Klocal1, Flocal, timeNow);

  cout.setf(ios::fixed);
  cout.setf(ios::showpoint);
  cout.precision(5);

  cout << "    Analytical Stiffness Matrix   " << endl;
  cout << endl;
  for(ii=0;ii<nsize;ii++)
  {
    cout << '\t' ;
    for(jj=0;jj<nsize;jj++)
    {
      cout  <<  fixed << Klocal1(ii,jj) << "   " ;
    }
    cout << endl;
  }
  cout << endl;
  cout << endl;


  cout << "    diff Stiffness Matrix   " << endl;
  cout << endl;
  cout << endl;
  for(ii=0;ii<nsize;ii++)
  {
    cout << '\t' ;
    for(jj=0;jj<nsize;jj++)
    {
      Kdiff(ii,jj) = Klocal1(ii,jj) - Klocal2(ii,jj);
      cout  <<  fixed <<  Klocal2(ii,jj) << "   " ;
    }
    cout << endl;
  }
  cout << endl;
  cout << endl;

  cout << "    difference of Stiffness Matrices   " << endl;
  cout << endl;
  cout << endl;
  for(ii=0;ii<nsize;ii++)
  {
    cout << '\t' ;
    for(jj=0;jj<nsize;jj++)
    {
      cout  <<  fixed << Kdiff(ii, jj) << "   " ;
    }
    cout << endl;
  }
  cout << endl;
  cout << endl;

  return 0;
}




int StabFEM::printComputerTimes()
{
    PetscPrintf(MPI_COMM_WORLD, "\n====================================================================\n");
    PetscPrintf(MPI_COMM_WORLD, "\n       Computer time in seconds \n");
    PetscPrintf(MPI_COMM_WORLD, "\n====================================================================\n");

    PetscPrintf(MPI_COMM_WORLD, "\n\n     (1) Matrix Pattern                  = %12.6f \n", computerTimePattern );
    PetscPrintf(MPI_COMM_WORLD, "\n\n     (2) Matrix Assembly                 = %12.6f \n", computerTimeAssembly );
    PetscPrintf(MPI_COMM_WORLD, "\n\n     (3) Matrix Solver                   = %12.6f \n", computerTimeSolver );
    PetscPrintf(MPI_COMM_WORLD, "\n\n     (4) Post process                    = %12.6f \n", computerTimePostprocess );
    PetscPrintf(MPI_COMM_WORLD, "\n\n     (5) Total time in time loop         = %12.6f \n", computerTimeTimeLoop );
    PetscPrintf(MPI_COMM_WORLD, "\n\n     (6) Difference (5-2-3-4)            = %12.6f \n", (computerTimeTimeLoop-computerTimeAssembly-computerTimeSolver-computerTimePostprocess) );

    PetscPrintf(MPI_COMM_WORLD, "\n\n====================================================================\n\n");
}





