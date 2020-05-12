
#include <algorithm>
#include <chrono>
#include "StabFEM.h"
#include "KimMoinFlow.h"
#include "UtilitiesGeneral.h"
//#include "SolutionData.h"
//#include "LagrangeElem2DNavierStokesTria3Node.h"
#include "LagrangeElem2DNavierStokesQuad4Node.h"
//#include "LagrangeElem3DNavierStokesTetra4Node.h"
//#include "LagrangeElem3DNavierStokesHexa8Node.h"


using namespace std;


StabFEM::StabFEM()
{
    ndof = 0; nElem = 0; nNode = 0; npElem = 0; fileCount = 0;
    totalDOF = 0;

    AlgoType = 2;

    elems = nullptr;

    solverPetsc = NULL;
}


StabFEM::~StabFEM()
{
    if(elems != NULL)
    {
      for(int ii=0;ii<nElem;++ii)
        delete elems[ii];

      delete [] elems;
      elems = NULL;
    }
}




void  StabFEM::readInputData(string&  fname)
{
    cout << " Reading input data \n\n " << endl;

    infilename = fname;

    std::ifstream  infile( string(infilename+".dat") );

    if(infile.fail())
    {
        cout << " Could not open the input nodes file " << endl;
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
    infile >> stringVal >> nNode;

    // read number of elements
    infile >> stringVal >> nElem;

    // read number of Dirichlet BCs
    infile >> stringVal >> nDBC;

    // read number of Force BCs
    infile >> stringVal >> nFBC;

    // read number of Output nodes
    infile >> stringVal >> nOutputFaceLoads;

    cout << " ndim              =  " << ndim << endl;
    cout << " nNode             =  " << nNode << endl;
    cout << " nElem             =  " << nElem << endl;
    cout << " npElem            =  " << npElem << endl;
    cout << " nDBC              =  " << nDBC << endl;
    cout << " nFBC              =  " << nFBC << endl;
    cout << " nOutputFaceLoads  =  " << nOutputFaceLoads << endl;

    // read nodal coordinates
    ////////////////////////////////////////////

    cout << " reading nodes " << endl;

    node_coords.resize(nNode);
    for(ii=0; ii<nNode; ++ii)
      node_coords[ii].resize(ndim);

    infile >> stringVal ;
    cout << " reading " << stringVal << endl;
    if(ndim == 2)
    {
      for(ii=0; ii<nNode; ++ii)
      {
        infile >> tempDbl >> node_coords[ii][0] >> node_coords[ii][1] ;
      }
    }
    else
    {
      for(ii=0; ii<nNode; ++ii)
      {
        infile >> tempDbl >> node_coords[ii][0] >> node_coords[ii][1] >> node_coords[ii][2];
      }
    }

    // read elements
    ////////////////////////////////////////////
    infile >> stringVal ;
    cout << " reading " << stringVal << '\t' << npElem << endl;

    elemConn.resize(nElem);

    if(npElem == 3)
    {
      for(int ee=0; ee<nElem; ++ee)
      {
        elemConn[ee].resize(npElem);

        infile >> arrayInt[0] >> arrayInt[1] >> arrayInt[2] >> arrayInt[3] >> arrayInt[4] >> arrayInt[5] >> arrayInt[6];

        //printf("%6d \t %6d \t %6d \t %6d \n", arrayInt[4], arrayInt[5], arrayInt[6], arrayInt[7]);

        for(ii=0; ii<npElem; ++ii)
          elemConn[ee][ii] = arrayInt[4+ii]-1;
      }
    }
    else if(npElem == 4)
    {
      for(int ee=0; ee<nElem; ++ee)
      {
        elemConn[ee].resize(npElem);

        infile >> arrayInt[0] >> arrayInt[1] >> arrayInt[2] >> arrayInt[3] >> arrayInt[4] >> arrayInt[5] >> arrayInt[6] >> arrayInt[7];

        for(ii=0; ii<npElem; ++ii)
          elemConn[ee][ii] = arrayInt[4+ii]-1;
      }
    }
    else if(npElem == 8)
    {
      for(int ee=0; ee<nElem; ++ee)
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
    cout << " reading " << stringVec[0] << stringVec[1] << stringVec[2] << endl;

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
    cout << " reading " << stringVal << endl;

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
       cout << " Could not open the Output file" << endl;
       exit(1);
    }

    fout_convdata.setf(ios::fixed);
    fout_convdata.setf(ios::showpoint);
    fout_convdata.precision(14);


    cout << " Input files have been read successfully \n\n " << endl;

    return;
}



void StabFEM::readControlParameters(string& fname)
{
    cout << " StabFEM::readControlParameters " << endl;

    //ifstream  infile("control-parameters.dat");
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

    cout << " Control parameters are successfully read " << endl;

    return;
}





void StabFEM::prepareInputData()
{
    printf("\n     StabFEM::prepareInputData()  .... STARTED ...\n");

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
    elems = new ElementBase* [nElem];

    for(ee=0;ee<nElem;++ee)
    {
      if(ndim == 2)
      {
        //if(npElem == 3)
          //elems[ee] = new LagrangeElem2DNavierStokesTria3Node;
        //else if(npElem == 4)
          elems[ee] = new LagrangeElem2DNavierStokesQuad4Node;
      }
      else
      {
        //if(npElem == 4)
          //elems[ee] = new LagrangeElem3DNavierStokesTetra4Node;
        //else if(npElem == 8)
          //elems[ee] = new LagrangeElem3DNavierStokesHexa8Node;
      }

      elems[ee]->SolnData = &(SolnData);

      elems[ee]->nodeNums = elemConn[ee];

      elems[ee]->prepareElemData(node_coords);
    }


    cout << " elements are created and prepared " << endl;

    ///////////////////////////////////////////////////////////////////
    //
    // set SolnData details
    //
    ///////////////////////////////////////////////////////////////////

    SolnData.initialise(nNode*ndof);

    double  xx, yy, zz, fact;

    cout << " aaaaaaaaaaaaaaa " << endl;

    Kovasznay analy;

    SolnData.solnApplied.setZero();
    for(ii=0; ii<nDBC; ++ii)
    {
        n1 = DirichletBCs[ii][0];
        n2 = DirichletBCs[ii][1];

        jj = n1*ndof+n2;

        //xx = node_coords[n1][0] ;
        //yy = node_coords[n1][1] ;

        //DirichletBCsVelo[ii][2] = analy.computeValue(n2, xx, yy);

        SolnData.solnApplied(jj) = DirichletBCs[ii][2];
    }
    //printVector(SolnData.solnApplied);

    printf("     StabFEM::prepareInputData()  .... FINISHED ...\n\n");

    return;
}




void StabFEM::assignBoundaryConditions(double timeCur, double dt, double timeFact)
{
    int ii, n1, n2, ind;
    for(ii=0; ii<nDBC; ++ii)
    {
        n1 = DirichletBCs[ii][0];
        n2 = DirichletBCs[ii][1];

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

    for(int ii=0; ii<nNode; ++ii)
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



