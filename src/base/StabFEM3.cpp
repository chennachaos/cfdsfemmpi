
#include "StabFEM.h"
#include "headersVTK.h"


void  StabFEM::postProcess()
{
    //cout << " StabFEM::postProcess " << endl;

    //
    // setup and write vtk data
    //
    //////////////////////////////////////////////


    vtkSmartPointer<vtkUnstructuredGrid>     uGridVTK     =  vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints>               pointsVTK    =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkVertex>               vertexVTK    =  vtkSmartPointer<vtkVertex>::New();

    vtkSmartPointer<vtkTriangle>             triaVTK      =  vtkSmartPointer<vtkTriangle>::New();
    vtkSmartPointer<vtkQuad>                 quadVTK      =  vtkSmartPointer<vtkQuad>::New();
    vtkSmartPointer<vtkTetra>                tetraVTK     =  vtkSmartPointer<vtkTetra>::New();
    vtkSmartPointer<vtkHexahedron>           hexaVTK      =  vtkSmartPointer<vtkHexahedron>::New();

    vtkSmartPointer<vtkFloatArray>           vecVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           scaVTK       =  vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writerUGridVTK =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();


    int  ee, ii, jj, kk, nn, n1, n2;
    double  val, vec[3]={0.0,0.0,0.0};

    vecVTK->SetNumberOfTuples(nNode);
    vecVTK->SetNumberOfComponents(3);
    scaVTK->SetNumberOfTuples(nNode);

    vecVTK->SetName("velocity");
    scaVTK->SetName("pressure");

    vtkIdType pt[10];

    if(ndim == 2)
    {
      for(ii=0; ii<nNode; ii++)
      {
        pt[0] = pointsVTK->InsertNextPoint(node_coords[ii][0], node_coords[ii][1], 0.0);

        kk = ii*ndof;

        vec[0] = soln[kk];
        vec[1] = soln[kk+1];
        val    = soln[kk+2];

        vecVTK->InsertTuple(ii, vec);
        scaVTK->SetTuple1(ii, val);
      }

      for(ee=0; ee<nElem; ee++)
      {
        npElem = elemConn[ee].size();

        if(npElem == 3)
        {
          for(ii=0; ii<npElem; ii++)
            triaVTK->GetPointIds()->SetId(ii, elemConn[ee][ii] );

          uGridVTK->InsertNextCell(triaVTK->GetCellType(), triaVTK->GetPointIds());
        }
        else if(npElem == 4)
        {
          for(ii=0; ii<npElem; ii++)
            quadVTK->GetPointIds()->SetId(ii, elemConn[ee][ii] );

          uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
        }
      }
    }
    else // (ndim == 3)
    {
      for(ii=0; ii<nNode; ii++)
      {
        pt[0] = pointsVTK->InsertNextPoint(node_coords[ii][0], node_coords[ii][1], node_coords[ii][2]);

        kk = ii*ndof;

        vec[0] = soln[kk];
        vec[1] = soln[kk+1];
        vec[2] = soln[kk+2];
        val    = soln[kk+3];

        vecVTK->InsertTuple(ii, vec);
        scaVTK->SetTuple1(ii, val);
      }

      for(ee=0; ee<nElem; ee++)
      {
        npElem = elemConn[ee].size();

        if(npElem == 4)
        {
          for(ii=0; ii<npElem; ii++)
            tetraVTK->GetPointIds()->SetId(ii, elemConn[ee][ii] );

          uGridVTK->InsertNextCell(tetraVTK->GetCellType(), tetraVTK->GetPointIds());
        }
      }
    }

    uGridVTK->SetPoints(pointsVTK);
    uGridVTK->GetPointData()->SetScalars(scaVTK);
    uGridVTK->GetPointData()->SetVectors(vecVTK);

    //Write the file.

    char VTKfilename[200];

    sprintf(VTKfilename,"%s%s%06d%s", infilename.c_str(), "-",fileCount,".vtu");

    writerUGridVTK->SetFileName(VTKfilename);

    writerUGridVTK->SetInputData(uGridVTK);
    writerUGridVTK->Write();


    fileCount = fileCount+1;

    return;
}








