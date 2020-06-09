
#ifndef incl_LagrangeElem3DNavierStokesTetra4Node_h
#define incl_LagrangeElem3DNavierStokesTetra4Node_h


#include "ElementBase.h"


class  LagrangeElem3DNavierStokesTetra4Node : public ElementBase
{
  public:

    LagrangeElem3DNavierStokesTetra4Node();

    virtual ~LagrangeElem3DNavierStokesTetra4Node();

    void prepareElemData(vector<vector<double> >& node_coods);

    virtual int  calcStiffnessAndResidual(vector<vector<double> >& node_coords, double* elemData, MatrixXd& Klocal, VectorXd& Flocal, double timeCur);

    virtual  int calcError(int index);
};










#endif

