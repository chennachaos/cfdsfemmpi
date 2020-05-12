
#ifndef incl_LagrangeElem2DNavierStokesQuad4Node_h
#define incl_LagrangeElem2DNavierStokesQuad4Node_h


#include "ElementBase.h"


class  LagrangeElem2DNavierStokesQuad4Node : public ElementBase
{
  public:

    LagrangeElem2DNavierStokesQuad4Node();

    virtual ~LagrangeElem2DNavierStokesQuad4Node();

    void prepareElemData(vector<vector<double> >& node_coods);


    virtual int  calcStiffnessAndResidual(vector<vector<double> >& node_coords, double* elemData, MatrixXdRM& Klocal, VectorXd& Flocal, double timeCur);

    virtual  int calcError(int index);
};










#endif

