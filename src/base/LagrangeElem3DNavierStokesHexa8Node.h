
#ifndef incl_LagrangeElem3DNavierStokesHexa8Node_h
#define incl_LagrangeElem3DNavierStokesHexa8Node_h


#include "ElementBase.h"


class  LagrangeElem3DNavierStokesHexa8Node : public ElementBase
{
  public:

    LagrangeElem3DNavierStokesHexa8Node();

    virtual ~LagrangeElem3DNavierStokesHexa8Node();

    void prepareElemData(vector<vector<double> >& node_coods);

    virtual int  calcStiffnessAndResidual(vector<vector<double> >& node_coords, double* elemData, MatrixXd& Klocal, VectorXd& Flocal, double timeCur);

    virtual  int calcError(int index);
};










#endif

