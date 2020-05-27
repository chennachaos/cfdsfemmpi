
#ifndef incl_LagrangeElem2DNavierStokesTria3Node_h
#define incl_LagrangeElem2DNavierStokesTria3Node_h


#include "ElementBase.h"


class  LagrangeElem2DNavierStokesTria3Node : public ElementBase
{
  public:

    LagrangeElem2DNavierStokesTria3Node();

    virtual ~LagrangeElem2DNavierStokesTria3Node();

    void prepareElemData(vector<vector<double> >& node_coods);

    virtual int  calcStiffnessAndResidual(vector<vector<double> >& node_coords, double* elemData, MatrixXdRM& Klocal, VectorXd& Flocal, double timeCur);

    virtual  int calcError(int index);
};










#endif

