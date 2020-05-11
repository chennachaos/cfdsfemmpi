
#ifndef incl_BasisFunctionsLagrange_h
#define incl_BasisFunctionsLagrange_h


/**
  Computes function values and the first derivative (wrt to spatial coordinates)
  for basis functions for the 2D problem
*/
int computeBasisFunctions2D(bool flag, int ETYPE, int degree, double* param, double* xNode, double* yNode, double*  N, double*  dN_dx, double* dN_dy, double&  Jac);


/**
  Computes function values and the first derivative (wrt to spatial coordinates)
  for basis functions for the 3D problem
*/
int computeBasisFunctions3D(bool flag, int ETYPE, int degree, double* param, double* xNode, double* yNode, double* zNode, double*  N, double*  dN_dx, double* dN_dy, double* dN_dz, double&  Jac);

/**
  Computes function values, first derivative and second derivative of
  univariate Lagrange polynomials of a given order p for
*/
void Lagrange_BasisFuns1D(int p, double xi, double* N, double* dN_dxi, double* d2N_dxi2);

/**
  Computes function values and first derivative of
  univariate Lagrange polynomials of a given order p
*/
void Lagrange_BasisFuns1D(int p, double xi, double* N, double* dN_dxi);

/**
  Computes function values of
  univariate Lagrange polynomials of a given order p
*/
void Lagrange_BasisFuns1D(int p, double xi, double* N);


/**
  Computes function values of
  Lagrange polynomials of a given order p for Triangular elements
*/
void LagrangeBasisFunsTria(int p, double xi, double zeta, double* N);

/**
  Computes function values and first derivative of
  Lagrange polynomials of a given order p for Triangular elements
*/
void LagrangeBasisFunsTria(int p, double xi, double zeta, double* N, double* dN_dxi, double* dN_dzeta);


/**
  Computes function values of
  Lagrange polynomials of a given order p for Quadrilateral elements
*/
void LagrangeBasisFunsQuad(int p, double xi, double zeta, double* N);

/**
  Computes function values and first derivative of
  Lagrange polynomials of a given order p for Quadrilateral elements
*/
void LagrangeBasisFunsQuad(int p, double xi, double zeta, double* N, double* dN_dxi, double* dN_dzeta);


/**
  Computes function values of
  Lagrange polynomials of a given order p for Tetrahedral elements
*/
void LagrangeBasisFunsTetra(int p, double xi1, double xi2, double xi3, double* N);

/**
  Computes function values and first derivative of
  Lagrange polynomials of a given order p for Tetrahedral elements
*/
void LagrangeBasisFunsTetra(int p, double xi, double zeta, double eta, double* N, double* dN_dxi, double* dN_dzeta, double* dN_deta);

/**
  Computes function values of
  Lagrange polynomials of a given order p for Hexahedral elements
*/
void LagrangeBasisFunsHexa(int p, double xi, double zeta, double eta, double* N);

/**
  Computes function values and first derivative of
  Lagrange polynomials of a given order p for Hexahedral elements
*/
void LagrangeBasisFunsHexa(int p, double xi, double zeta, double eta, double* N, double* dN_dxi, double* dN_dzeta, double* dN_deta);

/**
  Computes function values and first derivative of
  Lagrange polynomials of a given order p for Penta/wedge elements
*/
void LagrangeBasisFunsPrism(int p, double xi, double zeta, double eta, double* N, double* dN_dxi, double* dN_dzeta, double* dN_deta);

/**
  Computes function values and first derivative of
  Lagrange polynomials of a given order p for Pyramid elements
*/
void LagrangeBasisFunsPyramid(int p, double xi, double zeta, double eta, double* N, double* dN_dxi, double* dN_dzeta, double* dN_deta);


void LagrangeBasisFunsLine1D(int p, double uu, double *xNode, double *N, double *dN_dx, double& Jac);


/**
  Computes function values of
  Lagrange polynomials of an edge in 2D.
  This function is for computing the required basis functions for calculating
  the contribution from traction forces and/or applying Dirichlet BCs using Penalty/Nitsche method
*/
void LagrangeBasisFunsEdge2D(int p, double* param, double *xNode, double* yNode, double *N, double *normal, double& Jac);

/**
  Computes function values of
  Lagrange polynomials of an edge in 3D
  This function is for computing the required basis functions for calculating
  the contribution from traction forces and/or applying Dirichlet BCs using Penalty/Nitsche method
*/
void LagrangeBasisFunsEdge3D(int p, double* param, double *xNode, double* yNode, double *zNode, double *N, double *normal, double& Jac);


/**
  Computes function values of
  Lagrange polynomials of a triangular face in 3D
  This function is for computing the required basis functions for calculating
  the contribution from traction forces and/or applying Dirichlet BCs using Penalty/Nitsche method
*/
void LagrangeBasisFunsFaceTria(int p, double* param, double *xNode, double* yNode, double *zNode, double *N, double *normal, double& Jac);


/**
  Computes function values of
  Lagrange polynomials of a quadrilateral face in 3D
  This function is for computing the required basis functions for calculating
  the contribution from traction forces and/or applying Dirichlet BCs using Penalty/Nitsche method
*/
void LagrangeBasisFunsFaceQuad(int p, double* param, double *xNode, double* yNode, double *zNode, double *N, double *normal, double& Jac);





#endif


