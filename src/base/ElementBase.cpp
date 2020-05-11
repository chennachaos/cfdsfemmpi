
#include "ElementBase.h"
//#include "BasisFunctionsLagrange.h"
#include "SolutionData.h"

#include "KimMoinFlow.h"

using namespace std;



ElementBase::ElementBase()
{
  // cout << "     ElementBase: constructor ...\n\n";

  nlbf = ndof = nsize = nGP = 0;
  subdomId = 0;
}



ElementBase::~ElementBase()
{
//  cout << "     ElementBase: destructor ...\n\n";
}



