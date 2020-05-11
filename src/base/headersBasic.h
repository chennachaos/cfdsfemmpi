
#ifndef incl_headersBasic_h
#define incl_headersBasic_h


#include <iostream>
#include <limits.h>
#include <float.h>
#include <vector>
#include <assert.h>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <set>
#include <fstream>
#include <string.h>
#include <iomanip>
#include <stdlib.h>
#include <algorithm>

#include <numeric>
#include <iterator>
#include <functional>


using namespace std;

const double PI = 3.141592653589793238462643383279502884197169399375105;



template<typename T>
void findUnique(vector<T>& vect)
{
  sort(vect.begin(), vect.end());
  vect.erase(unique(vect.begin(), vect.end()), vect.end());
}





#endif
