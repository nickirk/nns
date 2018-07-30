#include <stdio.h>
#include <iostream>

#include "../src/HilbertSpace/Basis.hpp"
#include "../src/HilbertSpace/Determinant.hpp"
#include "defaultSystem.hpp"
#include "testComponents.hpp"
//Basis.cxx and Basis.hpp have been tested and passed
using namespace std;
using namespace networkVMC;

int main(){
  int numSites(4);
  auto basis = generateDefaultBasis(numSites);
  testBasis(basis);
}
