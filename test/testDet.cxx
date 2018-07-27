#include <stdio.h>
#include <iostream>

#include "../src/HilbertSpace/Basis.hpp"
#include "../src/HilbertSpace/Determinant.hpp"
//Determinant.hpp/.cxx have been tested and passed.

using namespace networkVMC;
int main(){
  detType det(5,0);
  detType det1(det);
  std::cout << "size det1= " << det1.size() << std::endl;
  for (int i = 0; i < getOccupiedPositions(det1).size(); i++)
  {
    std::cout << getOccupiedPositions(det1)[i] << "," ; 
  }
  std::cout << std::endl;
  std::cout << "verbatimCast of det1= " << verbatimCast(det1) << std::endl;
  create(det1, 3);
  create(det1, 2);
  create(det1, 0);
  for (int i = 0; i < getOccupiedPositions(det1).size(); i++)
  {
    std::cout << "det1 =" << getOccupiedPositions(det1)[i] << "," ; 
  }
  std::cout << std::endl;
  std::cout << "verbatimCast of det1= " << verbatimCast(det1) << std::endl;
  //annihilate
  annihilate(det1, 2);
  for (int i = 0; i < getOccupiedPositions(det1).size(); i++)
  {
    std::cout << getOccupiedPositions(det1)[i] << "," ; 
  }
  std::cout << std::endl;

  // test fermi sign
  detType a(5,false);
  create(a,1);
  create(a,3);
  create(a,4);
  std::cout << "Sign 3->5: " << excitationSign(a,3,5) << std::endl;
  std::cout << "Sign 4->2: " << excitationSign(a,4,2) << std::endl;
  auto b = excite(a,3,5);
  std::cout << "Sign 5->3: " << excitationSign(b,5,3) << std::endl;
  return 0;
}
