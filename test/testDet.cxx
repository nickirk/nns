#include <stdio.h>
#include <iostream>
#include "../src/Determinant.hpp"
#include "../src/Basis.hpp"
//Determinant.hpp/.cxx have been tested and passed.
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

  //create out of range
  //create(det1, 6);
  //annihilate out of range
  //annihilate(det1, 6);
  //create on occupied state
  //create(det1, 3);
  //annihilate on empty state
  annihilate(det1, 1);
  return 0;
}
