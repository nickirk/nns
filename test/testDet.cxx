#include <stdio.h>
#include <iostream>
#include "../src/Determinant.hpp"
#include "../src/Basis.hpp"

int main(){
  Determinant det(5);
  Determinant det1(det);
  std::cout << det1.getSize() << std::endl;
  for (int i = 0; i < det1.getOccupiedPositions().size(); i++)
  {
    std::cout << det1.getOccupiedPositions()[i] << "," ; 
  }
  std::cout << std::endl;
  std::cout << det1.intCast() << std::endl;
  det1.create(3);
  det1.create(2);
  det1.create(0);
  for (int i = 0; i < det1.getOccupiedPositions().size(); i++)
  {
    std::cout << det1.getOccupiedPositions()[i] << "," ; 
  }
  std::cout << std::endl;
  std::cout << det1.intCast() << std::endl;
  det1.create(4);
  for (int i = 0; i < det1.getOccupiedPositions().size(); i++)
  {
    std::cout << det1.getOccupiedPositions()[i] << "," ; 
  }
  std::cout << std::endl;
  std::cout << det1.intCast() << std::endl;
  det1.annihilate(2);
  for (int i = 0; i < det1.getOccupiedPositions().size(); i++)
  {
    std::cout << det1.getOccupiedPositions()[i] << "," ; 
  }
  std::cout << std::endl;
  std::cout << det1.intCast() << std::endl;
  
  return 0;
}
