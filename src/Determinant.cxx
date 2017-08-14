#include <vector>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include "Determinant.hpp"
#include "errors.hpp"
Determinant::Determinant(){
  size = 0;
  det.resize(size, 0);
}
Determinant::Determinant(int size_){
  size = size_;
  det.resize(size, 0); 
}

Determinant::Determinant(Determinant const &determinant_){
  size = determinant_.getSize();
  det.resize(size, 0); 
  std::vector<int> pos = determinant_.getOccupiedPositions();
  for (size_t i=0; i < pos.size(); i++){
    this->create(pos[i]);
  }
}

void Determinant::operator = (Determinant const &determinant_){
   size = determinant_.getSize();
   det.resize(size, 0);
   std::vector<int> pos(determinant_.getOccupiedPositions());
   for (int i=0; i<pos.size(); ++i){
     det[pos[i]]=1;
   }
}

void annihilate(detType &det, int pos){
  if (pos > det.size() || pos < 0){
    std::cerr << "Error! Create outside of range! "; 
    throw outOfRangeError(pos);
  }
  if (det[pos]){ 
    det[pos] = 0;
  }
  else{ 
    std::cerr << "Error! Cannot annihilate an unoccupied state! ";
    throw invalidAnnihilation(pos);
  }
}

void create(detType &det, int pos){
  if (pos > det.size() || pos < 0){
    std::cerr << "Error! Create outside of range! "; 
    throw outOfRangeError(pos);
  }
  if (det[pos]){
    std::cerr << "Error! Cannot create on an occupied state! ";
    throw invalidCreation(pos);
  }
  else{
    det[pos] = 1;
  }
}

std::vector<int> getOccupiedPositions(detType const &det) {
  std::vector<int> positions;
  for (int i=0; i < det.size(); i++){
    if (det[i]) positions.push_back(i);
  }
  return positions;
}
