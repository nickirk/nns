#include <vector>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include "Determinant.hpp"
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

void Determinant::annihilate(int pos){
  if (pos > size || pos < 0){
    std::cerr << "Error! Create outside of range! "; 
    std::abort();
  }
  if (det[pos]){ 
    det[pos] = 0;
  }
  else{ 
    std::cerr << "Error! Cannot annihilate an unoccupied state! ";
    std::abort();
  }
}

void Determinant::create(int pos){
  if (pos > size || pos < 0){
    std::cerr << "Error! Create outside of range! "; 
    std::abort();
  }
  if (det[pos]){
    std::cerr << "Error! Cannot create on an occupied state! ";
    std::abort();
  }
  else{
    det[pos] = 1;
  }
}

std::vector<int> Determinant::getOccupiedPositions() const {
  std::vector<int> positions;
  for (int i=0; i < size; i++){
    if (det[i]) positions.push_back(i);
  }
  return positions;
}

int Determinant::getSize() const {return size;}

int Determinant::intCast() const{
  int out=0;
  for(int i=0;i<size;++i){
    if(det[i]){
      out+=pow(2,i);
    }
  }
  return out;
}
