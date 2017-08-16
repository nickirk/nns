/*
 * Sampler.hpp
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */
#ifndef Sampler_DEFINED
#define Sampler_DEFINED

#include <vector>
#include "Determinant.hpp"
#include "Basis.hpp"
#include "Hamiltonian.hpp"

class sampler{
public:
  sampler(Hamiltonian const &H_, detType const &HF, Basis const &fullBasis_, int numStates_):H(H_),fullBasis(fullBasis_),
											     numStates(numStates_),cDet(std::vector<detType >(1,HF)){}
  // two functionalities: get a random coupled determinant and get an array of 
  // random coupled determinants
  void generateList(std::vector<detType > &list)const;
  detType getRandomDeterminant(detType const &startingPoint) const;
  // for ab-initio: introduce an overload of generateList for ab-initio hamiltonians

  // set the starting point
  void setReference(std::vector<detType > const &list){cDet = list;}
  void setNumStates(int newNumStates){numStates = newNumStates;}
  int getNumStates()const {return numStates;}
private:
  // Hamiltonian
  Hamiltonian &H;
  Basis &fullBasis;
  int numStates;
  // this is the reference space in terms of determinants
  std::vector<detType > cDet;
};

inline double dblAbs(double x){if(x>0) return x;return -x;}

#endif
