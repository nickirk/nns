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
#include "CoeffType.hpp"
#include "Nnw.hpp"

class Sampler{
public:
  Sampler(Hamiltonian const &H_, Basis const &fullBasis_, int numDets_, detType const &HF):H(H_),fullBasis(fullBasis_),
											     numDets(numDets_),cDet(HF){}
  virtual ~Sampler(){};
  // and functionalities: get a random coupled determinant
  detType getRandomConnection(detType const &startingPoint) const;
  virtual detType getDet() const{return cDet;}
  // This function only exists for compatibility with the markov implementation
  virtual void iterate(coeffType &cI, detType &dI, NeuralNetwork const &NNW) const
  	 	 {dI=getRandomConnection(cDet);cI=NNW.getCoeff(dI);}
  // set the starting point
  void setReference(detType const &start){cDet = start;}
  void setNumDets(int newNumDets){numDets = newNumDets;}
  int getNumDets()const {return numDets;}
private:
  // Hamiltonian
  Hamiltonian const &H;
  // and the corresponding basis
  Basis const &fullBasis;
protected:
  // Number of states to sample
  size_t numDets;
  // this is the current sample in terms of determinants
  detType cDet;
};

inline double dblAbs(double x){if(x>0) return x;return -x;}
void removeDuplicate(std::vector<detType> &list);
detType getRandomDeterminant(std::vector<int> const &spinConfig);
#endif
