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
#include "TypeDefine.hpp"
#include "Nnw.hpp"

class Sampler{
public:
  Sampler(Hamiltonian const &H_, Basis const &fullBasis_, int numDets_, detType const &HF):H(H_),fullBasis(fullBasis_),
											     numDets(numDets_),cDet(HF){}
  virtual ~Sampler(){};
  // and functionalities: get a random coupled determinant
  detType getRandomConnection(detType const &startingPoint) const;
  virtual detType getDet() const{return cDet;}
  virtual detType getDet(int i) const{return cDet;};
  // This function only exists for compatibility with the markov implementation
  virtual void iterate(coeffType &cI, detType &dI) const{
	  dI=getRandomConnection(cDet);
	  // WARNING: The default sampler does not take the coefficients into account
	  // IT HENCE CANNOT NOT OUTPUT THEM, they default to 0 and have to be set independently
	  cI = coeffType();
  }
  // set the starting point
  void setReference(detType const &start){cDet = start;}
  void setNumDets(int newNumDets){numDets = newNumDets;}
  virtual int getNumDets()const {return numDets;}
private:
  // Hamiltonian
  Hamiltonian const &H;
protected:
  // and the corresponding basis
  Basis const &fullBasis;
  // Number of states to sample
  size_t numDets;
  // this is the current sample in terms of determinants
  mutable detType cDet;
};

inline double dblAbs(double x){if(x>0) return x;return -x;}
void removeDuplicate(std::vector<detType> &list);
detType getRandomDeterminant(std::vector<int> const &spinConfig);
#endif
