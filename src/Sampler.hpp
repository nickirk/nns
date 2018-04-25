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
#include "SpinConfig.hpp"

// Base class for sampling, these objects take some input state and a Hamiltonian and generate
// lists of potentially relevant determinants
class Sampler{
public:
  Sampler(Hamiltonian const &H_, Basis const &fullBasis_, detType const &HF, int numDets_= 100):H(H_),numDets(numDets_),fullBasis(fullBasis_),
											     cDet(HF),sC(fullBasis_.getSpinConfig()){}
  virtual ~Sampler(){};
  // and functionalities: get a random coupled determinant
  detType getRandomConnection(detType const &startingPoint) const;
  virtual detType getDet() const{return cDet;};
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
  // set the number of sampled dets
  virtual void setNumDets(int newNum){numDets=newNum;}
  // get the number of dets
  virtual int getNumDets() const {return numDets;}
  Hamiltonian const& getH() {return H;}
private:
  // Hamiltonian
  Hamiltonian const &H;
protected:
  int numDets;
  // and the corresponding basis
  Basis const &fullBasis;
  // this is the current sample in terms of determinants
  mutable detType cDet;
  // the information on the number of electrons with a given spin
  SpinConfig sC;
};

inline double dblAbs(double x){if(x>0) return x;return -x;}
void removeDuplicate(std::vector<detType> &list);
detType getRandomDeterminant(SpinConfig const &spinConfig);
#endif
