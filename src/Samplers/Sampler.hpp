/*
 * Sampler.hpp
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */
#ifndef Sampler_DEFINED
#define Sampler_DEFINED

#include <vector>

#include "../HilbertSpace/Basis.hpp"
#include "../HilbertSpace/Determinant.hpp"
#include "../utilities/SpinConfig.hpp"
#include "../utilities/TypeDefine.hpp"
#include "../Hamiltonian/Hamiltonian.hpp"

namespace networkVMC{

// Base class for sampling, these objects take some input state and a Hamiltonian and generate
// lists of potentially relevant determinants
class Sampler{
public:
  Sampler(Hamiltonian const &H_, Basis const &fullBasis_, detType const &HF, int numDets_= 100):
	  H(H_),numDets(numDets_),fullBasis(fullBasis_),cDet(HF){}
  virtual ~Sampler(){};
  // and functionalities: get a random coupled determinant
  detType getRandomConnection(detType const &startingPoint) const;
  // This function is what samplers ought to do: Get a random determinant with some
  // coefficient
  virtual void iterate(coeffType &cI, detType &dI) const=0;

  // return either the current det or possibly some stored det (depending on implementation)
  virtual detType getDet() const{return cDet;};
  virtual detType getDet(int i) const{return cDet;};

  // Setters for various properties
  // set the starting point
  void setReference(detType const &start){cDet = start;}
  // set the number of sampled dets
  virtual void setNumDets(int newNum){numDets=newNum;}
  // get the number of dets
  virtual int getNumDets() const {return numDets;}
  Hamiltonian const& getH() const{return H;}
private:
  // Hamiltonian is needed because it defines what is connected to a given determinant
  // and only connections are sampled
  Hamiltonian const &H;
protected:
  int numDets;
  // and the corresponding basis including the information on the number of electrons
  // with a given spin
  Basis const &fullBasis;
  // this is the current sample in terms of determinants
  mutable detType cDet;
};

inline double dblAbs(double x){if(x>0) return x;return -x;}
void removeDuplicate(std::vector<detType> &list);
detType getRandomDeterminant(SpinConfig const &spinConfig);

}
#endif
