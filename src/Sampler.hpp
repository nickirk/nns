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
#include "Nnw.hpp"

class Sampler{
public:
  Sampler(Hamiltonian const &H_, Basis const &fullBasis_, NeuralNetwork const &NNW_, int numDets_, detType const &HF):H(H_),fullBasis(fullBasis_),NNW(NNW_),
											     numDets(numDets_),cDet(std::vector<detType >(1,HF)){}
  Sampler(Hamiltonian const &H_, Basis const &fullBasis_, NeuralNetwork const &NNW_, int numDets_, std::vector<detType > const &reference):H(H_),
														  fullBasis(fullBasis_),
														  NNW(NNW_),
														  numDets(numDets_),
														  cDet(reference)
  {};
  // two functionalities: get a random coupled determinant and get an array of 
  // random coupled determinants
  void generateList(std::vector<detType > &list) const;
  detType getRandomConnection(detType const &startingPoint) const;
  void initialiseList(std::vector<detType> &list, std::vector<int> const& spinConfig);
  void diffuse(std::vector<detType> &list, std::vector<int> const& spinConfig);
  // for ab-initio: introduce an overload of generateList for ab-initio hamiltonians

  // set the starting point
  void setReference(std::vector<detType > const &list){cDet = list;}
  void setNumDets(int newNumDets){numDets = newNumDets;}
  int getNumDets()const {return numDets;}
private:
  // Hamiltonian
  Hamiltonian const &H;
  Basis const &fullBasis;
  NeuralNetwork const &NNW;
  int numDets;
  // this is the reference space in terms of determinants
  std::vector<detType > cDet;
};

inline double dblAbs(double x){if(x>0) return x;return -x;}
void removeDuplicate(std::vector<detType> &list);
detType getRandomDeterminant(std::vector<int> const &spinConfig);
#endif
