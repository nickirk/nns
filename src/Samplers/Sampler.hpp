/*
 * Sampler.hpp
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */
#ifndef Sampler_DEFINED
#define Sampler_DEFINED

#include <vector>

#include "../Hamiltonian/ExcitationGenerators/ExcitationGenerator.hpp"
#include "../utilities/TypeDefine.hpp"
#include "../utilities/DeepCpyUniquePtr.hpp"
namespace networkVMC{

// Forward declaration to make compilation faster
class Hamiltonian;
class Basis;

// Base class for sampling, these objects take some input state and a Hamiltonian and generate
// lists of potentially relevant determinants
class Sampler{
  public:
  // explicit constructor
  Sampler(ExcitationGenerator const &eG_, int numDets_= 100);
  // construct the ExcitationGenerator implicitly from the Hamiltonian
  Sampler(Hamiltonian const &H_, detType const &HF_, int numDets_ = 100);
  Sampler(ExcitationGenerator const &eG_, detType const &HF, Basis const &fullBasis, 
      int numDets_ = 100);
  virtual ~Sampler(){};
  virtual Sampler* clone() const = 0;
  // This function is what samplers ought to do: Get a random determinant with some
  // coefficient
  // the only way to parallelize this is to pass the iteration count, too
  virtual void iterate(coeffType &cI, detType &dI, double& weight, int i)=0;

  // update the bias used for excitation generation in iterate()
  void updateBiases()const {excitGen->updateBiases();}

  // set the number of sampled dets
  virtual void setNumDets(int newNum){numDets=newNum;}
  // get the number of dets
  virtual int getNumDets() const {return numDets;}

  // for attributing CostFunctions
  virtual SamplerType type() const {return PreFetched;}
  friend void swap(Sampler &a, Sampler &b){
    std::swap(a.excitGen,b.excitGen);
    std::swap(a.numDets,b.numDets);
  }

protected:
  // and functionalities for implementation of sampling: get a random coupled determinant
  detType getRandomConnection(detType const &startingPoint, double &p) const;
  // get the probability of generating a when starting from b
  double getConnectionProb(detType const &source, detType const &target) const;

private:
  // excitation generator that does the sampling of determinants
  // it is owned by the sampler, so we can also create it in the
  // background, it does not need to be user-specified
  mutable DeepCpyUniquePtr<ExcitationGenerator> excitGen;
  // is mutable because the excitation generator changes its behaviour
  // over iterations, but this is not visible outisde
protected:
  int numDets;
  // this is the current sample in terms of determinants
  detType cDet;
};

inline double dblAbs(double x){if(x>0) return x;return -x;}

}
#endif
