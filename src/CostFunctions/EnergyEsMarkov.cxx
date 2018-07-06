/*
 * EnergyEsMarkov.cxx
 * based on EnergyCF.cxx
 *  Created on: Jan 29, 2018
 *      Author: Liao
 */

#include "EnergyEsMarkov.hpp"

#include <vector>
#include <complex>
#include <iostream>

#include "../utilities/State.hpp"
#include "../Hamiltonian/Hamiltonian.hpp"

namespace networkVMC{

double EnergyEsMarkov::evaluate(State const &input) const{
  double energyVal{0.0};
  normalizerCoeff=0.0;
  int numDets = input.size();
  #pragma omp parallel for reduction(+:energyVal)
  for (int i=0; i < numDets; ++i){
    double Hij(0.);
    coeffType c_i=input.coeff(i);
    std::vector<coeffType> coupledC_j = input.coupledCoeffs(i);
    std::vector<detType> coupledDets = input.coupledDets(i);
    std::vector<double> coupledWeights = input.coupledWeights(i);
    Hij = H(input.det(i), input.det(i));
    energyVal += std::real(std::conj(c_i) * Hij/std::conj(c_i));
    //std::cout << "ci=" << c_i << std::endl;
    for (size_t j=0; j < coupledC_j.size(); ++j){
        // don't forget to unbias using the Pgen. TODO
      Hij = H(input.det(i), coupledDets[j]);
      energyVal += std::real(std::conj(coupledC_j[j])* Hij / std::conj(c_i)/
          coupledC_j.size()/coupledWeights[j]);
    }
  }
  energyVal /= numDets;
  return energyVal;
}

nablaType EnergyEsMarkov::nabla(State const &input) const{
  energy = evaluate(input);
  int numDets = input.size();
  // spaceSize = size of sampled dets and their coupled ones
  int spaceSize = input.totalSize();

  nablaType dEdC(spaceSize, coeffType(0.,0.)) ;
  //not thread safe
  //assume we know the whole space size, reserve space
  #pragma omp parallel for
  for (int i=0; i < numDets; ++i){
    coeffType dEdCtmp;
    coeffType c_i = input.coeff(i);
//  put all the weighting step here instead of inside of RBM
    dEdCtmp = (H(input.det(i), input.det(i)) - energy)/std::conj(c_i);
    // add weights
    dEdC[i] = dEdCtmp / numDets;
    std::vector<coeffType> coupledC_j = input.coupledCoeffs(i);
    std::vector<detType> coupledDets = input.coupledDets(i);
    std::vector<double> coupledWeights = input.coupledWeights(i);
    std::cout << "EnergyEsMarkov.cxx: c_i=" << c_i << std::endl;
    for (size_t j=0; j < coupledDets.size(); ++j){
      // don't forget to unbias using the Pgen. TODO
      // in the excitgen, the weight should be updated with pgen
      dEdCtmp = H(input.det(i),coupledDets[j])/std::conj(c_i)/ coupledC_j.size()/coupledWeights[j];
      //std::cout << "EnergyEsMarkov.cxx: coupledC_j=" << coupledC_j[j] << std::endl;
      //std::cout << "EnergyEsMarkov.cxx: coupledC_j size=" << coupledC_j.size() << std::endl;
      //std::cout << "EnergyEsMarkov.cxx: coupledWeights j=" << coupledWeights[j] << std::endl;
      // unbias with numCoupledDets and Pgen
      //dEdC[i+j+1]=dEdCtmp / coupledDets.size() / numDets;
      dEdC[i+j+1]=dEdCtmp / numDets;
    }
  }
  return dEdC;
}

}

