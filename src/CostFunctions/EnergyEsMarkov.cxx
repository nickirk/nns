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

#include "../Hamiltonian/Hamiltonian.hpp"
#include "../utilities/State.hpp"

namespace networkVMC{

template <typename F, typename coeffType>
coeffType EnergyEsMarkov<F, coeffType>::evaluate(State<coeffType> const &input) const{
  coeffType energyVal=0.;
  //normalizerCoeff=input.getTotalWeights();
  int numDets = input.size();
  //#pragma omp parallel for reduction(+:energyVal)
  for (int i=0; i < numDets; ++i){
    double Hij(0.);
    coeffType c_i=input.coeff(i);
    std::vector<coeffType> coupledC_j = input.coupledCoeffs(i);
    std::vector<detType> coupledDets = input.coupledDets(i);
    std::vector<double> coupledWeights = input.coupledWeights(i);
    size_t coupledSize = coupledC_j.size();
    Hij = H(input.det(i), input.det(i)) * input.weight(i);
    energyVal += Hij;
    for (size_t j=0; j < coupledSize; ++j){
      // weight it with |T_i|^2, weight(i)
      Hij = H(input.det(i), coupledDets[j]) * input.weight(i);
      energyVal += coupledC_j[j]* Hij / (c_i * (coupledWeights[j] * coupledSize));
      //std::cout << "EnergyEsMarkov.cxx: cj/ci =" << coupledC_j[j]/c_i << std::endl;
      //std::cout << "EnergyEsMarkov.cxx: Hij =" << Hij << std::endl;
    }
  }
  // devide by total weights
  std::cout << "EnergyEsMarkov.cxx: totalweights=" << input.getTotalWeights() << std::endl;
  energyVal /= input.getTotalWeights();
  return energyVal;
}

template <typename F, typename coeffType>
EnergyEsMarkov<F, coeffType>::T EnergyEsMarkov<F, coeffType>::nabla(State<coeffType> const &input) const{
  coeffType energyM = evaluate(input);
  energy = energyM;
  int numDets = input.size();
  // spaceSize = size of sampled dets and their coupled ones
  int spaceSize = input.totalSize();
  // For each n, it has the same number of connected Dets as others.!!!
  T dEdC=T::Zero(spaceSize);
  //not thread safe
  //assume we know the whole space size, reserve space
  #pragma omp parallel for
  for (int i=0; i < numDets; ++i){
    coeffType c_i = input.coeff(i);
    // put all the weighting step here instead of inside of RBM
    F dEdCtmp = (H(input.det(i), input.det(i))) * input.weight(i);
    // divide by  totalWeights
    dEdC(i) = dEdCtmp / input.getTotalWeights();
    std::vector<coeffType> coupledC_j = input.coupledCoeffs(i);
    std::vector<detType> coupledDets = input.coupledDets(i);
    std::vector<double> coupledWeights = input.coupledWeights(i);
    size_t coupledSize = coupledDets.size();
    int pos=input.locate(i);
    for (size_t j=0; j < coupledSize; ++j){
      if (input.det(i)==coupledDets[j]){
        std::cout << "EnergyEsMarkov.cxx: Stop!" << std::endl;
        abort;
      }
      // weight
      dEdCtmp = H(input.det(i),coupledDets[j])* coupledC_j[j] * input.weight(i)
				/(c_i * (coupledWeights[j] * coupledSize ));
      // divided by totalWeights
      dEdC(numDets+pos+j)=dEdCtmp / input.getTotalWeights();
    }
   }
  // dEdC is stored as folllowing:
  // |c_0|c_1|...|c_{numDets-1}|c_0^{0}|c_0^{1}|...|c_0^{coupledSize-1}|c_1^{0}
  // |c..|c_{numDets-1}^{0}|...|c_{numDets-1}^{coupledSize-1}|
  return dEdC;
}

template class EnergyEsMarkov<double, double>;
template class EnergyEsMarkov<std::complex<double>, std::complex<double>>;
}

