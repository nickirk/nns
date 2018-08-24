/*
 * EnergyEstimator.cxx
 * based on EnergyCF.cxx
 *  Created on: Nov 08, 2017
 *      Author: Liao
 */

#include "EnergyEsPreFetched.hpp"

#include <iostream>
#include <vector>
#include <complex>

#include "../Hamiltonian/TwoBodyHamiltonian.hpp"
#include "../utilities/State.hpp"

namespace networkVMC{

template <typename F, typename coeffType>
coeffType EnergyEsPreFetched<F, coeffType>::evaluate(State<coeffType> const &input) const{
  //std::complex<double, Eigen::Matrix<double, Dynamic, 1> normalizerCoeffComplex(0.,0.);
  int numDets = input.size();
  //double energyVal = 0.0;
  coeffType energyVal = 0.0;
  normalizerCoeff=0.0;
  double norm = 0.0;
  {
    //#pragma omp parallel for reduction(+:energyVal, norm)
    for (int i=0; i < numDets; ++i){
      double Hij = 0.0;
      coeffType c_i=input.coeff(i);
      std::vector<coeffType> const &coupledC_j = input.coupledCoeffs(i);
      std::vector<detType> const &coupledDets = input.coupledDets(i);
      std::vector<double> const &coupledWeights = input.coupledWeights(i);
      norm += std::norm(c_i);
      //sign_i = (output_Cs[i]-0. < 1e-8)?-1:0;
      Hij = H(input.det(i), input.det(i));
      //energyVal += std::real(std::conj(c_i) * c_i * Hij);
      energyVal += std::conj(c_i) * c_i * Hij;
      for (size_t j=0; j < coupledC_j.size(); ++j){
        Hij = H(input.det(i), coupledDets[j]);
        //energyVal += std::real(std::conj(c_i) * coupledC_j[j] * Hij)/
        energyVal += std::conj(c_i) * coupledC_j[j] * Hij/
        		(coupledC_j.size()*coupledWeights[j]);
      }
    }
  }
  normalizerCoeff = norm;
  energyVal /= normalizerCoeff;
  return energyVal;
}

template <typename F, typename coeffType>
EnergyEsPreFetched<F, coeffType>::T EnergyEsPreFetched<F, coeffType>::nabla(State<coeffType> const &input) const{
  energy = evaluate(input);
  int numDets = input.size();
  int spaceSize = input.totalSize();
  T dEdC(spaceSize);

  #pragma omp parallel for
  for (int i=0; i < numDets; ++i){
    coeffType dEdC_i;
    coeffType A=0.;
    coeffType c_i = input.coeff(i);
    std::vector<coeffType> const &coupledC_j = input.coupledCoeffs(i);
    std::vector<detType> const &coupledDets = input.coupledDets(i);
    std::vector<double> const &coupledWeights = input.coupledWeights(i);
    dEdC[i] = c_i * (H(input.det(i), input.det(i)) - std::real(energy))
      /normalizerCoeff;
    int pos=input.locate(i);
    for (size_t j=0; j < coupledDets.size(); ++j){
      dEdC(numDets+pos+j) = coupledC_j[j] * H(input.det(i), coupledDets[j])
        /(coupledC_j.size()*coupledWeights[j])/normalizerCoeff;
      
    }
  }
  return dEdC;
}

template class EnergyEsPreFetched<double, double>;
template class EnergyEsPreFetched<std::complex<double>, std::complex<double>>;
}
