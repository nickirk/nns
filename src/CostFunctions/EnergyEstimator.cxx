/*
 * EnergyEstimator.cxx
 * based on EnergyCF.cxx
 *  Created on: Nov 08, 2017
 *      Author: Liao
 */

#include "EnergyEstimator.hpp"

#include <iostream>
#include <vector>
#include <complex>

#include "../Hamiltonian/Hamiltonian.hpp"
#include "../utilities/State.hpp"

namespace networkVMC{

double EnergyEstimator::evaluate(State const &input) const{
  //std::complex<double> normalizerCoeffComplex(0.,0.);
  int numDets = input.size();
  double energyVal{0.0};
  normalizerCoeff=0.0;
  double norm(0.);
  {
    #pragma omp parallel for reduction(+:energyVal, norm)
    for (int i=0; i < numDets; ++i){
      double Hij(0.);
      coeffType c_i=input.coeff(i);
      std::vector<coeffType> coupledC_j = input.coupledCoeffs(i);
      std::vector<detType> coupledDets = input.coupledDets(i);
      norm += std::norm(c_i);
      int tid = omp_get_thread_num();
      //sign_i = (output_Cs[i]-0. < 1e-8)?-1:0;
      Hij = H(input.det(i), input.det(i));
      energyVal += std::real(std::conj(c_i) * c_i * Hij);
      for (size_t j=0; j < coupledC_j.size(); ++j){
        Hij = H(input.det(i), coupledDets[j]);
        energyVal += std::real(std::conj(c_i) * coupledC_j[j] * Hij);
      }
    }
  }
  normalizerCoeff = norm;
  energyVal /= normalizerCoeff;
  return energyVal;
}

nablaType EnergyEstimator::nabla(State const &input) const{
  energy = evaluate(input);
  int numDets = input.size();
  std::vector<coeffType> dEdC(numDets, coeffType(0.,0.));

  #pragma omp parallel for
  for (int i=0; i < numDets; ++i){
    coeffType dEdC_i;
    std::complex<double> A(0.,0.);
    coeffType c_i = input.coeff(i);
    std::vector<coeffType> coupledC_j = input.coupledCoeffs(i);
    std::vector<detType> coupledDets = input.coupledDets(i);
    A += c_i * H(input.det(i), input.det(i));
    for (size_t j=0; j < coupledDets.size(); ++j){
      A += coupledC_j[j] * H(input.det(i),
                        coupledDets[j]);
    }
    A -=  energy * c_i;
    A /= normalizerCoeff;
    //dEdC_i.real(2. * std::real( A*std::conj(coeffType(1.,0.))));
    //dEdC_i.imag(2. * std::real( A*std::conj(coeffType(0.,1.))));
    //dEdC[i]=dEdC_i;
    // use the conjugate to make real and complex parametrization the same form
    // for real para, inside the para, use conj again to obtain the correct 
    // gradient. See Nnw.cxx backpropagate().
    dEdC[i]=(2.0 * std::conj(A));
  }
  return dEdC;
}

}
