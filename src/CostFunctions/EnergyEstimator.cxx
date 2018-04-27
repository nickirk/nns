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

namespace networkVMC{

double EnergyEstimator::evaluate(State const &input) const{
  double energyVal{0.0};
  normalizerCoeff=0.0;
  std::complex<double> normalizerCoeffComplex(0.,0.);
  double Hij(0.);
  int numDets = input.size();
  for (int i=0; i < numDets; ++i){
    coeffType c_i=input.coeff(i);
    std::vector<coeffType> coupledC_j = input.coupledCoeffs(i);
    std::vector<detType> coupledDets = input.coupledDets(i);
    normalizerCoeffComplex += std::norm(c_i);
    //sign_i = (output_Cs[i]-0. < 1e-8)?-1:0;
    Hij = H(input.det(i), input.det(i));
    energyVal += std::real(std::conj(c_i) * c_i * Hij);
    for (size_t j=0; j < coupledC_j.size(); ++j){
      Hij = H(input.det(i), coupledDets[j]);
      energyVal += std::real(std::conj(c_i) * coupledC_j[j] * Hij);
    }
  }
  normalizerCoeff = std::real(normalizerCoeffComplex);
  energyVal /= normalizerCoeff;
  return energyVal;
}

std::vector<Eigen::VectorXd> EnergyEstimator::nabla(State const &input) const{
  energy = evaluate(input);
  int numDets = input.size();
  std::vector<Eigen::VectorXd> dEdC(numDets);
  for (int i=0; i < numDets; ++i){
    Eigen::VectorXd dEdC_i=Eigen::VectorXd::Zero(2);
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
    dEdC_i[0] = 2. * std::real( A*std::conj(1.));
    dEdC_i[1] = 2. * std::real( A*std::conj(std::complex<double>(0.0,1.0)));
    dEdC[i]=(dEdC_i);
  }
  return dEdC;
}

}
