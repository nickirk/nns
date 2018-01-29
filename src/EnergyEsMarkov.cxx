/*
 * EnergyEsMarkov.cxx
 * based on EnergyCF.cxx
 *  Created on: Jan 29, 2018
 *      Author: Liao
 */

#include "EnergyEsMarkov.hpp"
#include <vector>
#include <complex>

double EnergyEsMarkov::evaluate(State const &input) const{
  double energyVal{0.0};
  normalizerCoeff=0.0;
  double Hij(0.);
  int numDets = input.size();
  for (int i=0; i < numDets; ++i){
    coeffType c_i=input.getCoeff(i);
    std::vector<coeffType> coupledC_j = input.getCoupledCoeffs(i);
    std::vector<detType> coupledDets = input.getCoupledDets(i);
    Hij = H(input.getDet(i), input.getDet(i));
    energyVal += std::real(c_i * Hij/c_i);
    for (size_t j=0; j < coupledC_j.size(); ++j){
      Hij = H(input.getDet(i), coupledDets[j]);
      energyVal += std::real(coupledC_j[j] * Hij / c_i );
    }
  }
  energyVal /= numDets;
  return energyVal;
}

std::vector<Eigen::VectorXd> EnergyEsMarkov::nabla(State const &input) const{
  energy = evaluate(input);
  int numDets = input.size();
  std::vector<Eigen::VectorXd> dEdC(numDets);
  for (int i=0; i < numDets; ++i){
    coeffType c_i = input.getCoeff(i);
    Eigen::VectorXd dEdC_i=Eigen::VectorXd::Zero(2);
    std::complex<double> A(0.,0.);
    std::vector<coeffType> coupledC_j = input.getCoupledCoeffs(i);
    std::vector<detType> coupledDets = input.getCoupledDets(i);
    A += H(input.getDet(i), input.getDet(i)) * c_i; 
    for (size_t j=0; j < coupledDets.size(); ++j){
      A += H(input.getDet(i),coupledDets[j]) * coupledC_j[j];
    }
    A /= std::norm(c_i);
    A -= energy/std::conj(c_i);
    A /= numDets;
    dEdC_i[0] = 2. * std::real( A);
    dEdC_i[1] = 2. * std::imag( A);
    dEdC[i]=(dEdC_i);
  }
  return dEdC;
}

