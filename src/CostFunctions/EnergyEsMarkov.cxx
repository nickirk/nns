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
  double Hij(0.);
  int numDets = input.size();
  for (int i=0; i < numDets; ++i){
    coeffType c_i=input.coeff(i);
    std::vector<coeffType> coupledC_j = input.coupledCoeffs(i);
    std::vector<detType> coupledDets = input.coupledDets(i);
    Hij = H(input.det(i), input.det(i));
    energyVal += std::real(c_i * Hij/c_i);
    for (size_t j=0; j < coupledC_j.size(); ++j){
      Hij = H(input.det(i), coupledDets[j]);
      energyVal += std::real(coupledC_j[j] * Hij / c_i );
    }
  }
  energyVal /= numDets;
  return energyVal;
}

std::vector<Eigen::VectorXd> EnergyEsMarkov::nabla(State const &input) const{
  energy = evaluate(input);
  int numDets = input.size();
  std::vector<Eigen::VectorXd> dEdC;
  Eigen::VectorXd dEdC_i=Eigen::VectorXd::Zero(1);
  for (int i=0; i < numDets; ++i){
    //coeffType c_i = input.getCoeff(i);
    if (i == 0 || input.det(i) != input.det(i-1)){
      double A(0.);
      A = (H(input.det(i), input.det(i)) - energy );
      dEdC_i(0) = 2. * (A);
      dEdC.push_back(dEdC_i);
      for (size_t j=0; j < input.coupledDets(i).size(); ++j){
        A = H(input.det(i),input.coupledDets(i)[j]);
        dEdC_i(0) = 2. * A;
        dEdC.push_back(dEdC_i);
      }
    }
  }
  return dEdC;
}

}

