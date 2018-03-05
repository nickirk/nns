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

double EnergyEsMarkov::evaluate(std::vector<State> const &input) const{
  double energyVal{0.0};
  normalizerCoeff=0.0;
  double Hij(0.);
  int numDets = input.size();
  for (int i=0; i < numDets; ++i){
    coeffType c_i=input[i].coeff;
    std::vector<coeffType> coupledC_j = input[i].coupledCoeffs;
    std::vector<detType> coupledDets = input[i].coupledDets;
    Hij = H(input[i].det, input[i].det);
    energyVal += std::real(c_i * Hij/c_i);
    for (size_t j=0; j < coupledC_j.size(); ++j){
      Hij = H(input[i].det, coupledDets[j]);
      energyVal += std::real(coupledC_j[j] * Hij / c_i );
    }
  }
  energyVal /= numDets;
  return energyVal;
}

std::vector<Eigen::VectorXd> EnergyEsMarkov::nabla(std::vector<State > const &input) const{
  energy = evaluate(input);
  int numDets = input.size();
  std::vector<Eigen::VectorXd> dEdC;
  Eigen::VectorXd dEdC_i=Eigen::VectorXd::Zero(1);
  for (int i=0; i < numDets; ++i){
    //coeffType c_i = input.getCoeff(i);
    if (i == 0 || input[i].det != input[i-1].det){
      double A(0.);
      A = (H(input[i].det, input[i].det) - energy ); 
      dEdC_i(0) = 2. * (A);
      dEdC.push_back(dEdC_i);
      for (size_t j=0; j < input[i].coupledDets.size(); ++j){
        A = H(input[i].det,input[i].coupledDets[j]);
        dEdC_i(0) = 2. * A;
        dEdC.push_back(dEdC_i);
      }
    }
  }
  return dEdC;
}

