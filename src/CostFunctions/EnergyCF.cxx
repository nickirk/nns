/*
 * EnergyCF.cxx
 *
 *  Created on: Oct 25, 2017
 *      Author: guther
 */

#include "EnergyCF.hpp"

#include <vector>
#include <complex>

namespace networkVMC{

// This gets the energy expectation value of the state 'input'
double EnergyCF::evaluate(std::vector<State> const &input) const{
  double energyVal{0.0};
  normalizerCoeff=0.0;
  std::complex<double> normalizerCoeffComplex(0.,0.);
  double Hij(0.);
  int numDets = input.size();
  for (int i=0; i < numDets; ++i){
	  // The entries of state are pairs of determinant/coefficient
    coeffType c_i=input[i].coeff;
    // also assign the denominator
    normalizerCoeffComplex += std::norm(c_i);
    //sign_i = (output_Cs[i]-0. < 1e-8)?-1:0;
    for (int j=0; j < numDets; ++j){
      coeffType c_j = input[j].coeff;
      Hij = H(input[i].det, input[j].det);
    // sum up the contributions to E
      energyVal += std::real(std::conj(c_i) * c_j * Hij);
    }
  }
  normalizerCoeff = std::real(normalizerCoeffComplex);
  energyVal /= normalizerCoeff;
  return energyVal;
}

// Here, we get the derivative of the energy with respect to the coefficients of input
std::vector<Eigen::VectorXd> EnergyCF::nabla(std::vector<State> const &input) const{
// To get the derivative, we also need the energy, so first evaluate
// This also assigns the normalizer
  energy = evaluate(input);
  std::vector<Eigen::VectorXd> dEdC;
  int numDets = input.size();
  // works similar to computing the energy
  for (int i=0; i < numDets; ++i){
    Eigen::VectorXd dEdC_i=Eigen::VectorXd::Zero(2);
    std::complex<double> A(0.,0.);
    for (int j=0; j < numDets; ++j){
      coeffType c_j=input[j].coeff;
      // We add up the contributions to the derivative
      A += c_j * H(input[i].det,
                        input[j].det);
    }
    coeffType c_i = input[i].coeff;
    A -=  energy * c_i;
    // re-use the normalizer from evaluate
    A /= normalizerCoeff;
    // And assign to the output
    dEdC_i[0] = 2. * std::real( A*std::conj(1.));
    dEdC_i[1] = 2. * std::real( A*std::conj(std::complex<double>(0.0,1.0)));
    dEdC.push_back(dEdC_i);
  }
  return dEdC;
}

}

