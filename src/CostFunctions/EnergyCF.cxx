/*
 * EnergyCF.cxx
 *
 *  Created on: Oct 25, 2017
 *      Author: guther
 */

#include "EnergyCF.hpp"
#include "../utilities/State.hpp"
#include <vector>
#include <complex>
#include "../Hamiltonian/Hamiltonian.hpp"

namespace networkVMC{

// This gets the energy expectation value of the state 'input'
template <typename F, typename coeffType>
coeffType EnergyCF<F, coeffType>::evaluate(State<coeffType> const &input) const{
  coeffType energyVal{0.0};
  normalizerCoeff=0.0;
  std::complex<double> normalizerCoeffComplex(0.,0.);
  double Hij(0.);
  int numDets = input.size();
  for (int i=0; i < numDets; ++i){
	  // The entries of state are pairs of determinant/coefficient
    coeffType c_i=input.coeff(i);
    // also assign the denominator
    normalizerCoeffComplex += std::norm(c_i);
    //sign_i = (output_Cs[i]-0. < 1e-8)?-1:0;
    for (int j=0; j < numDets; ++j){
      coeffType c_j = input.coeff(j);
      Hij = H(input.det(j), input.det(i));
    // sum up the contributions to E
      energyVal += std::real(std::conj(c_i) * c_j * Hij);
    }
  }
  normalizerCoeff = std::real(normalizerCoeffComplex);
  energyVal /= normalizerCoeff;
  return energyVal;
}

// Here, we get the derivative of the energy with respect to the coefficients of input
template <typename F, typename coeffType>
EnergyCF<F, coeffType>::T EnergyCF<F, coeffType>::nabla(State<coeffType> const &input) const{
// To get the derivative, we also need the energy, so first evaluate
// This also assigns the normalizer
  energy = evaluate(input);
  T dEdC;
  int numDets = input.size();
  // works similar to computing the energy
  for (int i=0; i < numDets; ++i){
    F dEdC_i;
    coeffType A=0.;
    for (int j=0; j < numDets; ++j){
      coeffType c_j=input.coeff(j);
      // We add up the contributions to the derivative
      A += c_j * H(input.det(i),
                        input.det(j));
    }
    coeffType c_i = input.coeff(i);
    A -=  energy * c_i;
    // re-use the normalizer from evaluate
    A /= normalizerCoeff;
    // And assign to the output
    // the way how it should implement depends on the feedforward nnw architecture, leave out for future when
    // this cost function is really needed.
    //dEdC_i.real(2. * std::real( A*std::conj(std::complex<double>(1.,0.))));
    //dEdC_i.imag(2. * std::real( A*std::conj(std::complex<double>(0.0,1.0))));
    dEdC(i)=(dEdC_i);
  }
  return dEdC;
}
template class EnergyCF<double, double>;
template class EnergyCF<std::complex<double>, std::complex<double>>;
}

