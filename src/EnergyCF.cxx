/*
 * EnergyCF.cxx
 *
 *  Created on: Oct 25, 2017
 *      Author: guther
 */

#include "EnergyCF.hpp"
#include <vector>
#include <complex>

double EnergyCF::evaluate(std::vector<State> const &input) const{
  double energyVal{0.0};
  normalizerCoeff=0.0;
  std::complex<double> normalizerCoeffComplex(0.,0.);
  double Hij(0.);
  int numDets = input.size();
  for (int i=0; i < numDets; ++i){
    coeffType c_i=input[i].coeff;
    normalizerCoeffComplex += std::norm(c_i);
    //sign_i = (output_Cs[i]-0. < 1e-8)?-1:0;
    for (int j=0; j < numDets; ++j){
      coeffType c_j = input[j].coeff;
      Hij = H(input[i].det, input[j].det);
      energyVal += std::real(std::conj(c_i) * c_j * Hij);
    }
  }
  //std::cout << "normE= " << normalizerCoeff << std::endl;
  normalizerCoeff = std::real(normalizerCoeffComplex);
  energyVal /= normalizerCoeff;
  return energyVal;
}

std::vector<Eigen::VectorXd> EnergyCF::nabla(std::vector<State> const &input) const{
  energy = evaluate(input);
  std::vector<Eigen::VectorXd> dEdC;
  int numDets = input.size();
  for (int i=0; i < numDets; ++i){
    Eigen::VectorXd dEdC_i=Eigen::VectorXd::Zero(2);
    std::complex<double> A(0.,0.);
    for (int j=0; j < numDets; ++j){
      coeffType c_j=input[j].coeff;
      A += c_j * H(input[i].det,
                        input[j].det);
    }
    coeffType c_i = input[i].coeff;
    A -=  energy * c_i;
    A /= normalizerCoeff;
    dEdC_i[0] = 2. * std::real( A*std::conj(1.));
    dEdC_i[1] = 2. * std::real( A*std::conj(std::complex<double>(0.0,1.0)));
    dEdC.push_back(dEdC_i);
  }
  return dEdC;
}

