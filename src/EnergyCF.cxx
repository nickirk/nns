/*
 * EnergyCF.cxx
 *
 *  Created on: Oct 25, 2017
 *      Author: guther
 */

#include "EnergyCF.hpp"
#include <vector>
#include <complex>

double EnergyCF::evaluate(State const &input) const{
  double energyVal{0.0};
  normalizerCoeff=0.0;
  std::complex<double> normalizerCoeffComplex(0.,0.);
  double Hij(0.);
  int numDets = input.size();
  double real(0.), imag(0.);
  for (int i=0; i < numDets; ++i){
    real = input.getCoeff(i)[0];
    imag = input.getCoeff(i)[1];
    std::complex<double> c_i(real, imag);
    normalizerCoeffComplex += fabs (std::conj(c_i)  * c_i);
    //sign_i = (output_Cs[i]-0. < 1e-8)?-1:0;
    for (int j=0; j < numDets; ++j){
      real = input.getCoeff(j)[0];
      imag = input.getCoeff(j)[1];
      std::complex<double> c_j(real, imag);
      Hij = H(input.getDet(i), input.getDet(j));
      energyVal += std::real(std::conj(c_i) * c_j * Hij);
    }
  }
  //std::cout << "normE= " << normalizerCoeff << std::endl;
  normalizerCoeff = std::real(normalizerCoeffComplex);
  energyVal /= normalizerCoeff;
  return energyVal;
}

std::vector<coeffType> EnergyCF::nabla(State const &input) const{
  energy = evaluate(input);
  std::vector<Eigen::VectorXd> dEdC;
  int numDets = input.size();
  for (int i=0; i < numDets; ++i){
    Eigen::Vector2d dEdC_i=Eigen::Vector2d::Zero();
    std::complex<double> A(0.,0.);
    for (int j=0; j < numDets; ++j){
      std::complex<double> c_j(input.getCoeff(j)[0], input.getCoeff(j)[1]);
      A += c_j * H(input.getDet(i),
                        input.getDet(j));
    }
    std::complex<double> c_i(input.getCoeff(i)[0], input.getCoeff(i)[1]);
    A -=  energy * c_i;
    A /= normalizerCoeff;
    dEdC_i[0] = 2. * std::real( A*std::conj(1.));
    dEdC_i[1] = 2. * std::real( A*std::conj(std::complex<double>(0.0,1.0)));
    dEdC.push_back(dEdC_i);
  }
  return dEdC;
}

