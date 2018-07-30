/*
 * StochasticReconfiguration.cxx
 *
 *  Created on: May 2, 2018
 *      Author: guther
 */

#include "StochasticReconfiguration.hpp"
#include <iostream>
namespace networkVMC {

template <typename T>
StochasticReconfiguration<T>::~StochasticReconfiguration() {
}

//---------------------------------------------------------------------------------------------------//

// Conjugate descent (I ported this from Nnw.cxx and will not touch it)
template <typename T>
void StochasticReconfiguration<T>::update(T &w, T const &force,
		  State const &input){

	// first, get the input vector's coefficients
  std::size_t numDets = input.size();
  int numPars=w.size();
  Eigen::MatrixXcd OkOkp = Eigen::MatrixXcd::Zero(numPars, numPars);
  Eigen::MatrixXcd S = Eigen::MatrixXcd::Zero(numPars, numPars);
  Eigen::VectorXcd Ok = Eigen::VectorXcd::Zero(numPars);
  Eigen::VectorXcd dCdWk= Eigen::VectorXcd::Zero(numPars);
  // do the mapping inside for loop, private
	for(std::size_t i = 0; i<numDets; ++ i){
    dCdWk=NNW.getSRDeriv(input.det(i));
    //dCdWk=NNW.getDeriv(input.det(i));
    OkOkp += (dCdWk*dCdWk.adjoint()).adjoint();
    Ok += dCdWk;
	}
  OkOkp /= static_cast<double>(numDets);
  Ok /= static_cast<double>(numDets);
  S = OkOkp - (Ok*Ok.adjoint()).adjoint();
  double lambda = std::max(100*std::pow(0.999,iteration), 10.);
  std::cout << "StochasticReconfiguration.cxx: lambda=" << lambda << std::endl;
  Eigen::VectorXcd I = Eigen::VectorXcd::Ones(numPars);
  S+=I.asDiagonal()*lambda;
  w-=Solver<T>::learningRate*S.inverse()*force;
	// increase the iteration counter
	iteration += 1;
}
//template class StochasticReconfiguration<VecType>;
template class StochasticReconfiguration<VecCType>;
} /* namespace networkVMC */
