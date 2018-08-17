/*
 * StochasticReconfiguration.cxx
 *
 *  Created on: May 2, 2018
 *      Author: guther and Liao
 */

#include "StochasticReconfiguration.hpp"
#include "../utilities/State.hpp"
#include <iostream>
namespace networkVMC {

template <typename T>
StochasticReconfiguration<T>::~StochasticReconfiguration() {
}

//---------------------------------------------------------------------------------------------------//

template <typename T>
void StochasticReconfiguration<T>::update(T &w, T const &force,
		  State const &input, SamplerType const &samplerType){

	// first, get the input vector's coefficients
  std::size_t numDets = input.size();
  int numPars=w.size();
  Eigen::MatrixXcd OkOkp = Eigen::MatrixXcd::Zero(numPars, numPars);
  Eigen::MatrixXcd S = Eigen::MatrixXcd::Zero(numPars, numPars);
  Eigen::VectorXcd Ok = Eigen::VectorXcd::Zero(numPars);
  Eigen::VectorXcd dCdWk= Eigen::VectorXcd::Zero(numPars);
  switch (samplerType){
    case Markov:
	  for(std::size_t i = 0; i<numDets; ++ i){
      // the derivatives differ in different sampling schemes
      dCdWk=NNW.getMarkovDeriv(input.det(i));
      OkOkp += (dCdWk*dCdWk.adjoint()).adjoint();
      Ok += dCdWk;
	  }
    OkOkp /= static_cast<double>(numDets);
    Ok /= static_cast<double>(numDets);
    break;
    case PreFetched:
	  for(std::size_t i = 0; i<numDets; ++ i){
      // the derivatives differ in different sampling schemes
      dCdWk=NNW.getDeriv(input.det(i));
      OkOkp += (dCdWk*dCdWk.adjoint()).adjoint();
      Ok += dCdWk;
	  }
    break;
    default:
    throw SamplerTypeDoesNotExist(samplerType);
  }

  S = OkOkp - (Ok*Ok.adjoint()).adjoint();
  // default hyperparameters tested with small Hubbard models,
  // no guarantee that they work for other systems. 
  // More tests should be run to observe the performace of these parameters
  double lambda = std::max(100*std::pow(0.99,iteration), 5.);
  std::cout << "StochasticReconfiguration.cxx: lambda=" << lambda << std::endl;
  Eigen::VectorXcd I = Eigen::VectorXcd::Ones(numPars);
  // Add \epsilon * I to S matrix to prevent ill inversion of the S matrix.
  S+=I.asDiagonal()*lambda;
  Eigen::MatrixXcd II = I.asDiagonal();
  std::cout << "StochasticReconfiguration.cxx: S=" << std::endl;
  std::cout << S << std::endl;
  w-=Solver<T>::learningRate*S.llt().solve(II)*force;
	// increase the iteration counter
	iteration += 1;
}
//template class StochasticReconfiguration<VecType>;
template class StochasticReconfiguration<VecCType>;
} /* namespace networkVMC */
