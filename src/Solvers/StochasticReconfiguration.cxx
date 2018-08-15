/*
 * StochasticReconfiguration.cxx
 *
 *  Created on: May 2, 2018
 *      Author: guther and Liao
 */

#include "StochasticReconfiguration.hpp"
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
      //dCdWk=NNW.getDeriv(input.det(i));
      OkOkp += (dCdWk*dCdWk.adjoint()).adjoint();
      Ok += dCdWk;
	  }
    OkOkp /= static_cast<double>(numDets);
    Ok /= static_cast<double>(numDets);
    //std::cout << "StochasticReconfiguration.cxx: OkOkp=" << std::endl;
    //std::cout << OkOkp << std::endl;
    break;
    case PreFetched:
	  for(std::size_t i = 0; i<numDets; ++ i){
      // the derivatives differ in different sampling schemes
      dCdWk=NNW.getDeriv(input.det(i));

      OkOkp += (dCdWk*dCdWk.adjoint()).adjoint();
      Ok += dCdWk;
	  }
    OkOkp /= static_cast<double>(numDets);
    Ok /= static_cast<double>(numDets);
    break;
    default:
    throw SamplerTypeDoesNotExist(samplerType);
  }

  S = OkOkp - (Ok*Ok.adjoint()).adjoint();
  // default hyperparameters tested with small Hubbard models,
  // no guarantee that they work for other systems. 
  // More tests should be run to observe the performace of these parameters
  double lambda = std::max(100*std::pow(0.9,iteration), 5.);
  Eigen::VectorXcd I = Eigen::VectorXcd::Ones(numPars);
  //std::cout << "I=" << I  << std::endl; 
  // Add \epsilon * I to S matrix to prevent ill inversion of the S matrix.
  double error=1.;
  Eigen::VectorXcd x;
  Eigen::ConjugateGradient<Eigen::MatrixXcd, Eigen::Lower|Eigen::Upper> cg;
  cg.setTolerance(1e-16);
  cg.setMaxIterations(10000);
  while (error>1e-16){
    //std::cout << "StochasticReconfiguration.cxx: lambda=" << lambda << std::endl;
    S+=I.asDiagonal()*lambda;
    //S+=S.diagonal().asDiagonal()*lambda;
    cg.compute(S);
    x = cg.solve(force);
    error = cg.error();
    lambda += 0.002;
    //std::cout << "StochasticReconfiguration.cxx: #iterations:     " 
    //  << cg.iterations() << std::endl;
    //std::cout << "StochasticReconfiguration.cxx: estimated error: " 
    //  << cg.error()      << std::endl;
  }
  x /= std::sqrt(x.dot(S * x).real()); 

  //std::cout << "StochasticReconfiguration.cxx: S=" << std::endl;
  //std::cout << S << std::endl;
  //S+=S.diagonal().asDiagonal()*lambda;
  //Eigen::MatrixXcd II = I.asDiagonal();
  //Eigen::MatrixXcd residual = S*(S.llt().solve(II)).adjoint();
  //std::cout << "StochasticReconfiguration.cxx: residual=" << std::endl;
  //std::cout << residual << std::endl;
  //w-=Solver<T>::learningRate*S.llt().solve(II)*force;
  //internalSolver.setLearningRate(Solver<T>::learningRate);
  //internalSolver.update(w, x, input, samplerType);
  w-=Solver<T>::learningRate*x;
	// increase the iteration counter
	iteration += 1;
}
//template class StochasticReconfiguration<VecType>;
template class StochasticReconfiguration<VecCType>;
} /* namespace networkVMC */
