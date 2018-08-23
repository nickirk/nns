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

template <typename F, typename coeffType>
StochasticReconfiguration<F, coeffType>::~StochasticReconfiguration() {
}

//---------------------------------------------------------------------------------------------------//

template <typename F, typename coeffType>
void StochasticReconfiguration<F, coeffType>::update(T &w, T const &force,
		  State<coeffType> const &input, SamplerType const &samplerType){

	// first, get the input vector's coefficients
  std::size_t numDets = input.size();
  int numPars=w.size();
  Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> OkOkp = Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic>::Zero(numPars, numPars);
  Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> S = Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic>::Zero(numPars, numPars);
  T Ok = T::Zero(numPars);
  T dCdWk= T::Zero(numPars);
  switch (samplerType){
    case Markov:
	  for(std::size_t i = 0; i<numDets; ++ i){
      // the derivatives differ in different sampling schemes
      dCdWk=NNW.getMarkovDeriv(input.det(i));
      // weight OkOk and Ok
      OkOkp += (dCdWk*dCdWk.adjoint()).adjoint() * input.weight(i);
      Ok += dCdWk * input.weight(i);
	  }
    //OkOkp /= static_cast<double>(numDets);
    //Ok /= static_cast<double>(numDets);
    OkOkp /= input.getTotalWeights();
    Ok /= input.getTotalWeights();
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
  double lambda = std::max(100*std::pow(0.9,iteration), 10.);
  T I = T::Ones(numPars);
  //std::cout << "I=" << I  << std::endl; 
  // Add \epsilon * I to S matrix to prevent ill inversion of the S matrix.
  double error=1.;
  T x;
  Eigen::ConjugateGradient<Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Lower|Eigen::Upper> cg;
  //Eigen::ConjugateGradient<Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Lower|Eigen::Upper> cg;
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
  x /= std::sqrt(std::real(x.dot(S * x)));

  //std::cout << "StochasticReconfiguration.cxx: S=" << std::endl;
  //std::cout << S << std::endl;
  //S+=S.diagonal().asDiagonal()*lambda;
  //Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> II = I.asDiagonal();
  //Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> residual = S*(S.llt().solve(II)).adjoint();
  //std::cout << "StochasticReconfiguration.cxx: residual=" << std::endl;
  //std::cout << residual << std::endl;
  //w-=Solver<F, coeffType>::learningRate*S.llt().solve(II)*force;
  //internalSolver.setLearningRate(Solver<F, coeffType>::learningRate);
  //internalSolver.update(w, x, input, samplerType);
  w-=Solver<F, coeffType>::learningRate*x;
	// increase the iteration counter
	iteration += 1;
}
template class StochasticReconfiguration<double, double>;

template class StochasticReconfiguration<std::complex<double>, std::complex<double>>;
} /* namespace networkVMC */
