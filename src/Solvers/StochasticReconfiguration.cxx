/*
 * StochasticReconfiguration.cxx
 *
 *  Created on: May 2, 2018
 *      Author: guther and Liao
 */

#include "StochasticReconfiguration.hpp"
#include "../utilities/State.hpp"
#include "../Network/Parametrization.hpp"
#include <Eigen/IterativeLinearSolvers>
namespace networkVMC {

StochasticReconfiguration::~StochasticReconfiguration() {
}

//---------------------------------------------------------------------------------------------------//

void StochasticReconfiguration::update(paraVector &w, paraVector const &force,
		  State const &input, SamplerType const &samplerType){

	// first, get the input vector's coefficients
  std::size_t numDets = input.size();
  int numPars=w.size();
  Eigen::Matrix<paraType, Eigen::Dynamic, Eigen::Dynamic> OkOkp = Eigen::Matrix<paraType, Eigen::Dynamic, Eigen::Dynamic>::Zero(numPars, numPars);
  Eigen::Matrix<paraType, Eigen::Dynamic, Eigen::Dynamic> S = Eigen::Matrix<paraType, Eigen::Dynamic, Eigen::Dynamic>::Zero(numPars, numPars);
  paraVector Ok = paraVector::Zero(numPars);
  paraVector dCdWk= paraVector::Zero(numPars);
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
    throw errors::SamplerTypeDoesNotExist(samplerType);
  }

  S = OkOkp - (Ok*Ok.adjoint()).adjoint();
  // default hyperparameters tested with small Hubbard models,
  // no guarantee that they work for other systems. 
  // More tests should be run to observe the performace of these parameters
  double lambda = std::max(100*std::pow(0.9,iteration), 10.);
  paraVector I = paraVector::Ones(numPars);
  // Add \epsilon * I to S matrix to prevent ill inversion of the S matrix.
  double error=1.;
  paraVector x;
  Eigen::ConjugateGradient<Eigen::Matrix<paraType, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Lower|Eigen::Upper> cg;
  cg.setTolerance(1e-16);
  cg.setMaxIterations(10000);
  while (error>1e-16){
    S+=I.asDiagonal()*lambda;
    cg.compute(S);
    x = cg.solve(force);
    error = cg.error();
    lambda += 0.002;
  }
  x /= std::sqrt(std::real(x.dot(S * x)));

  w-=Solver::learningRate*x;
	// increase the iteration counter
	iteration += 1;
}
} /* namespace networkVMC */
