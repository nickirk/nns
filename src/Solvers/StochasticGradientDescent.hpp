/*
 * StochasticGradientDescent.hpp
 *
 *  Created on: May 2, 2018
 *      Author: guther
 */

#ifndef SRC_SOLVERS_STOCHASTICGRADIENTDESCENT_HPP_
#define SRC_SOLVERS_STOCHASTICGRADIENTDESCENT_HPP_

#include "Solver.hpp"

namespace networkVMC {
class StochasticGradientDescent: public Solver {
public:
  using T=Eigen::Matrix<paraType, Eigen::Dynamic, 1>;
  StochasticGradientDescent(double gamma_):Solver(gamma_){};
  virtual ~StochasticGradientDescent(){};
  // gradient descent update scheme
  virtual void update(T  &w, T const &force,
		  State const &input, SamplerType const &samplerType=Markov){
  // just walk along the gradient with some stepsize gamma
	  w-=Solver::learningRate*force;
  }
};

} /* namespace networkVMC */

#endif /* SRC_SOLVERS_STOCHASTICGRADIENTDESCENT_HPP_ */
