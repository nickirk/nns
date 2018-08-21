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
template <typename F=std::complex<double>, typename coeffType=std::complex<double>>
class StochasticGradientDescent: public Solver<F, coeffType> {
public:
  using T=Eigen::Matrix<F, Eigen::Dynamic, 1>;
  StochasticGradientDescent(double gamma_):Solver<F, coeffType>(gamma_){};
  virtual ~StochasticGradientDescent(){};
  // gradient descent update scheme
  virtual void update(T  &w, T const &force,
		  State<coeffType> const &input, SamplerType const &samplerType=Markov){
  // just walk along the gradient with some stepsize gamma
	  w-=Solver<F>::learningRate*force;
  }
};

} /* namespace networkVMC */

#endif /* SRC_SOLVERS_STOCHASTICGRADIENTDESCENT_HPP_ */
