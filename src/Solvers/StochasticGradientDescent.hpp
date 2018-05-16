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
template <typename T=VecType>
class StochasticGradientDescent: public Solver<T> {
public:
  StochasticGradientDescent(double gamma_):Solver<T>(gamma_){};
  virtual ~StochasticGradientDescent(){};
  // gradient descent update scheme
  virtual void update(T  &w, T const &force,
		  State const &input=State()){
  // just walk along the gradient with some stepsize gamma
	  w-=Solver<T>::learningRate*force;
  }
};

} /* namespace networkVMC */

#endif /* SRC_SOLVERS_STOCHASTICGRADIENTDESCENT_HPP_ */
