/*
 * StochasticReconfiguration.hpp
 *
 *  Created on: May 2, 2018
 *      Author: guther
 */

#ifndef SRC_SOLVERS_STOCHASTICRECONFIGURATION_HPP_
#define SRC_SOLVERS_STOCHASTICRECONFIGURATION_HPP_

#include "Solver.hpp"
#include "ADAM.hpp"
#include "../Network/Parametrization.hpp"
#include <Eigen/IterativeLinearSolvers>

namespace networkVMC {
template <typename T=VecType>
class StochasticReconfiguration: public Solver<T> {
public:
  StochasticReconfiguration(Parametrization<T> &NNW_, double gamma_):
	  Solver<T>(gamma_),NNW(NNW_),iteration(0),internalSolver(gamma_) {};
  virtual ~StochasticReconfiguration();

  virtual void update(T &w, T const &force, State const &input, SamplerType 
      const &samplerType=Markov);
private:
  // need the parametrisation to get derivatives in order to construct  the 
  // S matrix
  Parametrization<T> &NNW;
  // and an iteration counter
  int iteration;
  ADAM<T> internalSolver;
};

} /* namespace networkVMC */

#endif /* SRC_SOLVERS_STOCHASTICRECONFIGURATION_HPP_ */
