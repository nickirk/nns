/*
 * StochasticReconfiguration.hpp
 *
 *  Created on: May 2, 2018
 *      Author: guther
 */

#ifndef SRC_SOLVERS_STOCHASTICRECONFIGURATION_HPP_
#define SRC_SOLVERS_STOCHASTICRECONFIGURATION_HPP_

#include "Solver.hpp"
#include "../Network/Parametrization.hpp"

namespace networkVMC {
template <typename T=VecType>
class StochasticReconfiguration: public Solver<T> {
public:
  StochasticReconfiguration(Parametrization<T> &NNW_, double gamma_):
	  Solver<T>(gamma_),NNW(NNW_),iteration(0) {};
  virtual ~StochasticReconfiguration();

  // thus, we only have this version
  virtual void update(T &w, T const &force, State const &input);
private:
  // this is a second order solver, so we need the parametrization
  Parametrization<T> &NNW;
  // and an iteration counter
  int iteration;
};

} /* namespace networkVMC */

#endif /* SRC_SOLVERS_STOCHASTICRECONFIGURATION_HPP_ */
