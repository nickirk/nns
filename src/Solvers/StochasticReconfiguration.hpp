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
class StochasticReconfiguration: public Solver<> {
public:
  StochasticReconfiguration(Parametrization<T> &par_, double gamma_):
	  Solver(gamma_),par(par_),iteration(0) {};
  virtual ~StochasticReconfiguration();

  // thus, we only have this version
  virtual void update(T &w, T const &force, State const &input);
private:
  // this is a second order solver, so we need the parametrization
  Parametrization<T> &par;
  // and an iteration counter
  int iteration;
};

} /* namespace networkVMC */

#endif /* SRC_SOLVERS_STOCHASTICRECONFIGURATION_HPP_ */
