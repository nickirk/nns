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

class StochasticReconfiguration: public Solver {
public:
  StochasticReconfiguration(Parametrization const &par_, double gamma_):
	  par(par_),gamma(gamma_),iteration(0) {};
  virtual ~StochasticReconfiguration();

  // thus, we only have this version
  virtual void update(VecType &w, VecType const &force, State const &input)
    const;
private:
  // this is a second order solver, so we need the parametrization
  Parametrization &par;
  // the stepsize
  double gamma;
  // and an iteration counter
  int iteration;
};

} /* namespace networkVMC */

#endif /* SRC_SOLVERS_STOCHASTICRECONFIGURATION_HPP_ */
