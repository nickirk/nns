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
#include "../utilities/MatrixTypeDefine.hpp"

namespace networkVMC {

class Parametrization;

class StochasticReconfiguration: public Solver {
public:
  StochasticReconfiguration(Parametrization &NNW_, double gamma_):
	  Solver(gamma_),NNW(NNW_),iteration(0),internalSolver(gamma_) {};
  virtual ~StochasticReconfiguration();

  virtual void update(paraVector &w, paraVector const &force, State const &input, SamplerType
      const &samplerType=Markov);
private:
  // need the parametrisation to get derivatives in order to construct  the 
  // S matrix
  Parametrization &NNW;
  // and an iteration counter
  int iteration;
  ADAM internalSolver;
};

} /* namespace networkVMC */

#endif /* SRC_SOLVERS_STOCHASTICRECONFIGURATION_HPP_ */
