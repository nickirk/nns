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
template <typename F=std::complex<double>, typename coeffType=std::complex<double>>
class StochasticReconfiguration: public Solver<F, coeffType> {
public:
  using T=Eigen::Matrix<F, Eigen::Dynamic, 1>;
  StochasticReconfiguration(Parametrization<F, coeffType> &NNW_, double gamma_):
	  Solver<F, coeffType>(gamma_),NNW(NNW_),iteration(0),internalSolver(gamma_) {};
  virtual ~StochasticReconfiguration();

  virtual void update(T &w, T const &force, State<coeffType> const &input, SamplerType
      const &samplerType=Markov);
private:
  // need the parametrisation to get derivatives in order to construct  the 
  // S matrix
  Parametrization<F, coeffType> &NNW;
  // and an iteration counter
  int iteration;
  ADAM<F, coeffType> internalSolver;
};

} /* namespace networkVMC */

#endif /* SRC_SOLVERS_STOCHASTICRECONFIGURATION_HPP_ */
