/*
 * AcceleratedGradientDescent.hpp
 *
 *  Created on: May 2, 2018
 *      Author: guther
 */

#ifndef SRC_SOLVERS_ACCELERATEDGRADIENTDESCENT_HPP_
#define SRC_SOLVERS_ACCELERATEDGRADIENTDESCENT_HPP_

#include "Solver.hpp"

namespace networkVMC {

// Nesterov's Accelerated Gradient Descent solver
class AcceleratedGradientDescent: public Solver {
public:
  AcceleratedGradientDescent(double learningRate);
  virtual ~AcceleratedGradientDescent();
  // implementation of the update method
  void update(VecType &w, VecType const &force) const;
private:
  double learningRate;
  // number of parameters to optimize (will be automatically determined)
  int numPars;
  // some obscure parameters of the accelerated gradient descent
  mutable double lambdaS1, lambdaS, gammaS, gammaS1;
  mutable VecType yS, yS1, Egz2;
  // auxiliary variable to check if we already know the number of parameters
  mutable bool uninitialized;
};

} /* namespace networkVMC */

#endif /* SRC_SOLVERS_ACCELERATEDGRADIENTDESCENT_HPP_ */
