/*
 * ADAM.hpp
 *
 *  Created on: May 2, 2018
 *      Author: guther
 */

#ifndef SRC_SOLVERS_ADAM_HPP_
#define SRC_SOLVERS_ADAM_HPP_

#include "Solver.hpp"

namespace networkVMC {

class ADAM: public Solver {
public:
  ADAM(double learningRate);
  virtual ~ADAM();
  // implementation of the update method
  void update(VecType &w, VecType const &force) const;
private:
  // the learningRate (something like a step size)
  double learningRate;
  // number of parameters to be optimized
  int numPars;
  // iteration counter, as the internal parameters are adjusted over the course
  // of optimization
  mutable int iteration;

  // some internal ADAM parameters
  mutable VecType m,m1,v,v1;
  mutable double beta1, beta2;
};

} /* namespace networkVMC */

#endif /* SRC_SOLVERS_ADAM_HPP_ */
