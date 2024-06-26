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
  void update(paraVector &w, paraVector const &force, State const &input, SamplerType
      const &samplerType=Markov);
private:
  // number of parameters to be optimized
  int numPars;
  // iteration counter, as the internal parameters are adjusted over the course
  // of optimization
  int iteration;

  // some internal ADAM parameters
  paraVector m,m1,v,v1;
  double beta1, beta2;
};

} /* namespace networkVMC */

#endif /* SRC_SOLVERS_ADAM_HPP_ */
