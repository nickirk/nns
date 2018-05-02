/*
 * Solver.hpp
 *
 *  Created on: Jan 10, 2018
 *      Author: guther
 */

#ifndef SRC_SOLVERS_SOLVER_HPP_
#define SRC_SOLVERS_SOLVER_HPP_

#include "../utilities/TypeDefine.hpp"
#include "../utilities/State.hpp"
#include <Eigen/Dense>

namespace networkVMC{

// Class doing the minimization in a force-field given by force by different schemes
class Solver {
public:
  Solver(){};
  virtual ~Solver(){};
  // For first-order solvers, we need the parameters and a force
  // Second order solvers also require the parameterization in the constructor
  // and get the second derivatives on the fly
  virtual void update(VecType &w, VecType const &force, State const &input=State()) const=0;
};

}

#endif /* SRC_SOLVERS_SOLVER_HPP_ */
