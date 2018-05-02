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
  // all solvers have a train rate/step size parameter
  Solver(double learningRate_):learningRate(learningRate_){};
  virtual ~Solver(){};
  // For first-order solvers, we need the parameters and a force
  // Second order solvers also require the parameterization in the constructor
  // and get the second derivatives on the fly
  virtual void update(VecType &w, VecType const &force, State const &input=State())=0;

  // set/get the learning rate (the learning rate is an internal parameter,
  // changing it does not change the solver itself)
  double getLearningRate() const {return learningRate;}
  void setLearningRate(double newRate) {learningRate = newRate;}
protected:
  mutable double learningRate;
};

}

#endif /* SRC_SOLVERS_SOLVER_HPP_ */
