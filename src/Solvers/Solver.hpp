/*
 * Solver.hpp
 *
 *  Created on: Jan 10, 2018
 *      Author: guther
 */

#ifndef SRC_SOLVERS_SOLVER_HPP_
#define SRC_SOLVERS_SOLVER_HPP_

#include "../utilities/TypeDefine.hpp"
#include <Eigen/Dense>

namespace networkVMC{

class State;

// Class doing the minimization in a force-field given by force by different schemes
template<typename T=VecType>
class Solver {
public:
  // all solvers have a train rate/step size parameter
  Solver(double learningRate_):learningRate(learningRate_){};
  virtual ~Solver(){};
  // For first-order solvers, we need the parameters and a force
  // Second order solvers also require the parameterization in the constructor
  // and get the second derivatives on the fly
  // Need samplerType especially in StochasticReconfiguration, because while
  // constructing the S matrix, need to know which getDeriv to call.
  virtual void update(T &w, T const &force, State const &input, SamplerType
    const &samplerType=Markov)=0;

  // set/get the learning rate (the learning rate is an internal parameter,
  // changing it does not change the solver itself)
  double getLearningRate() const {return learningRate;}
  void setLearningRate(double newRate) {learningRate = newRate;}
protected:
  mutable double learningRate;
};

}

#endif /* SRC_SOLVERS_SOLVER_HPP_ */
