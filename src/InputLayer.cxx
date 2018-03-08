/*
 * InputLayer.hpp
 * Created on 06.3.2018
 * Author: Ke Liao 
 */

#include <vector>
#include <Eigen/Dense>
#include "InputLayer.hpp"
#include "Determinant.hpp"

InputLayer(std::vector<Eigen::VectorXd> inputs_, int size_): 
    Layer(inputs_), numNrn(size_){
  activations.push_back(Eigen::VectorXd::Zero(numNrn));
}

void InputLayer::processSignal(detType const &det){
  activations[0]=Eigen::VectorXd::Zero(numNrn);
  int numStates=det.size();
  for (int state=0; state < numStates; ++state){
    activations[0](state) = det[state]?1.0:-1.0;
    //here the inputs might not be needed.
    inputs[0](state) = det[state]?1.0:-1.0;
  }
}
