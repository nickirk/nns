/*
 * InputLayer.cxx
 * Created on 06.3.2018
 * Author: Ke Liao 
 */

#include <vector>
#include <Eigen/Dense>
#include "InputLayer.hpp"
#include "Determinant.hpp"

InputLayer::InputLayer(std::vector<Eigen::VectorXd> const &inputs_, int size_):
    Layer(inputs_, "Linear"), numNrn(size_){
  activations.resize(1,Eigen::VectorXd::Zero(numNrn));
}

InputLayer::~InputLayer(){

}

void InputLayer::processSignal(detType const det){
  int numStates=det.size();
  for (int state=0; state<numStates; ++state){
    inputs[0](state) = det[state]?1.0:-1.0;
  }
  activations=inputs;
}
