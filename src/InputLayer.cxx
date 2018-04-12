/*
 * InputLayer.cxx
 * Created on 06.3.2018
 * Author: Ke Liao 
 */

#include <vector>
#include <Eigen/Dense>
#include "InputLayer.hpp"
#include "Determinant.hpp"

InputLayer(std::vector<Eigen::VectorXd> const &inputs_, int size_): 
    Layer(inputs_, &Linear(double)), numNrn(size_){
  activations(1,Eigen::VectorXd::Zero(numNrn));
}

void InputLayer::processSignal(){
  activations=inputs;
}
