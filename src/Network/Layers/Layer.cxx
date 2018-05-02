/*
 * Layer.cxx
 * Created on 27.3.2018
 * Author: Ke Liao 
 */


#include <vector>
#include <Eigen/Dense>
#include <string>
#include "Layer.hpp"

namespace networkVMC{

Layer::Layer(std::vector<Eigen::VectorXd> const &inputs_, std::string actFunc_):
inputs(inputs_){
  //default activation function
  numPara=0;
  if (actFunc_ == "Linear"){
    actFunc=&Linear;
    actFuncPrime=&LinearPrime;
  }
  else if (actFunc_ == "Tanh"){
    actFunc=&Tanh;
    actFuncPrime=&TanhPrime;
  }
  else if (actFunc_ == "Sigmoid"){
    actFunc=&Sigmoid;
    actFuncPrime=&SigmoidPrime;
  }
  else if (actFunc_ == "Rectifier"){
    actFunc=&Rectifier;
    actFuncPrime=&RectifierPrime;
  }
  else throw ActFuncDoNotExist(actFunc_);
}

Layer::~Layer(){

}

// Processing an input signal can be done by any Layer, given that the
// number of neurons is equal to the number of orbitals
void Layer::processSignal(detType const &det) const{
  int numStates=det.size();
  //set the activations into the determinants
  for (int state=0; state<numStates; ++state){
    activations[0](state) = det[state]?1.0:-1.0;
  }
}

}
