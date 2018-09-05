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

template <typename F, typename coeffType>
Layer<F, coeffType>::Layer(std::vector<Eigen::VectorXd> const &inputs_, std::string actFunc_):
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
  else if (actFunc_=="LogCosh"){
    actFunc=&LogCosh;
    actFuncPrime=&LogCoshPrime;
  }
  else throw errors::ActFuncDoNotExist(actFunc_);
}

template <typename F, typename coeffType>
Layer<F, coeffType>::~Layer(){
}

// Processing an input signal can be done by any Layer, given that the
// number of neurons is equal to the number of orbitals
template <typename F, typename coeffType>
void Layer<F, coeffType>::processSignal(detType const &det) const{
  int numStates=det.size();
  // This only works if the fed determinant is valid
  if(activations[0].size()!=numStates) throw errors::InvalidDeterminantError(det);
  //set the activations into the determinants
  for (int state=0; state<numStates; ++state){
    activations[0](state) = det[state]?1.0:-1.0;
  }
}

}
