/*
 * Layer.cxx
 * Created on 27.3.2018
 * Author: Ke Liao 
 */


#include "Layer.hpp"

#include <vector>
#include <Eigen/Dense>
#include <string>

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
