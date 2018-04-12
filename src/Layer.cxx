/*
 * Layer.cxx
 * Created on 27.3.2018
 * Author: Ke Liao 
 */


#include <vector>
#include <Eigen/Dense>
#include <string>
#include "math/MathFunctions.hpp"
#include "Layer.hpp"

Layer::Layer(std::vector<Eigen::VectorXd> const &inputs_, int actFunc_):
inputs(inputs_){
  //default activation function
  actFunc=&Tanh;
  actFuncPrime=&TanhPrime;
  numPara=0;
  if (actFunc_ == Linear){
    actFunc(&Linear(double));
    actFuncPrime(&LinearPrime(double));
  }
  else if (actFunc_ == Tanh){
    actFunc(&Tanh(double));
    actFuncPrime(&TanhPrime(double));
  }
  else if (actFunc_ == Sigmoid){
    actFunc(&Sigmoid(double));
    actFuncPrime(&SigmoidPrime(double));
  }
  else if (actFunc_ == Sigmoid){
    actFunc(&Sigmoid(double));
    actFuncPrime(&SigmoidPrime(double));
  }
  else if (actFunc_ == Rectifier){
    actFunc(&Rectifier(double));
    actFuncPrime(&RectifierPrime(double));
  }
}
