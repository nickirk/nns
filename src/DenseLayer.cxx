/*
 * DenseLayer.hpp
 * Created on 06.3.2018
 * Author: Ke Liao 
 */

#include <vector>
#include <Eigen/Dense>
#include "DenseLayer.hpp"

DenseLayer(std::vector<Eigen::MatrixXd> const &inputs_, 
           double &(actFunc_)(double), int size_)
  Layer(inputs_, actFunc_), numNrn(size_){
  for (size_t i(0); i<inputs.size(); i++){
    z.push_back(Eigen::VectorXd::Zero(numNrn));  
  }
  activations.push_back(Eigen::VectorXd::Zero(numNrn));
}

void DenseLayer::mapPara(double *adNNP, int &startPoint){
  for(size_t i(0); i<inputs.size(); i++){
    weights.push_back(Eigen::Map<Eigen::MatrixXd>(adNNP+startPoint, numNrn, 
      inputs[i].size());
    startPoint+=numNrn*inputs[i].size();
  }
  biases.push_back(Eigen::Map(Eigen::VectorXd)(adNNP+startPoint, numNrn)); 
  startPoint+=numNrn*inputs[i].size();
}

void DenseLayer::processSignal(){
  activations[0]=Eigen::VectorXd::Zero(numNrn);
  for(size_t i(0); i<inputs.size(); i++){
    z[i] = weights[i]*inputs[i]+biases[0]; 
    activations[0]+=z[i];
  }
  activations[0]=activations[0].unaryExpr(actFunc);
}
