/*
 * DenseLayer.hpp
 * Created on 06.3.2018
 * Author: Ke Liao 
 */

#include <vector>
#include <Eigen/Dense>
#include "DenseLayer.hpp"

DenseLayer(std::vector<Eigen::MatrixXd> const &inputs_, 
           double &(actFunc_)(double), int size_):
  Layer(inputs_) actFunc(actFunc_), numNrn(size_){
  z(inputs.size(),Eigen::VectorXd::Zero(numNrn));  
  activations(1,Eigen::VectorXd::Zero(numNrn));
}

void DenseLayer::mapPara(double *adNNP, int &startPoint){
  for(size_t i(0); i<inputs.size(); i++){
    Eigen::Map<Eigen::MatrixXd> weightsTmp(adNNP+startPoint,numNrn,
		  input[i].size());
    //weightsTmp /= weightsTmp.size();
    weightsTmp = weightsTmp.unaryExpr(&NormalDistribution);
    weights.push_back(weightsTmp);
    startPoint+=numNrn*inputs[i].size();
  }
  Eigen::Map<Eigen::VectorXd> biaseTmp(adNNP+startPoint,numNrn); 
  //biaseTmp /= biaseTmp.size();
  biaseTmp = biaseTmp.unaryExpr(&NormalDistribution);
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
