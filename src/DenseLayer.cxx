/*
 * DenseLayer.cxx
 * Created on 06.3.2018
 * Author: Ke Liao 
 */

#include <vector>
#include <Eigen/Dense>
#include "DenseLayer.hpp"

DenseLayer(std::vector<Eigen::VectorXd> const &inputs_, 
           double &(actFunc_)(double), int size_):
  Layer(inputs_, &(actFunc_)(double)) numNrn(size_){
  z(inputs.size(),Eigen::VectorXd::Zero(numNrn));  
  activations(1,Eigen::VectorXd::Zero(numNrn));
  deltas(1,Eigen::VectorXd::Zero(numNrn));
  numPara = numNrn*inputs[0].size()*inputs.size()+numNrn;
}

void DenseLayer::mapPara(double *adNNP, double *adNablaNNP, int &startPoint){
  std::vector<Eigen::Map<Eigen::MatrixXd>> weightsTmp;
  std::vector<Eigen::Map<Eigen::MatrixXd>> NablaWeightsTmp;
  for(size_t i(0); i<inputs.size(); i++){
    Eigen::Map<Eigen::MatrixXd> weightTmp(adNNP+startPoint,numNrn,
		  input[i].size());
    Eigen::Map<Eigen::MatrixXd> NablaWeightTmp(adNablaNNP+startPoint,numNrn,
		  input[i].size());
    //weightsTmp /= weightsTmp.size();
    weightTmp = weightTmp.unaryExpr(&NormalDistribution);
    weightsTmp.push_back(weightsTmp);
    NablaWeightsTmp.push_back(NablaWeightsTmp);
    startPoint += numNrn*inputs[i].size();
  }
  weights.push_back(weightsTmp);
  NablaWeights.push_back(NablaWeightsTmp);
  Eigen::Map<Eigen::VectorXd> biaseTmp(adNNP+startPoint,numNrn); 
  Eigen::Map<Eigen::VectorXd> NablaBiaseTmp(adNablaNNP+startPoint,numNrn); 
  //biaseTmp /= biaseTmp.size();
  biaseTmp = biaseTmp.unaryExpr(&NormalDistribution);
  biases.push_back(biaseTmp); 
  NablaBiases.push_back(NablaBiaseTmp); 
  startPoint+=numNrn;
}

void DenseLayer::processSignal(){
  activations[0]=Eigen::VectorXd::Zero(numNrn);
  for(size_t i(0); i<inputs.size(); i++){
    z[i] = weights[0][i]*inputs[i]+biases[0]; 
    activations[0]+=z[i];
  }
  activations[0]=activations[0].unaryExpr(actFunc);
}


void DenseLayer::backProp(
    std::vector<Eigen::VectorXd> prevDelta, 
    weightType &prevWeights;
    ){
    deltas[0] = prevWeights[0][0].transpose() * deltaPreviousLayer; 
    deltas[0] = deltas[0].array()* z[0].unaryExpr(&actFuncPrime).array();
    nablaBiases[0] = delta[0];
    //get a weight matrix. \partial C/\partial w^{l}_{jk} = a^{l-1}_k \delta_j^l 
    //the layer here refers to the lth layer of Biases and weights, so for
    //activation layer refers to the l-1th layer.
    for (size_t i(0); i < inputs.size(); i++)
      nablaWeights[0][i] = delta[0] * inputs[i].transpose();
}
