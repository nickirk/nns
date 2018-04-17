/*
 * DenseLayer.cxx
 * Created on 06.3.2018
 * Author: Ke Liao 
 */

#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include "DenseLayer.hpp"

DenseLayer::DenseLayer(std::vector<Eigen::VectorXd> const &inputs_, 
    std::string actFunc_, int size_):
  Layer(inputs_, actFunc_), numNrn(size_){
  z.resize(inputs.size(),Eigen::VectorXd::Zero(numNrn));
  activations.resize(1,Eigen::VectorXd::Zero(numNrn));
  deltas.resize(1,Eigen::VectorXd::Zero(numNrn));
  numPara = numNrn*inputs[0].size()*inputs.size()+numNrn;
}

DenseLayer::~DenseLayer(){
}

void DenseLayer::mapPara(double *adNNP, double *adNablaNNP, int &startPoint){
  std::vector<Eigen::Map<Eigen::MatrixXd>> weightsTmp;
  std::vector<Eigen::Map<Eigen::MatrixXd>> NablaWeightsTmp;
  for(size_t i(0); i<inputs.size(); i++){
    Eigen::Map<Eigen::MatrixXd> weightTmp(adNNP+startPoint,numNrn,
		  inputs[i].size());
    Eigen::Map<Eigen::MatrixXd> NablaWeightTmp(adNablaNNP+startPoint,numNrn,
		  inputs[i].size());
    //weightsTmp /= weightsTmp.size();
    weightTmp = weightTmp.unaryExpr(&NormalDistribution);
    weightsTmp.push_back(weightTmp);
    NablaWeightsTmp.push_back(NablaWeightTmp);
    startPoint += numNrn*inputs[i].size();
  }
  weights.push_back(weightsTmp);
  nablaWeights.push_back(NablaWeightsTmp);
  Eigen::Map<Eigen::VectorXd> biaseTmp(adNNP+startPoint,numNrn); 
  Eigen::Map<Eigen::VectorXd> NablaBiaseTmp(adNablaNNP+startPoint,numNrn); 
  //biaseTmp /= biaseTmp[0].size();
  biaseTmp = biaseTmp.unaryExpr(&NormalDistribution);
  biases.push_back(biaseTmp); 
  nablaBiases.push_back(NablaBiaseTmp);
  startPoint+=numNrn;
}

void DenseLayer::processSignal(){
  activations[0]=Eigen::VectorXd::Zero(numNrn);
  for(size_t i(0); i<inputs.size(); i++){
    z[i] = weights[0][i]*inputs[i];
    activations[0]+=z[i];
  }
  activations[0]+=biases[0];
  activations[0]=activations[0].unaryExpr(actFunc);
}


void DenseLayer::backProp(std::vector<Eigen::VectorXd> const &prevDelta,
    weightType const &prevWeights){
    deltas[0] = prevWeights[0][0].transpose() * prevDelta[0];
    deltas[0] = deltas[0].array()* (z[0].unaryExpr(actFuncPrime)).array();
    nablaBiases[0] = deltas[0];
    //get a weight matrix. \partial C/\partial w^{l}_{jk} = a^{l-1}_k \delta_j^l 
    //the layer here refers to the lth layer of Biases and weights, so for
    //activation layer refers to the l-1th layer.
    for (size_t i(0); i < inputs.size(); i++)
      nablaWeights[0][i] = deltas[0] * inputs[i].transpose();
}

void DenseLayer::backProp(
    Eigen::VectorXd const &prevDelta
    ){
    deltas[0] = 
    (prevDelta.array() * 
    (z[0].unaryExpr(actFuncPrime)).array()).matrix();
    nablaBiases[0] = deltas[0];
    //get a weight matrix. \partial C/\partial w^{l}_{jk} = a^{l-1}_k \delta_j^l 
    //the layer here refers to the lth layer of Biases and weights, so for
    //activation layer refers to the l-1th layer.
    for (size_t i(0); i < inputs.size(); i++)
      nablaWeights[0][i] = deltas[0] * inputs[i].transpose();
}
