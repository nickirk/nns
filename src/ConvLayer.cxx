/*
 * ConvLayer.cxx
 * Created on 27.3.2018
 * Author: Ke Liao 
 */
#include <vector>
#include <Eigen/Dense>
#include "ConvLayer.hpp"


ConvLayer(
          std::vector<Eigen::VectorXd> const &inputs_, 
          double &(actFunc_)(double), int size_,
          int numFilters_, int stride_
          ):
            Layer(inputs_, &(actFunc_)(double)), 
            numNrn(size_),
            numFilters(numFilters_),
            stride(stride_){
  z(inputs.size(),Eigen::VectorXd::Zero(numNrn));  
  activations(1,Eigen::VectorXd::Zero(numNrn));
  numPara = numNrn*numFilters+numFilters;

}

void ConvLayer::mapPara(double *adNNP, int &startPoint){
  for(int i(0); i<numFilters; i++){
    Eigen::Map<Eigen::MatrixXd> weightsTmp(adNNP+startPoint,numNrn, 1);
    //weightsTmp /= weightsTmp.size();
    weightsTmp = weightsTmp.unaryExpr(&NormalDistribution);
    weights.push_back(weightsTmp);
    startPoint+=numNrn;
  }
  Eigen::Map<Eigen::VectorXd> biaseTmp(adNNP+startPoint,1); 
  //biaseTmp /= biaseTmp.size();
  biaseTmp = biaseTmp.unaryExpr(&NormalDistribution);
  biases.push_back(Eigen::Map(Eigen::VectorXd)(adNNP+startPoint, 1)); 
  startPoint+=1;
}

std::vector<Eigen::VectorXd> ConvLayer::convolve(
    Eigen::VectorXd const &input, int startPt
    ){
  for (int i(0); i < numFilters; i++){
    
  } 
  
}
