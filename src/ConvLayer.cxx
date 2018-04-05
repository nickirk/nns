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
  //remember that filters always extend the full depth of the input volume
  //in our case, a filter is always a 2 dimensional matrix, i.e. length and
  //depth(=inputs.size())
  //Damn, it requires higher dimension than 2...
  for(int i(0); i<numFilters; i++){
    Eigen::Map<Eigen::MatrixXd> weightsTmp(adNNP+startPoint,numNrn, 
        inputs.size());
    //weightsTmp /= weightsTmp.size();
    weightsTmp = weightsTmp.unaryExpr(&NormalDistribution);
    weights.push_back(weightsTmp);
    startPoint+=numNrni*inputs.size();
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
  //need to determine the output size. 
  //calculate the output vector size, 
  //(N_in-sizeFilter)/stride+1;

  for (int i(0); i < numFilters; i++){
    for (int j(0); j < sizeAct; j++){		
   //map the input into a vector which is the same size 
   //as the filter and then take dot product. 
   //To map the input into a vector, need the stride size and 
   //the address of the first element of the input.
   //How to get the address of the reference?
   //applying the address-of operator to the reference is the same as taking 
   //the address of the original object. So just add & in front of input
    Eigen::Map<Eigen::VectorXd> inputTmp(&input+startPt, sizeFilter);
    Eigen::VectorXd actTmp(sizeAct);
    actTmp(j)=inputTmp.dot(Filters(i));
    }
  } 
}
