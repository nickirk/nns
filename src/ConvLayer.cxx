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
          int numFilters_, int lengthFilter_, int stride_
          ):
            Layer(inputs_, &(actFunc_)(double)), 
            numNrn(size_),
            numFilters(numFilters_),
            depthFilter(inputs_.size()),
            lengthFilter(lengthFilter_),
            stride(stride_){
  z(inputs.size(),Eigen::VectorXd::Zero(numNrn));  
  activations(1,Eigen::VectorXd::Zero(numNrn));
  numPara = numNrn*numFilters+numFilters;

}

void ConvLayer::mapPara(double *adNNP, int &startPoint){
  //remember that filters always extend the full depth of the input volume
  //in our case, a filter consists of several 2 dimensional matrix, i.e. 
  //std::vecotr<Eigen::MatrixXd> (we restrict ourselves for now to one D 
  //system, so the matrix have dimension lx1). 
  //The size of the vector is the same as the 
  //inputs.size().
  for(int i(0); i<numFilters; i++){
    std::vector<Eigen::Map<Eigen::MatrixXd>> filterTmp;
    for (int j(0); j < depthFilter; j++){
      Eigen::Map<Eigen::MatrixXd> weightsTmp(adNNP+startPoint,numNrn, 
          1);
      //weightsTmp /= weightsTmp.size();
      weightsTmp = weightsTmp.unaryExpr(&NormalDistribution);
      filterTmp.push_back(weightsTmp);
      startPoint+=numNrn;
    }
    weights.push_back(filterTmp);
    Eigen::Map<Eigen::VectorXd> biaseTmp(adNNP+startPoint,1); 
    //biaseTmp /= biaseTmp.size();
    biaseTmp = biaseTmp.unaryExpr(&NormalDistribution);
    biases.push_back(biaseTmp); 
    startPoint+=1;
  }
}

std::vector<Eigen::VectorXd> ConvLayer::convolve(){
  //need to determine the output size. 
  //calculate the output vector size, 
  //(N_in-sizeFilter)/stride+1;

  for (int i(0); i < numFilters; i++){
    int startPt(0);
    for (int j(0); j < sizeAct; j++){
      for (int k(0); k< depthFilter; k++){    
        //map the input into a vector which is the same size 
        //as the filter and then take dot product. 
        //To map the input into a vector, need the stride size and 
        //the address of the first element of the input.
        //How to get the address of the reference?
        //applying the address-of operator to the reference is 
        //the same as taking 
        //the address of the original object. So just add & in front of input
         Eigen::Map<Eigen::VectorXd> inputTmp(&inputs[k]+startPt,lengthFilter);
         Eigen::VectorXd actTmp(sizeAct);
         actTmp(j)+=inputTmp.dot(Filters(i));
      }
      startPt+=stride;
    }
  } 
}
