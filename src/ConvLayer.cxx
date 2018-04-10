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
          double &(actFunc_)(double), 
          int lengthFilter_, int numFilters_, int stride_
          ):
            Layer(inputs_, &(actFunc_)(double)), 
            numFilters(numFilters_),
            depthFilter(inputs_.size()),
            lengthFilter(lengthFilter_),
            stride(stride_){
  //sizeAct is the size of the activations of each filter.
  sizeAct=(inputs[0].size()-lengthFilter)/stride+1;
  z(numFilters,Eigen::VectorXd::Zero(sizeAct));  
  activations(numFilters,Eigen::VectorXd::Zero(sizeAct));
  deltas(numFilters,Eigen::VectorXd::Zero(sizeAct));
  //checked
  numPara = depthFilter*lengthFilter*numFilters+numFilters;
}

void ConvLayer::mapPara(double *adNNP, int &startPoint){
  //remember that filters always extend the full depth of the input volume
  //in our case, a filter consists of several 2 dimensional matrix, i.e. 
  //std::vecotr<Eigen::MatrixXd> (we restrict ourselves for now to one D 
  //system, so the matrix have dimension lx1). 
  for(int i(0); i<numFilters; i++){
    std::vector<Eigen::Map<Eigen::MatrixXd>> filterTmp;
    for (int j(0); j < depthFilter; j++){
      Eigen::Map<Eigen::MatrixXd> weightsTmp(adNNP+startPoint,lengthFilter, 
          1);
      //weightsTmp /= weightsTmp.size();
      weightsTmp = weightsTmp.unaryExpr(&NormalDistribution);
      filterTmp.push_back(weightsTmp);
      startPoint+=lengthFilter;
    }
    weights.push_back(filterTmp);
    Eigen::Map<Eigen::VectorXd> biaseTmp(adNNP+startPoint,1); 
    //biaseTmp /= biaseTmp.size();
    biaseTmp = biaseTmp.unaryExpr(&NormalDistribution);
    biases.push_back(biaseTmp); 
    startPoint+=1;
  }
}

std::vector<Eigen::VectorXd> ConvLayer::processSignal(){
  //need to determine the output size. 
  //calculate the output vector size, 
  //(N_in-lengthFilter)/stride+1;
  for (int i(0); i < numFilters; i++){
    int startPt(0);
    for (int j(0); j < sizeAct; j++){
      for (int k(0); k < depthFilter; k++){    
        //map the input into a vector which is the same size 
        //as the filter and then take dot product. 
        //To map the input into a vector, need the stride size and 
        //the address of the first element of the input.
        //How to get the address of a reference?
        //applying the address-of operator to the reference is 
        //the same as taking 
        //the address of the original object. So just add & in front of input
        Eigen::Map<Eigen::VectorXd> inputTmp(&inputs[k]+startPt,lengthFilter);
        // dot product with a reversed weight matrix, due to the mathematical 
        // definition of convolution. Since here the weight matrix is nx1, 
        // which is column major, so reverse colwise.
        z[i](j)+=inputTmp.transpose() * weights[i][k].colwise().reverse()
                  +biases[i];
      }
      startPt+=stride;
    }
    activations[i] = z[i];
    activations[i].unaryExpr(&actFunc);
  } 

}

