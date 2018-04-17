/*
 * ConvLayer.cxx
 * Created on 27.3.2018
 * Author: Ke Liao 
 */
#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include "ConvLayer.hpp"


ConvLayer::ConvLayer(
          std::vector<Eigen::VectorXd> const &inputs_, 
          std::string actFunc_, 
          int numFilters_,
          int lengthFilter_,
          int stride_
          ):
            Layer(inputs_, actFunc_), 
            numFilters(numFilters_),
            depthFilter(inputs_.size()),
            lengthFilter(lengthFilter_),
            stride(stride_){
  //sizeAct is the size of the activations of each filter.
  sizeAct=(inputs[0].size()-lengthFilter)/stride+1;
  z.resize(numFilters,Eigen::VectorXd::Zero(sizeAct));  
  activations.resize(numFilters,Eigen::VectorXd::Zero(sizeAct));
  deltas.resize(numFilters,Eigen::VectorXd::Zero(sizeAct));
  //checked
  numPara = depthFilter*lengthFilter*numFilters+numFilters;
}

ConvLayer::~ConvLayer(){

}

void ConvLayer::mapPara(double *adNNP, double *adNablaNNP, int &startPoint){
  //remember that filters always extend the full depth of the input volume
  //in our case, a filter consists of several 2 dimensional matrix, i.e. 
  //std::vecotr<Eigen::MatrixXd> (we restrict ourselves for now to one D 
  //system, so the matrix have dimension lx1). 
  for(int i(0); i<numFilters; i++){
    std::vector<Eigen::Map<Eigen::MatrixXd>> filterTmp;
    std::vector<Eigen::Map<Eigen::MatrixXd>> nablaFilterTmp;
    for (int j(0); j < depthFilter; j++){
      Eigen::Map<Eigen::MatrixXd> weightsTmp(adNNP+startPoint,lengthFilter, 1);
      Eigen::Map<Eigen::MatrixXd> nablaWeightsTmp(adNablaNNP+startPoint,lengthFilter,
          1);
      //weightsTmp /= weightsTmp.size();
      weightsTmp = weightsTmp.unaryExpr(&NormalDistribution);
      //nablaWeightsTmp = nablaWeightsTmp.unaryExpr(&NormalDistribution);
      filterTmp.push_back(weightsTmp);
      nablaFilterTmp.push_back(nablaWeightsTmp);
      startPoint+=lengthFilter;
    }
    weights.push_back(filterTmp);
    nablaWeights.push_back(nablaFilterTmp);
    Eigen::Map<Eigen::VectorXd> biaseTmp(adNNP+startPoint,1); 
    Eigen::Map<Eigen::VectorXd> nablaBiaseTmp(adNablaNNP+startPoint,1);
    //biaseTmp /= biaseTmp.size();
    biaseTmp = biaseTmp.unaryExpr(&NormalDistribution);
    biases.push_back(biaseTmp); 
    nablaBiases.push_back(nablaBiaseTmp);
    startPoint+=1;
  }
}

void ConvLayer::processSignal(){
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
        Eigen::Map<const Eigen::VectorXd> inputTmp(&inputs[k](0)+startPt,lengthFilter);
        // dot product with a reversed weight matrix, due to the mathematical 
        // definition of convolution. Since here the weight matrix is nx1, 
        // which is column major, so reverse colwise.
        z[i](j)+=(weights[i][k].colwise().reverse().transpose()*inputTmp
                  +biases[i])(0);
      }
      startPt+=stride;
    }
    activations[i] = z[i];
    activations[i]=activations[i].unaryExpr(actFunc);
  } 
}

void ConvLayer::backProp(
    std::vector<Eigen::VectorXd> const &prevDelta,
    weightType const &prevWeights
    ){
  for (int i(0); i < numFilters; i++){
    deltas[i] = prevWeights[i][0].transpose() * prevDelta[0];
    deltas[i] = deltas[i].array()* (z[i].unaryExpr(actFuncPrime)).array();
    //bias of the convLayer is just delta. 
    nablaBiases[i](0) = deltas[i].sum();
    //the complicated part is the weights. We need to retrieve the 
    //nablaWeigths for each depth. Since we did a convolution, the 
    //nablaWeights must contain a sum over a certain region of the 
    //inputs from last layer. 
    for (int j(0); j < depthFilter; j++){
      for (int row(0); row < nablaWeights[i][j].rows(); row++){
        for (int col(0); col < nablaWeights[i][j].cols(); col++){
        	std::cout << "ConvLayer.cxx : row=" << row << std::endl;
          Eigen::Map<const Eigen::VectorXd, 0, Eigen::InnerStride<Eigen::Dynamic>>
		  	  inputTmp(&inputs[j](0)+row,deltas[i].size(), Eigen::InnerStride<Eigen::Dynamic>(stride));
          std::cout << "ConvLayer.cxx: deltas[i].size()=" << deltas[i].size() << std::endl;
          std::cout << "ConvLayer.cxx: deltas[i]=\n" << deltas[i] << std::endl;
          std::cout << "ConvLayer.cxx: inputTmp=\n" << inputTmp.transpose() << std::endl;
          std::cout << "ConvLayer.cxx: nablaWeights.size()=" << nablaWeights[i][j].size() << std::endl;
          nablaWeights[i][j](row,col) = (inputTmp.transpose()*(deltas[i].colwise().reverse()))(0);
        }
      }
    }
  }
}
