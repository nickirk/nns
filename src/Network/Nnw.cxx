/*
 * Nnw.cxx
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */
#include <vector>
#include <iostream>
#include <random>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <complex>
#include "Nnw.hpp"

namespace networkVMC{
NeuralNetwork::NeuralNetwork():Parametrization(){

  //initial value for NNW para
  numLayers = 0;
  //initial para for Nesterov's accelerated gradient descent
  numNNP = 0;

}

NeuralNetwork::~NeuralNetwork(){
}

//---------------------------------------------------------------------------//

//initialise the network after construction functions are called.
void NeuralNetwork::initialiseNetwork(){
  //calculate the size of the NNP array:
  numNNP=0;
  for (int layer(0); layer<numLayers; ++layer){
    numNNP+=Layers[layer]->getNumPara();
  }
  NNP = Eigen::VectorXd::Ones(numNNP);
  double *adNNP = &NNP(0);
  nablaNNP = Eigen::VectorXd::Zero(numNNP);
  //get the address of nablaNNP
  double *adNablaNNP = &nablaNNP(0);
  int startPoint(0);
  for (int layer(0); layer<numLayers; ++layer){
    Layers[layer]->mapPara(adNNP, adNablaNNP, startPoint);
  }

}
//---------------------------------------------------------------------------//
// construction function of the NNW
void NeuralNetwork::constrInputLayer(int numStates){

  feedIns.resize(1, Eigen::VectorXd::Zero(numStates));
  Layers.addInputLayer(feedIns,numStates);
  numLayers++;
}

//---------------------------------------------------------------------------------------------------//
void NeuralNetwork::constrDenseLayer(
    std::vector<Eigen::VectorXd> const &inputs_, std::string actFunc_,
    int size_
    ){
  Layers.addDenseLayer(inputs_,actFunc_,size_);
  numLayers++;
}

//---------------------------------------------------------------------------------------------------//

void NeuralNetwork::constrConvLayer(
    std::vector<Eigen::VectorXd> const &inputs_,
    std::string actFunc_,
    int numFilters_,
    int lengthFilter_,
    int stride_
    ){
 Layers.addConvLayer(inputs_,actFunc_,numFilters_,lengthFilter_,stride_);
 numLayers++;
}

//---------------------------------------------------------------------------------------------------//

coeffType NeuralNetwork::getCoeff(detType const &det) const{
	//Run the network
  Eigen::VectorXd output=feedForward(det);
	// and extract the coefficient from the last layer
  //std::cout << "Nnw.cxx:getCoeff(): output= " << output << std::endl;
	return coeffType(output(0),output(1));
}

//---------------------------------------------------------------------------//

VecType NeuralNetwork::feedForward(detType const& det) const{
	if(Layers.size()==0) throw EmptyNetworkError();
  // Note that the first layer always needs to have a number
  // of neurons equal to the number of orbitals
  Layers[0]->processSignal(det);
  //std::cout << "Acts layer " << 0 << " =" << std::endl;
  //std::cout<< Layers[0]->getActs()[0] << std::endl;
  //0th layer has no weights!!! The getWeights() for InputLayer is not defined..
  //FixMe, return a 0 vector and throw error.
  //std::cout << "Weights layer " << 0 << " =" << std::endl;
  //std::cout<< Layers[0]->getWeights()[0] << std::endl;
  for (int layer(1); layer < numLayers; ++layer){
      Layers[layer]->processSignal();
  }

  return Layers[numLayers-1]->getActs()[0];
}

//---------------------------------------------------------------------------//

VecType NeuralNetwork::backPropagate(
       VecType const &lastLayerFeedBack
     ){
  //everytime the backPropagate is called, we should reset nabla* to zero.
  nablaNNP *= 0.;
  // Does not work with empty networks
  if(Layers.size()==0) throw EmptyNetworkError();
  Layers[numLayers-1]->backProp(lastLayerFeedBack);
  for (size_t layer(numLayers-2); layer > 0; layer--){
    Layers[layer]->backProp(Layers[layer+1]->getDeltas(),
                            Layers[layer+1]->getWeights());
  }

  return nablaNNP;
}

//---------------------------------------------------------------------------//

Eigen::VectorXd NeuralNetwork::calcNablaPars(
	   State const &inputState,
	   nablaType const &dEdC
     ){
  int numDets = inputState.size();
  Eigen::VectorXd deltaNNP(Eigen::VectorXd::Zero(numNNP));
  for (int epoch=0; epoch < numDets; ++epoch){
    // obtain inputSignals and activations of all layers
    feedForward(inputState.det(epoch));
    // calculate the derivatives of this determinant
    deltaNNP += backPropagate(dEdC[epoch]);
  }
  return deltaNNP;
}

//---------------------------------------------------------------------------------------------------//

VecType NeuralNetwork::calcNablaParsConnected(
	   State const &inputState,
	   nablaType const &dEdC
     ){
  int numDets = inputState.size();
  Eigen::VectorXd deltaNNP(Eigen::VectorXd::Zero(numNNP));
  Eigen::VectorXcd deltaNNPc(numNNP);
  Eigen::Vector2d realMask;
  realMask << 1, 0;
  Eigen::Vector2d imagMask;
  imagMask << 0, 1;
  int pos = 0;

  Eigen::VectorXd deltaNNPTmpPrev(Eigen::VectorXd::Zero(numNNP));
  for (int epoch(0); epoch < numDets; ++epoch){
    Eigen::VectorXd deltaNNPTmp(Eigen::VectorXd::Zero(numNNP));
    if (epoch == 0 || (inputState.det(epoch) != inputState.det(epoch-1))){
      feedForward(inputState.det(epoch));
      deltaNNPc.real() = backPropagate(realMask) ;
      deltaNNPc.imag() = -backPropagate(imagMask) ;
      deltaNNPc *=  dEdC[pos](0);
      deltaNNPc /= std::conj(inputState.coeff(epoch));
      deltaNNPTmp += 2 * deltaNNPc.real();
      pos++;
      std::vector<detType> coupledDets = inputState.coupledDets(epoch);
      std::vector<coeffType > coupledCoeffs = inputState.coupledCoeffs(epoch);
      for (size_t i(0); i < coupledDets.size(); ++i){
        feedForward(coupledDets[i]);
        deltaNNPc.real() = backPropagate(realMask) ;
        deltaNNPc.imag() = backPropagate(imagMask) ;
        deltaNNPc *=  dEdC[pos](0);
        deltaNNPc /= inputState.coeff(epoch);
        deltaNNPTmp += 2 * deltaNNPc.real();
        pos++;
      }
      deltaNNPTmpPrev = deltaNNPTmp;
    }
    else  {deltaNNPTmp = deltaNNPTmpPrev;}
    deltaNNP += deltaNNPTmp;
  }
  deltaNNP /= numDets;
  std::cout << deltaNNP << std::endl;
  return deltaNNP;
}

//---------------------------------------------------------------------------------------------------//
Eigen::MatrixXcd NeuralNetwork::calcdCdwSR(
  State  const &outputState
  ){
//This step produce a complex matrix dCdw. It is done via the same backPropagate
//algorithm, with dEdC set to a vector of (1,0) (real part) and a vector of (0,1)
//(imaginary part). 
  int numDets = outputState.size();
  //Here, not like in backPropagate,
  //the dEdC is just a parameter which we choose to produce the real
  //or imaginary part of dCdw, doesn't mean dEdC itself.
  Eigen::MatrixXcd dCdw(numNNP,numDets);
  Eigen::Vector2d dedc;
  Eigen::Vector2d mask;
  mask << -1., 1.;
  for (int i(0); i <=1; ++i){
  //loop through dedc = (1,0) and (0,1) which will give us the real and imaginary part of
  //dCdw, respectivly.
    Eigen::MatrixXd dCdwTmp(numNNP,numDets);
    dedc << 1, 0;
    dedc += mask * i;
    std::vector<Eigen::VectorXd> dEdC(numDets,dedc);
    //create w_k|--C_i matrix
    for (int epoch(0); epoch < numDets; ++epoch){
      feedForward(outputState.det(i));
      backPropagate(dEdC[epoch]);
      //fill up the dCdwTmp matrix column by column. Column index is C_i, row index is w_k
      dCdwTmp.col(epoch) << nablaNNP; 
    }
    if (i==0) dCdw.real()=dCdwTmp;
    if (i==1) dCdw.imag()=dCdwTmp;
  }
  return dCdw;
}

}
//---------------------------------------------------------------------------------------------------//

