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

template<typename T>
NeuralNetwork<T>::NeuralNetwork(){
  //initial value for NNW para
  numLayers = 0;
  //initial para for Nesterov's accelerated gradient descent
  numNNP = 0;

}

template <typename T>
NeuralNetwork<T>::~NeuralNetwork(){
}

//---------------------------------------------------------------------------//

//initialise the network after construction functions are called.
template <typename T>
void NeuralNetwork<T>::initialiseNetwork(){
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
template <typename T>
void NeuralNetwork<T>::constrInputLayer(int numStates){

  feedIns.resize(1, Eigen::VectorXd::Zero(numStates));
  Layers.addInputLayer(feedIns,numStates);
  numLayers++;
}

//---------------------------------------------------------------------------------------------------//
template <typename T>
void NeuralNetwork<T>::constrDenseLayer(
    std::vector<Eigen::VectorXd> const &inputs_, std::string actFunc_,
    int size_
    ){
  Layers.addDenseLayer(inputs_,actFunc_,size_);
  numLayers++;
}

//---------------------------------------------------------------------------------------------------//

template <typename T>
void NeuralNetwork<T>::constrConvLayer(
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

template <typename T>
coeffType NeuralNetwork<T>::getCoeff(detType const &det) const{
	//Run the network
	Eigen::VectorXd output;
#pragma omp critical
	{
  output=feedForward(det);
	}
	// and extract the coefficient from the last layer
  //std::cout << "Nnw.cxx:getCoeff(): output= " << output << std::endl;
	return coeffType(output(0),output(1));
}

//---------------------------------------------------------------------------//

template <typename T>
VecType NeuralNetwork<T>::feedForward(detType const& det) const{
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

template <typename T>
VecType NeuralNetwork<T>::backPropagate(
       coeffType const &lastLayerFeedBack
     ){
  //everytime the backPropagate is called, we should reset nabla* to zero.
  nablaNNP *= 0.;
  // Does not work with empty networks
  if(Layers.size()==0) throw EmptyNetworkError();
  // use the conjugate of lastLayerFeedBack
  Layers[numLayers-1]->backProp(std::conj(lastLayerFeedBack));
  for (size_t layer(numLayers-2); layer > 0; layer--){
    Layers[layer]->backProp(Layers[layer+1]->getDeltas(),
                            Layers[layer+1]->getWeights());
  }

  return nablaNNP;
}

//---------------------------------------------------------------------------//

template <typename T>
Eigen::VectorXd NeuralNetwork<T>::calcNablaPars(
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

template <typename T>
VecType NeuralNetwork<T>::calcNablaParsConnected(
	   State const &inputState,
	   nablaType const &dEdC
     ){
  int numDets = inputState.size();
  Eigen::VectorXd deltaNNP(Eigen::VectorXd::Zero(numNNP));
  Eigen::VectorXcd deltaNNPc(numNNP);
  coeffType realMask(1.,0.);
  coeffType imagMask(0.,1.);
  int pos = 0;

  Eigen::VectorXd deltaNNPTmp(Eigen::VectorXd::Zero(numNNP));
  for (int i(0); i < numDets; ++i){
    feedForward(inputState.det(i));
    deltaNNPc.real() = backPropagate(realMask) ;
    deltaNNPc.imag() = -backPropagate(imagMask) ;
    deltaNNPc *=  dEdC[i].real();
    //deltaNNPc /= std::conj(inputState.coeff(i));
    deltaNNPTmp += 2 * deltaNNPc.real();
    std::vector<detType> coupledDets = inputState.coupledDets(i);
    std::vector<coeffType > coupledCoeffs = inputState.coupledCoeffs(i);
    int pos = inputState.locate(i);
    for (size_t j(0); j < coupledDets.size(); ++j){
      feedForward(coupledDets[j]);
      deltaNNPc.real() = backPropagate(realMask) ;
      deltaNNPc.imag() = backPropagate(imagMask) ;
      deltaNNPc *=  dEdC[numDets+pos+j].real();
      //deltaNNPc /= inputState.coeff(i);
      deltaNNPTmp += 2 * deltaNNPc.real();
    }
  }
  deltaNNP += deltaNNPTmp;
  //deltaNNP /= numDets;
  return deltaNNP;
}

//---------------------------------------------------------------------------------------------------//
/*
Eigen::MatrixXcd NeuralNetwork<T>::calcdCdwSR(
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
    std::vector<coeffType> dEdC(numDets,dedc);
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
*/

//instantiate class
template class NeuralNetwork<VecType>;
}
//---------------------------------------------------------------------------------------------------//

