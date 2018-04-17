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


NeuralNetwork::NeuralNetwork(Hamiltonian const &H_, 
  CostFunction const &externalCF):H(H_), cf(&externalCF), sl(Solver(0.5)){

  //initial value for NNW para
  numLayers = 0;

  //inital para for ADAm
  beta1=0.9;
  beta2=0.999;
  //initial value for momentum algo
  momentumDamping = 0.6;
  momentum = false;

  //initial para for Nesterov's accelerated gradient descent
  lamdaS1 = 0.;
  lamdaS =0.;
  gammaS = 0.;
  gammaS1 = 0.;

}

NeuralNetwork::~NeuralNetwork(){

}
//initialise the network after construction functions are called.
void NeuralNetwork::initialiseNetwork(){
  //calculate the size of the NNP array:
  numNNP=0;
  for (int layer(0); layer<numLayers; ++layer){
    numNNP+=Layers[layer]->getNumPara();
  }
  NNP = Eigen::VectorXd::Ones(numNNP);
  adNNP = &NNP(0);
  nablaNNP = Eigen::VectorXd::Zero(numNNP);
  //get the address of nablaNNP
  adNablaNNP = &nablaNNP(0);
  int startPoint(0);
  for (int layer(0); layer<numLayers; ++layer){
    Layers[layer]->mapPara(adNNP, adNablaNNP, startPoint);
  }

  //-----------------------------------

  generlisedForcePrev = Eigen::VectorXd::Zero(numNNP);
  //vectors for RMSprop
  yS = Eigen::VectorXd::Zero(numNNP);
  yS1 = Eigen::VectorXd::Zero(numNNP);
  Egz2 = Eigen::VectorXd::Zero(numNNP);
  m  = Eigen::VectorXd::Zero(numNNP);
  m1  = Eigen::VectorXd::Zero(numNNP);
  v  = Eigen::VectorXd::Zero(numNNP);
  v1  = Eigen::VectorXd::Zero(numNNP);


  //vectors for ADAM
  m  = Eigen::VectorXd::Zero(numNNP);
  m1  = Eigen::VectorXd::Zero(numNNP);
  v  = Eigen::VectorXd::Zero(numNNP);
  v1  = Eigen::VectorXd::Zero(numNNP);

}
//---------------------------------------------------------------------------//
// construction function of the NNW
void NeuralNetwork::constrInputLayer(int numStates){

  feedIns.resize(1, Eigen::VectorXd::Zero(numStates));
  Layers.push_back(new InputLayer(feedIns, numStates));
  numLayers++;
}

void NeuralNetwork::constrDenseLayer(
    std::vector<Eigen::VectorXd> const &inputs_, std::string actFunc_,
    int size_
    ){
  Layers.push_back(new DenseLayer(inputs_, actFunc_, size_));
  numLayers++;
}


void NeuralNetwork::constrConvLayer(
    std::vector<Eigen::VectorXd> const &inputs_,
    std::string actFunc_,
    int numFilters_,
    int lengthFilter_,
    int stride_
    ){
 Layers.push_back(new ConvLayer(inputs_, actFunc_, numFilters_, lengthFilter_, 
  stride_));
 numLayers++;
}


coeffType NeuralNetwork::getCoeff(detType const &det) const{
	//Run the network
  Eigen::VectorXd output=feedForward(det);
	// and extract the coefficient from the last layer
  //std::cout << "Nnw.cxx:getCoeff(): output= " << output << std::endl;
	return coeffType(output(0),output(1));
}


void NeuralNetwork::updateParameters(
    int method, std::vector<State> const &outputState, double learningRate, 
    int iteration
    ){
  //method corresponds to
  // 0: Stochastic gradiend desend
  // 1: Stochastic reconfiguration
  // 2: Nesterov's Accelerated Gradient Descent
  // 3: ADAM
  Eigen::VectorXd generlisedForce=calcNablaNNP(outputState);
  if (method == 0){
    sl.update(NNP,generlisedForce);
    //update weights and biases
  }
  //Stochastic reconfiguration
  else if (method == 1){
    Eigen::MatrixXcd dCdw=calcdCdwSR(outputState);
    Eigen::VectorXcd ci(outputState.size());
    for (size_t i(0); i<outputState.size(); ++i) ci(i) = outputState[i].coeff;
    sl.update(NNP,generlisedForce,ci,dCdw, iteration);
  }
  else if (method ==2){
     lamdaS1 = (1+std::sqrt(1+4*lamdaS*lamdaS))/2.;
     gammaS = (1-lamdaS)/(lamdaS1);//*std::exp(-1./100*iteration);
     double rho=0.9;//*std::exp(-1./100*iteration);
     Egz2 = rho*Egz2.matrix()+(1-rho)*generlisedForce.array().square().matrix();
     Egz2 += Eigen::VectorXd::Ones(numNNP)*1e-4;
     Eigen::VectorXd RMS = (Egz2).array().sqrt();
     Eigen::VectorXd tau = learningRate * RMS.array().inverse();
     //yS1 = NNP - learningRate * generlisedForce;
     yS1 = (NNP.array() -  tau.array() * generlisedForce.array()).matrix();
     if (iteration == 0) yS=yS1;
     NNP = (1-gammaS)*yS1 + gammaS*yS;
     lamdaS = lamdaS1;
     yS = yS1; 
  }
  else if (method == 3){
    m1 = beta1*m + (1-beta1)*generlisedForce; 
    v1 = beta2*v + (1-beta2)*generlisedForce.array().square().matrix(); 
    m = m1;
    v = v1;
    m1 = m1/(1-std::pow(beta1,iteration+1));
    v1 = v1/(1-std::pow(beta2,iteration+1));
    NNP -= learningRate * (m1.array()/(v1.array().sqrt()+1e-8)).matrix();
  }
  else if (method == 4){
    if (iteration==0) generlisedForcePrev = generlisedForce;
    double beta = (
                    generlisedForce.transpose()
                    *(generlisedForce-generlisedForcePrev)
                    /(generlisedForcePrev.transpose()*generlisedForcePrev)
                  )(0);
    generlisedForce += beta*generlisedForcePrev; 
    //Eigen::MatrixXcd dCdw=calcdCdwSRSR(inputSignalsEpochs, activationsEpochs, outputState);
    //std::vector<coeffType> outputCs = outputState.getAllCoeff();
    //Eigen::Map<Eigen::VectorXcd> ci(&(outputCs[0]),outputCs.size());
    //Eigen::Map<Eigen::VectorXcd> ci(&(outputState[0].coeff),outputState.size());
    Eigen::VectorXcd ci(outputState.size());
    for (size_t i(0); i<outputState.size(); ++i) ci(i) = outputState[i].coeff;
    NNP -= learningRate * generlisedForce;
    generlisedForcePrev = generlisedForce;
  }
}

//---------------------------------------------------------------------------//

Eigen::VectorXd NeuralNetwork::feedForward(detType const& det) const {
  Layers[0]->processSignal(det);
  //std::cout << "Acts layer " << 0 << " =" << std::endl;
  //std::cout<< Layers[0]->getActs()[0] << std::endl;
  //0th layer has no weights!!! The getWeights() for InputLayer is not defined..
  //FixMe, return a 0 vector and throw error.
  //std::cout << "Weights layer " << 0 << " =" << std::endl;
  //std::cout<< Layers[0]->getWeights()[0] << std::endl;
  for (int layer(1); layer < numLayers; ++layer){
      Layers[layer]->processSignal();
      //for (size_t nf(0); nf<Layers[layer]->getActs().size(); nf++){
	//	  std::cout << "Acts layer " << layer << " =" << std::endl;
    //	  std::cout<< "nf=" <<nf << "\n" << Layers[layer]->getActs()[nf] << std::endl;
      //	  for (size_t d(0); d<Layers[layer]->getWeights()[nf].size(); d++){
      //		  std::cout << "Weights layer " << layer << " nf= " << nf << " d=" << d << "\n" << std::endl;
      //		  std::cout<< Layers[layer]->getWeights()[nf][d] << std::endl;
      //	  }
      //	  std::cout << "Nnw.cxx: Biases nf=" << nf << "\n" << std::endl;
   	//	  std::cout<< Layers[layer]->getBiases()[nf] << std::endl;
      //}
  }

  return Layers[numLayers-1]->getActs()[0];
}

//---------------------------------------------------------------------------//

Eigen::VectorXd NeuralNetwork::backPropagate(
       Eigen::VectorXd const &lastLayerFeedBack
     ){
  //everytime the backPropagate is called, we should reset nabla* to zero.
  nablaNNP *= 0.;
  Layers[numLayers-1]->backProp(lastLayerFeedBack);
  for (size_t layer(numLayers-2); layer > 0; layer--){
    Layers[layer]->backProp(Layers[layer+1]->getDeltas(),
                            Layers[layer+1]->getWeights());
  }

  return nablaNNP;
}

//---------------------------------------------------------------------------//

Eigen::VectorXd NeuralNetwork::calcNablaNNP(
	   std::vector<State > const &outputState
     ){
  int numDets = outputState.size();
  //everytime the backPropagate is called, we should reset nabla* to zero.
  Eigen::VectorXd deltaNNP(Eigen::VectorXd::Zero(numNNP));
  std::vector<Eigen::VectorXd> dEdC = cf->nabla(outputState);
  for (int epoch=0; epoch < numDets; ++epoch){
    // obtain inputSignals and activations of all layers
    feedForward(outputState[epoch].det);
    // calculate the derivatives of this determinant
    deltaNNP += backPropagate(dEdC[epoch]);
  }
  return deltaNNP;
}

//---------------------------------------------------------------------------------------------------//

Eigen::VectorXd NeuralNetwork::calcNablaNNPMk(
	   std::vector<State > const &outputState
     ){
  int numDets = outputState.size();
  Eigen::VectorXd deltaNNP(Eigen::VectorXd::Zero(numNNP));
  Eigen::VectorXcd deltaNNPc(numNNP);
  // initialise it to 0!!!
  std::vector<Eigen::VectorXd> dEdC = cf->nabla(outputState);
  Eigen::Vector2d realMask;
  realMask << 1, 0;
  Eigen::Vector2d imagMask;
  imagMask << 0, 1;
  int pos = 0;

  Eigen::VectorXd deltaNNPTmpPrev(Eigen::VectorXd::Zero(numNNP));
  for (int epoch(0); epoch < numDets; ++epoch){
    Eigen::VectorXd deltaNNPTmp(Eigen::VectorXd::Zero(numNNP));
    if (epoch == 0 || (outputState[epoch].det != outputState[epoch-1].det)){  
      feedForward(outputState[epoch].det);
      deltaNNPc.real() = backPropagate(realMask) ;
      deltaNNPc.imag() = -backPropagate(imagMask) ;
      deltaNNPc *=  dEdC[pos](0);
      deltaNNPc /= std::conj(outputState[epoch].coeff); 
      deltaNNPTmp += 2 * deltaNNPc.real();
      pos++;
      std::vector<detType> coupledDets = outputState[epoch].coupledDets;
      std::vector<coeffType > coupledCoeffs = outputState[epoch].coupledCoeffs;
      for (size_t i(0); i < coupledDets.size(); ++i){
        feedForward(coupledDets[i]);
        deltaNNPc.real() = backPropagate(realMask) ;
        deltaNNPc.imag() = backPropagate(imagMask) ;
        deltaNNPc *=  dEdC[pos](0);
        deltaNNPc /= outputState[epoch].coeff; 
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
  std::vector<State > const &outputState
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
      feedForward(outputState[i].det);
      backPropagate(dEdC[epoch]);
      //fill up the dCdwTmp matrix column by column. Column index is C_i, row index is w_k
      dCdwTmp.col(epoch) << nablaNNP; 
    }
    if (i==0) dCdw.real()=dCdwTmp;
    if (i==1) dCdw.imag()=dCdwTmp;
  }
  return dCdw;
}

//---------------------------------------------------------------------------------------------------//

