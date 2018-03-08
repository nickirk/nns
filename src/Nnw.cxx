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
#include "State.hpp"
#include "Nnw.hpp"
#include "utilities/Errors.hpp"
#include "NormCF.hpp"
#include "Solver.hpp"
#include "Layer.hpp"


NeuralNetwork::NeuralNetwork(Hamiltonian const &H_, std::vector<int> const &sizes_, 
CostFunction const &externalCF):H(H_), sizes(sizes_), cf(&externalCF), sl(Solver(0.5)){
  //initial value for NNW para
  numLayers = 0;
  generlisedForcePrev = Eigen::VectorXd::Zero(numNNP);

  //initial value for momentum algo
  momentumDamping = 0.6;
  momentum = false;

  //initial para for Nesterov's accelerated gradient descent
  lamdaS1 = 0.;
  lamdaS =0.;
  gammaS = 0.;
  gammaS1 = 0.;
  //-----------------------------------
  //inital para for ADAm
  beta1=0.9;
  beta2=0.999;

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
  //map parameters to matrices and initalize them with normal distribution.
  //map nablaPara to matrices and vectors with 0.
}
//initialise the network after construction functions are called.
void NeuralNetwork::initialiseNetwork(){
  //calculate the size of the NNP array:
  numNNP=0;
  for (int layer(0); layer<numLayers; ++layer){
    numNNP+=network[layer].numPara;
  }
  NNP = Eigen::VectorXd::Ones(numNNP);
  adNNP = &NNP(0);
  int startPoint(0);
  for (int layer(0); layer<numLayers; ++layer){
    network[layer].mapPara(adNNP, startPoint);
  /*
  //initial activity signals. 
  activations.push_back(Eigen::VectorXd::Zero(sizes[0]));
  inputSignals.push_back(Eigen::VectorXd::Zero(sizes[0]));
  for (int layer=0; layer < numLayers; ++layer){
    activations.push_back(Eigen::VectorXd::Zero(sizes[layer+1]));
    inputSignals.push_back(Eigen::VectorXd::Zero(sizes[layer+1]));
  }
  */
  }

  nablaNNP = Eigen::VectorXd::Zero(numNNP);
  //get the address of nablaNNP
  adNablaNNP = &nablaNNP(0);
}
//---------------------------------------------------------------------------------------------------//
// construction function of the NNW
void NeuralNetwork::constrInputLayer(int numNrn){


}

void NeuralNetwork::constrDenseLayer(std::vector<Eigen::MatrixXd> const 
                                     &inputs_, double &(actFunc_)(double),
                                     int size_){
 DenseLayer denseLayer(inputs_, actFunc_(double), size_);
 network.push_back(denseLayer);
 numLayers++;
}
/*
void NeuralNetwork::constrConvLayer(std::vector<Eigen::MatrixXd> const 
                                     &inputs_,double &(actFunc_)(double),
                                     int size_){
 ConvLayer convLayer(inputs_, actFunc_(double));
 network.push_back(convLayer);
 numLayers++;
}
*/


coeffType NeuralNetwork::getCoeff(detType const &det) const{
	//Run the network
	feedForward(det);
	// and extract the coefficient from the last layer
	return outputLayer();
}


void NeuralNetwork::updateParameters(int method, std::vector<State> const &outputState, double learningRate, int iteration){
  //method corresponds to
  // 0: Stochastic gradiend desend
  // 1: Stochastic reconfiguration
  // 2: Nesterov's Accelerated Gradient Descent
  // 3: ADAM
  Eigen::VectorXd generlisedForce=calcNablaNNPMk(outputState);
  if (method == 0){
    sl.update(NNP,generlisedForce);
    //update weights and biases
  }
  //Stochastic reconfiguration
  else if (method == 1){
    Eigen::MatrixXcd dCdw=calcdCdwSR(outputState);
    Eigen::VectorXcd ci(outputState.size());
    for (size_t i(0); i<outputState.size(); ++i) ci(i) = outputState[i].coeff;
    //Eigen::Map<Eigen::VectorXcd> ci(&(outputState[0].coeff),outputState.size());

    sl.update(NNP,generlisedForce,ci,dCdw, iteration);
  }
  else if (method ==2){
     lamdaS1 = (1+std::sqrt(1+4*lamdaS*lamdaS))/2.;
     gammaS = (1-lamdaS)/(lamdaS1);//*std::exp(-1./100*iteration);
     double rho=0.9;//*std::exp(-1./100*iteration);
     Egz2 = rho*Egz2.matrix() + (1-rho)*generlisedForce.array().square().matrix();
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
    double beta = (generlisedForce.transpose()*(generlisedForce-generlisedForcePrev)/
                  (generlisedForcePrev.transpose()*generlisedForcePrev))(0);
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

//---------------------------------------------------------------------------------------------------//

Eigen::VectorXd NeuralNetwork::feedForward(detType const& det) const{
  int numStates=det.size();
  int numLayersNeuron(sizes.size());
  for (int state=0; state<numStates; ++state){
    activations[0][state] = det[state]?1.0:-1.0;
    inputSignals[0][state] = det[state]?1.0:-1.0;
  }

  for (int layer=1; layer < numLayersNeuron; ++layer){
    //Here the 0th layer of weights correspond to the connections between
    //the 0th and 1st layer. Here layer refers to the Neuron layer. When use 
    //it to refer the Biases and weights' layer, we need conversion.
    activations[layer] = weights[layer-1]*activations[layer-1]+biases[layer-1];
    inputSignals[layer] = activations[layer];
    if (layer == numLayersNeuron-1)
    activations[layer] = activations[layer].unaryExpr(&Linear);
    else 
    activations[layer] = activations[layer].unaryExpr(&Tanh);
  }
  return activations[numLayersNeuron-1];
}

//---------------------------------------------------------------------------------------------------//

Eigen::VectorXd NeuralNetwork::backPropagate(
       Eigen::VectorXd	            const &lastLayerFeedBack
     ){
  int numLayers(sizes.size()-1);
  int numLayersNeuron(sizes.size());
  //everytime the backPropagate is called, we should reset nabla* to zero.
  nablaNNP *= 0.;

  Eigen::VectorXd deltaTheLastLayer;
  deltaTheLastLayer = 
  (lastLayerFeedBack.array() * 
  inputSignals[numLayersNeuron-1].unaryExpr(&Linear_prime).array()).matrix();
  //adding up all nablaBiases and nablaWeights, in the end use the average
  //of them to determine the final change of the weights and biases.
  //Treat them as Vectors and Matrices.
  //Starting from the last layer of neuron and backpropagate the error.
  //nablaWeights have the same structure as weights, only numLayers-2 layers.
  std::vector<Eigen::VectorXd> deltaAllLayers;
  deltaAllLayers.push_back(deltaTheLastLayer);
  nablaBiases[numLayers-1] = deltaTheLastLayer;
  nablaWeights[numLayers-1] = 
    deltaTheLastLayer * activations[numLayersNeuron-2].transpose();
  for (int layer=numLayers-2; layer >= 0; --layer){
    //Calculating the error from the second last layer to
    //the 1st layer (0th layer is the input neurons, have no biases.)
    //Remember weights have only 0th -- numLayers-2 layers.
    Eigen::VectorXd deltaPreviousLayer=deltaAllLayers[numLayers-2-layer];
    Eigen::VectorXd deltaThisLayer;
    //MatrixXd deltaAsDiagonal;
    //\delta_l = ((weights_{l+1, l}^T\delta^{l+1})) .* tanh'(z^l);
    //where ^T means transpose, l means lth layer, tanh' means the derivative
    //of tanh and z^l= tanh(weights_{l,l-1}activations^{l-1}) is the input signal 
    //on the lth layer neurons, where weights_{l, l-1} corresponds to the 
    //connections between the lth and l-1 th layers of neurons. Notice
    // .* means elementwise multiply.
    
    deltaThisLayer = weights[layer+1].transpose() * deltaPreviousLayer; 
    //To achive elementwise multiply, form a diagonal vector.
    //deltaAsDiagonal = deltaThisLayer.asDiagonal();

    deltaThisLayer = deltaThisLayer.array()
       * inputSignals[layer+1].unaryExpr(&Tanh_prime).array();
    deltaAllLayers.push_back(deltaThisLayer);
    nablaBiases[layer] = deltaThisLayer;
    //get a weight matrix. \partial C/\partial w^{l}_{jk} = a^{l-1}_k \delta_j^l 
    //the layer here refers to the lth layer of Biases and weights, so for
    //activation layer refers to the l-1th layer.
    nablaWeights[layer] = deltaThisLayer * activations[layer].transpose();
  }
  return nablaNNP;
}

//---------------------------------------------------------------------------------------------------//

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

