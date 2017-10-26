/*
 * Nnw.cxx
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */
#include <vector>
#include <iostream>
#include <math.h>
#include <random>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <complex>
#include "Nnw.hpp"
#include "errors.hpp"
//using namespace Eigen;
NeuralNetwork::NeuralNetwork(std::vector<int> const &sizes_, Hamiltonian const&H_,
  Basis const&fullBasis_):sizes(sizes_), H(H_), fullBasis(fullBasis_),normalizerCoeff(0.0){
  energy = 0;
  sampleEnergy = 0.;
  momentumDamping = 0.5;
  momentum = true;
  int numLayersBiasesWeights = sizes.size()-1;
  activations.push_back(Eigen::VectorXd::Zero(sizes[0]));
  inputSignal.push_back(Eigen::VectorXd::Zero(sizes[0]));
  //listDetsToTrain.push_back(fullBasis.getDetByIndex(0));
  for (int layer=0; layer < numLayersBiasesWeights; ++layer){
    //0th layer for biases and weights starting from the 1st layer of neurons;
    //e.g biases[0] represents the biases of neurons on the 1st layer and
    // weights[0] represents the weights of connections between the 0th and 1st
    // layers of neurons.
    activations.push_back(Eigen::VectorXd::Zero(sizes[layer+1]));
    inputSignal.push_back(Eigen::VectorXd::Zero(sizes[layer+1]));
    
    biases.push_back(Eigen::VectorXd::Random(sizes[layer+1]));
    nablaBiases.push_back(Eigen::VectorXd::Zero(sizes[layer+1]));
    gFactorBiases.push_back(Eigen::VectorXd::Ones(sizes[layer+1]));
    gFactorBiasesPrev.push_back(Eigen::VectorXd::Ones(sizes[layer+1]));
    //Pay special attention to weights, it has sizes.size()-1 layers,
    //instead of sizes.size() layers. Especially when reference to which
    //layer, should be careful.
    weights.push_back(Eigen::MatrixXd::Random(sizes[layer+1], sizes[layer]));
    //weights.push_back(1./sizes[layer]*MatrixXd::Random(sizes[layer+1], sizes[layer]));
    nablaWeights.push_back(Eigen::MatrixXd::Zero(sizes[layer+1], sizes[layer]));
    gFactorWeights.push_back(Eigen::MatrixXd::Ones(sizes[layer+1], sizes[layer]));
    gFactorWeightsPrev.push_back(Eigen::MatrixXd::Ones(sizes[layer+1], sizes[layer]));
  }
  nablaWeightsPrev = nablaWeights;
  nablaBiasesPrev = nablaBiases;
}

std::vector<detType> NeuralNetwork::train(std::vector<detType> const &listDetsToTrain, double eta){
  //clear the outputCs vector for each new training
  outputCs.clear();
  //listDetsToTrain = listDetsToTrain_; 
  int numDets(listDetsToTrain.size());
  //std::vector<detType> listDetsToTrainFiltered;
  std::vector<std::vector<Eigen::VectorXd>> inputSignal_Epochs;
  std::vector<std::vector<Eigen::VectorXd>> activations_Epochs;
  std::vector<detType> listDetsToTrainFiltered;
  int numLayersNeuron(sizes.size());
  int numLayersBiasesWeights(sizes.size()-1);
  //reset nablaWeights and nablaBiases to 0 values;
  for (int layer=0; layer < numLayersBiasesWeights; ++layer){
    nablaBiases[layer] *= 0.; 
    nablaWeights[layer] *=0.;
  } 
  double prandom{0.0};
  std::random_device rng;
  double const normalizer = static_cast<double>(rng.max());
  double prob{1.0};
  double coef(0.);
  double coefPrev = feedForward(listDetsToTrain[0])(0);
  inputSignal_Epochs.push_back(inputSignal);
  activations_Epochs.push_back(activations);
  outputCs.push_back(activations[numLayersNeuron-1]);
  listDetsToTrainFiltered.push_back(listDetsToTrain[0]);
  for (int epoch=1; epoch < numDets; ++epoch){
    //initial input layer
    prandom=rng()/normalizer;
    coef=feedForward(listDetsToTrain[epoch])(0);
    prob=fabs(coef/coefPrev);
    coefPrev = coef;
    //if ((prob-1)>1.e-8 || (prandom-prob) < -1.e-8){
    if (true){
      inputSignal_Epochs.push_back(inputSignal);
      activations_Epochs.push_back(activations);
      outputCs.push_back(activations[numLayersNeuron-1]);
      listDetsToTrainFiltered.push_back(listDetsToTrain[epoch]);
    } 
  }
  
  //calculating variational energy
  energy = calcEnergy(listDetsToTrainFiltered);
  sampleEnergy = calcSampleEnergy(listDetsToTrainFiltered);
  //calcLocalEnergy(listDetsToTrainFiltered); 
  nablaWeightsPrev = nablaWeights;
  nablaBiasesPrev = nablaBiases;
  //gFactorWeightsPrev = gFactorWeights;
  //gFactorBiasesPrev = gFactorBiases;
  backPropagate(listDetsToTrainFiltered, inputSignal_Epochs, activations_Epochs); 
  //update weights and biases
  //0th layer is the input layer,
  //the biases should not change.
  //weights have only 0 -- numLayers-2 layers.
  
  for (int layer=0; layer < numLayersBiasesWeights; ++layer){
    if(momentum){
      gFactorBiases[layer] = ((nablaBiases[layer].array() 
                            * nablaBiasesPrev[layer].array()) > 1e-8).select(
                            (gFactorBiasesPrev[layer].array()+0.05).matrix(), gFactorBiasesPrev[layer]*0.95);
      gFactorWeights[layer] = ((nablaWeights[layer].array()
                       * nablaWeightsPrev[layer].array()) > 1e-8).select(
                         (gFactorWeightsPrev[layer].array()+0.05).matrix(), gFactorWeightsPrev[layer]* 0.95);
      nablaBiases[layer] = (nablaBiases[layer].array() * gFactorBiases[layer].array()).matrix();
      nablaWeights[layer] = (nablaWeights[layer].array() * gFactorWeights[layer].array()).matrix();
    }
    nablaBiases[layer] = -eta/numDets  * nablaBiases[layer];
    nablaWeights[layer] = -eta/numDets  * nablaWeights[layer]; 
    if (momentum){
      nablaBiases[layer] += momentumDamping * nablaBiasesPrev[layer];
      nablaWeights[layer] += momentumDamping * nablaWeightsPrev[layer]; 
    }
    biases[layer] += nablaBiases[layer];
    weights[layer] += nablaWeights[layer];
  }
  //double HFCoeff(outputCs.back());
  double probAmp(0.);
  double max(0);
  for (size_t i=0; i < outputCs.size(); ++i){
    probAmp = sqrt(pow(outputCs[i][0],2)+pow(outputCs[i][1],2)); //probility amplitude
    if (fabs(probAmp) -fabs(max) > 1e-8) max = probAmp;
  }
  std::vector<detType> seeds;
  for (size_t i=0; i < outputCs.size(); ++i){
    probAmp = sqrt(pow(outputCs[i][0],2)+pow(outputCs[i][1],2)); //probility amplitude
    if (fabs(probAmp)-0.1*fabs(max) > 1e-8){
      seeds.push_back(listDetsToTrain[i]);
    }
  }
  return seeds;
}

Eigen::VectorXd NeuralNetwork::feedForward(detType const& det){
  int numStates=det.size();
  int numLayersNeuron(sizes.size());
  for (int state=0; state<numStates; ++state){
    //std::cout << "det[" << state << "]=" << det[state] << std::endl;
    activations[0][state] = det[state]?1.0:0.0;
    inputSignal[0][state] = det[state]?1.0:0.0;
  }

  for (int layer=1; layer < numLayersNeuron; ++layer){
    //Here the 0th layer of weights correspond to the connections between
    //the 0th and 1st layer. Here layer refers to the Neuron layer. When use 
    //it to refer the Biases and weights' layer, we need conversion.
    activations[layer] = weights[layer-1]*activations[layer-1]+biases[layer-1];
    inputSignal[layer] = activations[layer];
    if (layer == numLayersNeuron-1)
    activations[layer] = activations[layer].unaryExpr(&Tanh);
    else 
    activations[layer] = activations[layer].unaryExpr(&Tanh);
  }
  return activations[numLayersNeuron-1];
}

void NeuralNetwork::backPropagate(
       std::vector<detType> const& listDetsToTrain,
       std::vector<std::vector<Eigen::VectorXd>> const &inputSignal_Epochs,
       std::vector<std::vector<Eigen::VectorXd>> const &activations_Epochs
     ){
  int numDets = listDetsToTrain.size();
  int numLayersNeuron(sizes.size());
  int numLayersBiasesWeights(sizes.size()-1);
  // catch here
  std::vector<Eigen::VectorXd> dEdC(NablaE_C(listDetsToTrain));
  for (int epoch=0; epoch < numDets; ++epoch){
    Eigen::VectorXd deltaTheLastLayer;
    //std::cout <<"num Det= " << epoch<< std::endl;
    //std::cout << "inputsignal= "<<inputSignal_Epochs[epoch][numLayersNeuron-1][0] << std::endl; 
    //std::cout << "inputsignal size= "<<inputSignal_Epochs.size() << std::endl; 
    //std::cout << "activation= " <<activations_Epochs[epoch][numLayersNeuron-1][0] << std::endl; 
    //std::cout << "activation size= " <<activations_Epochs.size() << std::endl; 
    //std::cout << "Output coeff= " << outputCs[epoch] << std::endl;
    //std::cout << "Output coeff size= " << outputCs.size() << std::endl;
    deltaTheLastLayer = 
    (dEdC[epoch].array() * 
    inputSignal_Epochs[epoch][numLayersNeuron-1].unaryExpr(&Tanh_prime).array()).matrix();
    //adding up all nablaBiases and nablaWeights, in the end use the average
    //of them to determine the final change of the weights and biases.
    //Treat them as Vectors and Matrices.
    //Starting from the last layer of neuron and backpropagate the error.
    //nablaWeights have the same structure as weights, only numLayers-2 layers.
    std::vector<Eigen::VectorXd> deltaAllLayers;
    deltaAllLayers.push_back(deltaTheLastLayer);
    nablaBiases[numLayersBiasesWeights-1] += deltaTheLastLayer;
    nablaWeights[numLayersBiasesWeights-1] += 
      deltaTheLastLayer * activations_Epochs[epoch][numLayersNeuron-2].transpose();
    for (int layer=numLayersBiasesWeights-2; layer >= 0; --layer){
      //Calculating the error from the second last layer to
      //the 1st layer (0th layer is the input neurons, have no biases.)
      //Remember weights have only 0th -- numLayers-2 layers.
      Eigen::VectorXd deltaPreviousLayer=deltaAllLayers[numLayersBiasesWeights-2-layer];
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

      //check if it right
      deltaThisLayer = deltaThisLayer.array()
         * inputSignal_Epochs[epoch][layer+1].unaryExpr(&Tanh_prime).array();
      //deltaThisLayer = deltaAsDiagonal * 
      // VectorXd::Ones(deltaThisLayer.size());
      deltaAllLayers.push_back(deltaThisLayer);
      nablaBiases[layer] += deltaThisLayer;
      //get a weight matrix. \partial C/\partial w^{l}_{jk} = a^{l-1}_k \delta_j^l 
      //the layer here refers to the lth layer of Biases and weights, so for
      //activation layer refers to the l-1th layer.
      nablaWeights[layer] += deltaThisLayer * activations_Epochs[epoch][layer].transpose();
    }
  }
}

double NeuralNetwork::calcSampleEnergy(std::vector<detType> const&listDetsToTrain) {
  double energyVal{0.0};
  normalizerCoeff=0.;
  std::complex<double> normalizerCoeffComplex(0.,0.);
  double Hij(0.);
  int numDets = outputCs.size();
  std::vector<detType> coupledList;
  for (int i=0; i < numDets; ++i){
    std::complex<double> c_i(outputCs[i][0], outputCs[i][1]);
    normalizerCoeffComplex += fabs (std::conj(c_i)  * c_i);
    //sign_i = (outputCs[i]-0. < 1e-8)?-1:0; 
    coupledList.clear();
    coupledList=getCoupledStates(listDetsToTrain[i]);
    for (int j=0; j < coupledList.size(); ++j){
    //for (int j=0; j < numDets; ++j){
      Eigen::VectorXd outputCoupled=feedForward(coupledList[j]);
      std::complex<double> c_j(outputCoupled[0], outputCoupled[1]);
      //std::complex<double> c_j(outputCs[j][0], outputCs[j][1]);
      //std::cout << "j=" << j << " C_j=" << c_j << std::endl;
      Hij = H(listDetsToTrain[i], coupledList[j]);
      energyVal += std::real(std::conj(c_i) * c_j * Hij);
    }
  }
  //std::cout << "normE= " << normalizerCoeff << std::endl;
  normalizerCoeff = std::real(normalizerCoeffComplex);
  energyVal /= normalizerCoeff;
  return energyVal;
}

double NeuralNetwork::calcEnergy(std::vector<detType> const&listDetsToTrain) {
  double energyVal{0.0};
  normalizerCoeff=0.;
  std::complex<double> normalizerCoeffComplex(0.,0.);
  double Hij(0.);
  int numDets = outputCs.size();
  double real(0.), imag(0.);
  std::vector<detType> coupledList;
  for (int i=0; i < numDets; ++i){
    std::complex<double> c_i(outputCs[i][0], outputCs[i][1]);
    normalizerCoeffComplex += fabs (std::conj(c_i)  * c_i);
    //sign_i = (outputCs[i]-0. < 1e-8)?-1:0; 
    coupledList.clear();
    coupledList=getCoupledStates(listDetsToTrain[i]);
    //for (int j=0; j < coupledList.size(); ++j){
    for (int j=0; j < numDets; ++j){
      //Eigen::VectorXd outputCoupled=feedForward(coupledList[j]);
      //std::complex<double> c_j(outputCoupled[0], outputCoupled[1]);
      std::complex<double> c_j(outputCs[j][0], outputCs[j][1]);
      //std::cout << "j=" << j << " C_j=" << c_j << std::endl;
      Hij = H(listDetsToTrain[i], listDetsToTrain[j]);
      energyVal += std::real(std::conj(c_i) * c_j * Hij);
    }
  }
  //std::cout << "normE= " << normalizerCoeff << std::endl;
  normalizerCoeff = std::real(normalizerCoeffComplex);
  energyVal /= normalizerCoeff;
  return energyVal;
}

std::vector<Eigen::VectorXd> NeuralNetwork::NablaE_C(
                     std::vector<detType> const&listDetsToTrain                        
                    ){
  std::vector<Eigen::VectorXd> dEdC;
  int numDets = listDetsToTrain.size();
  // If the input list does have a different size than the number of trained coefficients,
  // throw a corresponding error
  if(numDets != static_cast<int>(outputCs.size())) throw sizeMismatchError(numDets,outputCs.size());
  std::vector<detType> coupledList;
  for (int i=0; i < numDets; ++i){
    Eigen::Vector2d dEdC_i=Eigen::Vector2d::Zero();
    std::complex<double> A(0.,0.);
    coupledList.clear();
    coupledList=getCoupledStates(listDetsToTrain[i]);
    //for (int j=0; j < coupledList.size(); ++j){
    for (int j=0; j < numDets; ++j){
      //Eigen::VectorXd outputCoupled=feedForward(coupledList[j]);
      //std::complex<double> c_j(outputCoupled[0], outputCoupled[1]);
      std::complex<double> c_j(outputCs[j][0], outputCs[j][1]);
      A += c_j * H(listDetsToTrain[i], 
                        listDetsToTrain[j]);
    }
    std::complex<double> c_i(outputCs[i][0], outputCs[i][1]);
    A -=  energy * c_i;
    A /= normalizerCoeff;
    dEdC_i[0] = 2. * std::real( A*std::conj(1.));
    dEdC_i[1] = 2. * std::real( A*std::conj(1. * ii));
    dEdC.push_back(dEdC_i);
  }
  return dEdC;
}

double Tanh_prime(double in){return 1-tanh(in)*tanh(in);};
double Tanh(double in){return tanh(in);};
double Linear(double in) {return in;};
double Linear_prime(double in){return 1;};
double Gaussian(double in){return exp(-pow(in,2));};
double Gaussian_prime(double in){return -2*in*Gaussian(in);};
double GaussianAntiSym(double in){
  double value(0.);
  if (in > 1e-8) value = 1.-Gaussian(in);
  else if (in < -1e-8) value = Gaussian(in)-1.;
  return value;
}
double GaussianAntiSym_prime(double in){
  double value(0.);
  if (in > 1e-8) value = -Gaussian_prime(in);
  else if (in < -1e-8) value = Gaussian_prime(in);
  return value;
}
double Sigmoid(double in){return 1./(1+exp(-in));};
double Sigmoid_prime(double in){return Sigmoid(in)*(1-Sigmoid(in));};
