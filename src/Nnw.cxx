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
  momentumDamping = 0.8;
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
    nabla_biases.push_back(Eigen::VectorXd::Zero(sizes[layer+1]));
    g_biases.push_back(Eigen::VectorXd::Ones(sizes[layer+1]));
    g_biasesPrev.push_back(Eigen::VectorXd::Ones(sizes[layer+1]));
    //Pay special attention to weights, it has sizes.size()-1 layers,
    //instead of sizes.size() layers. Especially when reference to which
    //layer, should be careful.
    weights.push_back(Eigen::MatrixXd::Random(sizes[layer+1], sizes[layer]));
    //weights.push_back(1./sizes[layer]*MatrixXd::Random(sizes[layer+1], sizes[layer]));
    nabla_weights.push_back(Eigen::MatrixXd::Zero(sizes[layer+1], sizes[layer]));
    g_weights.push_back(Eigen::MatrixXd::Ones(sizes[layer+1], sizes[layer]));
    g_weightsPrev.push_back(Eigen::MatrixXd::Ones(sizes[layer+1], sizes[layer]));
  }
  nabla_weightsPrev = nabla_weights;
  nabla_biasesPrev = nabla_biases;
}

std::vector<detType> NeuralNetwork::train(std::vector<detType> const &listDetsToTrain, double eta){
  //clear the output_Cs vector for each new training
  output_Cs.clear();
  //listDetsToTrain = listDetsToTrain_; 
  int numDets(listDetsToTrain.size());
  //std::vector<detType> listDetsToTrainFiltered;
  std::vector<std::vector<Eigen::VectorXd>> inputSignal_Epochs;
  std::vector<std::vector<Eigen::VectorXd>> activations_Epochs;
  std::vector<detType> listDetsToTrainFiltered;
  int numLayersNeuron(sizes.size());
  int numLayersBiasesWeights(sizes.size()-1);
  //reset nabla_weights and nabla_biases to 0 values;
  for (int layer=0; layer < numLayersBiasesWeights; ++layer){
    nabla_biases[layer] *= 0.; 
    nabla_weights[layer] *=0.;
  } 
  double prandom{0.0};
  std::random_device rng;
  double const normalizer = static_cast<double>(rng.max());
  double prob{1.0};
  double coef(0.);
  double coefPrev = feedForward(listDetsToTrain[0])(0);
  inputSignal_Epochs.push_back(inputSignal);
  activations_Epochs.push_back(activations);
  output_Cs.push_back(activations[numLayersNeuron-1]);
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
      output_Cs.push_back(activations[numLayersNeuron-1]);
      listDetsToTrainFiltered.push_back(listDetsToTrain[epoch]);
    } 
  }
  
  //calculating variational energy
  energy = calcEnergy(listDetsToTrainFiltered);
  //calcLocalEnergy(listDetsToTrainFiltered); 
  nabla_weightsPrev = nabla_weights;
  nabla_biasesPrev = nabla_biases;
  //g_weightsPrev = g_weights;
  //g_biasesPrev = g_biases;
  backPropagate(listDetsToTrainFiltered, inputSignal_Epochs, activations_Epochs); 
  //update weights and biases
  //0th layer is the input layer,
  //the biases should not change.
  //weights have only 0 -- numLayers-2 layers.
  
  for (int layer=0; layer < numLayersBiasesWeights; ++layer){
    if(momentum){
      g_biases[layer] = ((nabla_biases[layer].array() 
                            * nabla_biasesPrev[layer].array()) > 1e-8).select(
                            (g_biasesPrev[layer].array()*1.05).matrix(), g_biasesPrev[layer]*0.95);
      g_weights[layer] = ((nabla_weights[layer].array()
                       * nabla_weightsPrev[layer].array()) > 1e-8).select(
                         (g_weightsPrev[layer].array()*1.05).matrix(), g_weightsPrev[layer]* 0.95);
      nabla_biases[layer] = (nabla_biases[layer].array() * g_biases[layer].array()).matrix();
      nabla_weights[layer] = (nabla_weights[layer].array() * g_weights[layer].array()).matrix();
    }
    nabla_biases[layer] = -eta / numDets * nabla_biases[layer];
    nabla_weights[layer] = -eta / numDets * nabla_weights[layer]; 
    if (!momentum){
      nabla_biases[layer] += momentumDamping * nabla_biasesPrev[layer];
      nabla_weights[layer] += momentumDamping * nabla_weightsPrev[layer]; 
    }
    biases[layer] += nabla_biases[layer];
    weights[layer] += nabla_weights[layer];
  }
  //double HFCoeff(output_Cs.back());
  double probAmp(0.);
  double max(0);
  for (size_t i=0; i < output_Cs.size(); ++i){
    probAmp = pow(output_Cs[i][0],2)-pow(output_Cs[i][1],2); //probility amplitude
    if (fabs(probAmp) -fabs(max) > 1e-8) max = probAmp;
  }
  std::vector<detType> seeds;
  for (size_t i=0; i < output_Cs.size(); ++i){
    probAmp = pow(output_Cs[i][0],2)+pow(output_Cs[i][1],2); //probility amplitude
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
    //std::cout << weights[layer-1] << std::endl;
    //std::cout << activations[layer-1] << std::endl;
    activations[layer] = weights[layer-1]*activations[layer-1]+biases[layer-1];
    inputSignal[layer] = activations[layer];
    //if (layer == numLayersNeuron-1){
    //  std::cout << "inputSignal act= " << activations[layer] << std::endl;        
    //  std::cout << "inputSignal size= " << inputSignal.size() << std::endl;        
    //  std::cout << "inputSignal= " << inputSignal[layer] << std::endl;        
    //}
    if (layer == numLayersNeuron-1)
    activations[layer] = activations[layer].unaryExpr(&Tanh);
    else 
    activations[layer] = activations[layer].unaryExpr(&Tanh);
    //if (layer == numLayersNeuron-1){
    //  std::cout << "activation= " << activations[layer] << std::endl;        
    //}
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
    //std::cout << "Output coeff= " << output_Cs[epoch] << std::endl;
    //std::cout << "Output coeff size= " << output_Cs.size() << std::endl;
    deltaTheLastLayer = 
    (dEdC[epoch].array() * 
    inputSignal_Epochs[epoch][numLayersNeuron-1].unaryExpr(&Tanh_prime).array()).matrix();
    //adding up all nabla_biases and nabla_weights, in the end use the average
    //of them to determine the final change of the weights and biases.
    //Treat them as Vectors and Matrices.
    //Starting from the last layer of neuron and backpropagate the error.
    //nabla_weights have the same structure as weights, only numLayers-2 layers.
    std::vector<Eigen::VectorXd> deltaAllLayers;
    deltaAllLayers.push_back(deltaTheLastLayer);
    nabla_biases[numLayersBiasesWeights-1] += deltaTheLastLayer;
    nabla_weights[numLayersBiasesWeights-1] += 
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
      nabla_biases[layer] += deltaThisLayer;
      //get a weight matrix. \partial C/\partial w^{l}_{jk} = a^{l-1}_k \delta_j^l 
      //the layer here refers to the lth layer of Biases and weights, so for
      //activation layer refers to the l-1th layer.
      nabla_weights[layer] += deltaThisLayer * activations_Epochs[epoch][layer].transpose();
    }
  }
}

double NeuralNetwork::calcEnergy(std::vector<detType> const&listDetsToTrain) const{
  double energyVal{0.0};
  normalizerCoeff=0.;
  std::complex<double> normalizerCoeffComplex(0.,0.);
  double Hij(0.);
  int numDets = output_Cs.size();
  double real(0.), imag(0.);
  for (int i=0; i < numDets; ++i){
    real = output_Cs[i][0];
    imag = output_Cs[i][1];
    std::complex<double> c_i(real, imag);
    normalizerCoeffComplex += fabs (std::conj(c_i)  * c_i);
    //sign_i = (output_Cs[i]-0. < 1e-8)?-1:0; 
    for (int j=0; j < numDets; ++j){
      real = output_Cs[j][0];
      imag = output_Cs[j][1];
      std::complex<double> c_j(real, imag);
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
  if(numDets != static_cast<int>(output_Cs.size())) throw sizeMismatchError(numDets,output_Cs.size());
  for (int i=0; i < numDets; ++i){
    Eigen::Vector2d dEdC_i=Eigen::Vector2d::Zero();
    std::complex<double> A(0.,0.);
    for (int j=0; j < numDets; ++j){
      std::complex<double> c_j(output_Cs[j][0], output_Cs[j][1]);
      A += c_j * H(listDetsToTrain[i], 
                        listDetsToTrain[j]);
    }
    std::complex<double> c_i(output_Cs[i][0], output_Cs[i][1]);
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
double Sigmoid(double in){return 1./(1+exp(-in));};
double Sigmoid_prime(double in){return Sigmoid(in)*(1-Sigmoid(in));};
