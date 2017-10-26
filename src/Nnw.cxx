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
#include "State.hpp"
#include "Nnw.hpp"
#include "errors.hpp"
//using namespace Eigen;
NeuralNetwork::NeuralNetwork(std::vector<int> const &sizes_, CostFunction const &externalCF):sizes(sizes_), cf(&externalCF){
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
  nabla_weightsPrev = nabla_weights;
  nabla_biasesPrev = nabla_biases;

  // Assign the output state and the evaluator to 0
  outputState = State();
  // Initialize the energy with NaN, so we see if there is a problem
  energy = 0.0/0.0;
}

std::vector<detType> NeuralNetwork::train(std::vector<detType> const &listDetsToTrain, double eta){
  //clear the output_Cs vector for each new training
  std::vector<coeffType > output_Cs;
  output_Cs.clear();
  nablaWeightsPrev = nablaWeights;
  nablaBiasesPrev = nablaBiases;
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
  //calcLocalEnergy(listDetsToTrainFiltered); 
  // State generation fails if number of determinants and coefficients are not equal
  try{
    outputState = State(listDetsToTrainFiltered, output_Cs);
  }
  catch(sizeMismatchError const &ex){
	  // If the two are not equal in size, shrink them
	  int mSize = (listDetsToTrainFiltered.size() > output_Cs.size())?output_Cs.size():listDetsToTrain.size();
	  output_Cs.resize(mSize);
	  listDetsToTrainFiltered.resize(mSize);
	  outputState = State(listDetsToTrainFiltered,output_Cs);
  }
  //g_weightsPrev = g_weights;
  //g_biasesPrev = g_biases;
  backPropagate(inputSignal_Epochs, activations_Epochs);
  //update weights and biases
  //0th layer is the input layer,
  //the biases should not change.
  //weights have only 0 -- numLayers-2 layers.
  
  for (int layer=0; layer < numLayersBiasesWeights; ++layer){
    if(momentum){
      gFactorBiases[layer] = ((nablaBiases[layer].array() 
                            * nablaBiasesPrev[layer].array()) > 1e-8).select(
                            (gFactorBiasesPrev[layer].array()*1.05).matrix(), gFactorBiasesPrev[layer]*0.95);
      gFactorWeights[layer] = ((nablaWeights[layer].array()
                       * nablaWeightsPrev[layer].array()) > 1e-8).select(
                         (gFactorWeightsPrev[layer].array()*1.05).matrix(), gFactorWeightsPrev[layer]* 0.95);
      nablaBiases[layer] = (nablaBiases[layer].array() * gFactorBiases[layer].array()).matrix();
      nablaWeights[layer] = (nablaWeights[layer].array() * gFactorWeights[layer].array()).matrix();
    }
    nablaBiases[layer] = -eta / numDets * nablaBiases[layer];
    nablaWeights[layer] = -eta / numDets * nablaWeights[layer]; 
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
    probAmp = pow(outputCs[i][0],2)-pow(outputCs[i][1],2); //probility amplitude
    if (fabs(probAmp) -fabs(max) > 1e-8) max = probAmp;
  }
  std::vector<detType> seeds;
  for (size_t i=0; i < outputCs.size(); ++i){
    probAmp = pow(outputCs[i][0],2)+pow(outputCs[i][1],2); //probility amplitude
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
       std::vector<std::vector<Eigen::VectorXd>> const &inputSignal_Epochs,
       std::vector<std::vector<Eigen::VectorXd>> const &activations_Epochs
     ){
  int numDets = outputState.size();
  int numLayersNeuron(sizes.size());
  int numLayersBiasesWeights(sizes.size()-1);
  // catch here
  Evaluator cfWrapper(*cf, outputState);
  std::vector<Eigen::VectorXd> dEdC(cfWrapper.nabla());
  energy = cfWrapper();
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

double Tanh_prime(double in){return 1-tanh(in)*tanh(in);};
double Tanh(double in){return tanh(in);};
double Linear(double in) {return in;};
double Linear_prime(double in){return 1;};
double Gaussian(double in){return exp(-pow(in,2));};
double Gaussian_prime(double in){return -2*in*Gaussian(in);};
double Sigmoid(double in){return 1./(1+exp(-in));};
double Sigmoid_prime(double in){return Sigmoid(in)*(1-Sigmoid(in));};
