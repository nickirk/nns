/*
 * Nnw.cxx
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */
#include <vector>
#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "Nnw.hpp"
using namespace Eigen;
NeuralNetwork::NeuralNetwork(std::vector<int> const &sizes_, Hamiltonian const&H_,
  Basis const&fullBasis_):sizes(sizes_), H(H_), fullBasis(fullBasis_){
  energy = 0;
  int numLayersNeuron = sizes.size();
  int numLayersBiasesWeights = sizes.size()-1;
  //listDetsToTrain.push_back(fullBasis.getDetByIndex(0));
  for (int layer=0; layer < numLayersBiasesWeights; ++layer){
    //0th layer for biases and weights starting from the 1st layer of neurons;
    //e.g biases[0] represents the biases of neurons on the 1st layer and
    // weights[0] represents the weights of connections between the 0th and 1st
    // layers of neurons.
    biases.push_back(0.1*VectorXd::Random(sizes[layer+1]));
    nabla_biases.push_back(VectorXd::Zero(sizes[layer+1]));
    //Pay special attention to weights, it has sizes.size()-1 layers,
    //instead of sizes.size() layers. Especially when reference to which
    //layer, should be careful.
    weights.push_back(MatrixXd::Random(sizes[layer+1], sizes[layer]));
    //weights.push_back(1./sizes[layer]*MatrixXd::Random(sizes[layer+1], sizes[layer]));
    nabla_weights.push_back(MatrixXd::Zero(sizes[layer+1], sizes[layer]));
  }
}

void NeuralNetwork::calcEnergy(std::vector<detType> const&listDetsToTrain){
  energy = 0.;
  sign = 0;
  int sign_i=0;
  double normalizerCoeff(0.);
  double Hij(0.);
  int numDets = listDetsToTrain.size();
  for (int i=0; i < numDets; ++i){
    normalizerCoeff += output_Cs[i] * output_Cs[i];
    sign_i = (output_Cs[i]-0. < 1e-8)?-1:0; 
    sign +=sign_i;
    for (int j=0; j < numDets; ++j){
      Hij = H(listDetsToTrain[i], listDetsToTrain[j]);
      energy += output_Cs[i] * output_Cs[j] * Hij;
    }
  }
  //std::cout << "normE= " << normalizerCoeff << std::endl;
  energy /= normalizerCoeff;
  std::cout << "energy=" << energy <<std::endl;

}

std::vector<detType> NeuralNetwork::train(std::vector<detType> const &listDetsToTrain, double eta){
  //clear the output_Cs vector for each new training
  output_Cs.clear();
  //listDetsToTrain = listDetsToTrain_; 
  int numDets(listDetsToTrain.size());
  std::vector<std::vector<VectorXd>> inputSignal_Epochs;
  std::vector<std::vector<VectorXd>> activations_Epochs;
  std::vector<VectorXd> inputSignal; 
  std::vector<VectorXd> activations; 
  //reset nabla_weights and nabla_biases to 0 values;
  int numLayersNeuron(sizes.size());
  int numLayersBiasesWeights(sizes.size()-1);
  //initialize activations, nabla_biases and nabla_weights to 0
  activations.push_back(VectorXd::Zero(sizes[0]));
  inputSignal.push_back(VectorXd::Zero(sizes[0]));
  for (int layer=0; layer < numLayersBiasesWeights; ++layer){
    activations.push_back(VectorXd::Zero(sizes[layer+1]));
    inputSignal.push_back(VectorXd::Zero(sizes[layer+1]));
    nabla_biases[layer] *= 0.; 
    nabla_weights[layer] *=0.;
  } 
  for (int epoch=0; epoch < numDets; ++epoch){
    //initial input layer
    int numStates = listDetsToTrain[epoch].size();
    for (int state=0; state<numStates; ++state){
      activations[0][state] = listDetsToTrain[epoch][state]?1.0:0.0;
      inputSignal[0][state] = listDetsToTrain[epoch][state]?1.0:0.0;
    }
    //activations[0] = Eigen::Map<Eigen::VectorXd>(listDetsToTrain[epoch].data());
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
    inputSignal_Epochs.push_back(inputSignal);
    activations_Epochs.push_back(activations);
    output_Cs.push_back(activations[numLayersNeuron-1][0]);
   
    //std::cout <<"1 num Det= " << epoch<< std::endl;
    //std::cout << "inputSignalLastLayer= " << inputSignal[numLayersNeuron-1] <<std::endl;
    //std::cout << "inputsignal_epochLL= "<<inputSignal_Epochs[epoch][numLayersNeuron-1] << std::endl; 
    //std::cout << "C_" << epoch << "= " << output_Cs[epoch] << std::endl;
    //std::cout << "1 inputsignal_epoch size= "<<inputSignal_Epochs.size() << std::endl; 
    ////std::cout << "activations= " << activations[numLayersNeuron-1][0] <<std::endl;
    //std::cout << "1 activation_epoch= " <<activations_Epochs[epoch][numLayersNeuron-1][0] << std::endl; 

    //std::cout << "1 activation_epoch size= " <<activations_Epochs.size() << std::endl; 
    //std::cout << "1 Output coeff= " << output_Cs[epoch] << std::endl;
    //std::cout << "1 Output coeff size= " << output_Cs.size() << std::endl;
    //std::cout << "1 last bias= " << biases[1] << std::endl;
  }
  
  //calculating variational energy
  calcEnergy(listDetsToTrain); 
  backPropagate(listDetsToTrain, inputSignal_Epochs, activations_Epochs); 
  //update weights and biases
  //0th layer is the input layer,
  //the biases should not change.
  //weights have only 0 -- numLayers-2 layers.
  for (int layer=0; layer < numLayersBiasesWeights; ++layer){
    biases[layer] -= eta / numDets * nabla_biases[layer];
    weights[layer] -= eta / numDets * nabla_weights[layer];
  }
  double HFCoeff(output_Cs[0]);
  std::vector<detType> seeds;
  //std::vector<double> dEdC(NablaE_C(listDetsToTrain));
  for (int i=0; i < output_Cs.size(); ++i){
    //std::cout << "HFCoeff= " << HFCoeff << std::endl;
    std::cout << "C_" << i << "= " << output_Cs[i] << std::endl;
    std::cout << "intCast= " << verbatimCast(listDetsToTrain[i]) << std::endl;
    //std::cout << "yes or no" << "= " <<fabs(fabs(output_Cs[i])-0.2*fabs(HFCoeff))  << std::endl;
    if (fabs(output_Cs[i])-0.5*fabs(HFCoeff) > 1e-8){
      //std::cout << "yes" << std::endl;
      seeds.push_back(listDetsToTrain[i]);
    }
  }
  return seeds;
}

//void NeuralNetwork::feedForward();
void NeuralNetwork::backPropagate(
       std::vector<detType> const& listDetsToTrain,
       std::vector<std::vector<VectorXd>> const &inputSignal_Epochs,
       std::vector<std::vector<VectorXd>> const &activations_Epochs
     ){
  int numDets = listDetsToTrain.size();
  int numLayersNeuron(sizes.size());
  int numLayersBiasesWeights(sizes.size()-1);
  std::vector<double> dEdC(NablaE_C(listDetsToTrain));
  for (int epoch=0; epoch < numDets; ++epoch){
    VectorXd deltaTheLastLayer;
    //std::cout <<"num Det= " << epoch<< std::endl;
    //std::cout << "inputsignal= "<<inputSignal_Epochs[epoch][numLayersNeuron-1][0] << std::endl; 
    //std::cout << "inputsignal size= "<<inputSignal_Epochs.size() << std::endl; 
    //std::cout << "activation= " <<activations_Epochs[epoch][numLayersNeuron-1][0] << std::endl; 
    //std::cout << "activation size= " <<activations_Epochs.size() << std::endl; 
    //std::cout << "Output coeff= " << output_Cs[epoch] << std::endl;
    //std::cout << "Output coeff size= " << output_Cs.size() << std::endl;
    deltaTheLastLayer = 
    dEdC[epoch]*inputSignal_Epochs[epoch][numLayersNeuron-1].unaryExpr(
                                    &Tanh_prime);
    //adding up all nabla_biases and nabla_weights, in the end use the average
    //of them to determine the final change of the weights and biases.
    //Treat them as Vectors and Matrices.
    //Starting from the last layer of neuron and backpropagate the error.
    //nabla_weights have the same structure as weights, only numLayers-2 layers.
    std::vector<VectorXd> deltaAllLayers;
    deltaAllLayers.push_back(deltaTheLastLayer);
    nabla_biases[numLayersBiasesWeights-1] += deltaTheLastLayer;
    nabla_weights[numLayersBiasesWeights-1] += 
      deltaTheLastLayer * activations_Epochs[epoch][numLayersNeuron-2].transpose();
    for (int layer=numLayersBiasesWeights-2; layer >= 0; --layer){
      //Calculating the error from the second last layer to
      //the 1st layer (0th layer is the input neurons, have no biases.)
      //Remember weights have only 0th -- numLayers-2 layers.
      VectorXd deltaPreviousLayer=deltaAllLayers[numLayersBiasesWeights-2-layer];
      VectorXd deltaThisLayer;
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

std::vector<double> NeuralNetwork::NablaE_C(
                     std::vector<detType> const&listDetsToTrain                        
                    ){
  std::vector<double> dEdC;
  int numDets = listDetsToTrain.size();
  for (int i=0; i < numDets; ++i){
    double dEdC_i(0.);
    double normalizerCoeff(0.);
    for (int j=0; j < numDets; ++j){
      normalizerCoeff += output_Cs[j]*output_Cs[j];
      dEdC_i += 2 * output_Cs[j] * H(listDetsToTrain[i], 
                                                   listDetsToTrain[j]);
    }
    dEdC_i -= 2*energy*output_Cs[i];
    dEdC_i /= normalizerCoeff;
    dEdC.push_back(dEdC_i);
  }
  return dEdC;
}

double Tanh_prime(double in){return 1-tanh(in)*tanh(in);};
double Tanh(double in){return tanh(in);};
double Linear(double in) {return in;};
double Linear_prime(double in){return 1;};
