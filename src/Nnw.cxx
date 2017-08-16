/*
 * Nnw.cxx
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */
#include <vector>
#include <stdio.h>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Core>
using namespace Eigen;
NeuralNetwork::NeuralNetwork(std::vector<int> sizes_, Hamiltonian const&H_,
  Basis const&fullBasis_):sizes(sizes_), H(H_), fullBasis(fullBasis_){
  energy = 0;
  int numLayersNeuron = sizes.size();
  int numLayersBaisesWeights = sizes.size()-1;
  for (int layer=0; layer < numLayersNeuron; ++layer){
    activations.push_back(VectorXd::Zero(sizes[i]));
    //0th layer for biases and weights starting from the 1st layer of neurons;
    //e.g biases[0] represents the biases of neurons on the 1st layer and
    // weights[0] represents the weights of connections between the 0th and 1st
    // layers of neurons.
    if (layer < numLayersBaisesWeights) {
      biases.push_back(VectorXd::Random(sizes[layer+1]));
      nabla_biases.push_back(VectorXd::Zero(sizes[layer+1]));
      //Pay special attention to weights, it has sizes.size()-1 layers,
      //instead of sizes.size() layers. Especially when reference to which
      //layer, should be careful.
      weights.push_back(MatrixXd::Random(sizes[layer+1], sizes[layer]));
      nabla_weights.push_back(MatrixXd::Zero(sizes[layer+1], sizes[layer]));
    }
  }
}

void NeuralNetwork::calcEnergy(std::vector<detType> &list){
  energy = 0.;
  double normalizerCoeff(0.);
  for (int i=0; i < list.size(); ++i){
    for (int j=i; j < list.size(); ++j){
      Hij = 2 * H(list[i], list[j]);
      energy += output_Cs[i] * output_Cs[j] * Hij;
      normalizerCoeff += 2* output_Cs[i] * output_Cs[j];
    }
  }
  energy /= normalizerCoeff;
}

void NeuralNetwork::train(std::vector<detType> &list, double eta){
  int numDets(list.size());
  std::vector<vector<VectorXd>> inputSignal_Epochs;
  std::vector<VectorXd> inputSignal; 
  //reset nabla_weights and nabla_biases to 0 values;
  int numLayersNeuron(sizes.size());
  int numLayersBaisesWeights(sizes.size()-1);
  for (int layer=0; layer < numLayersBaisesWeights; ++layer){
    nabla_biases[layer] *= 0.; 
    nabla_weights[layer] *=0.;
  } 
  for (int epoch=0; epoch < numDets; ++epoch){
    //initial input layer
    Map<VectorxXd> activations[0](*list[epoch][0], list[epoch][0].size());
    for (int layer=1; layer < numLayersNeuron; ++layer){
      //Here the 0th layer of weights correspond to the connections between
      //the 0th and 1st layer. Here layer refers to the Neuron layer. When use 
      //it to refer the Biases and weights' layer, we need conversion.
      activations[layer] = weights[layer-1]*activations[layer-1]+biases[layer];
      inputSignal.push_back(activations[layer]);
      activations[layer] = activations[layer].unaryExpr(&Tanh);
    }
    inputSignal_Epochs.push_back(inputSignal);
    output_Cs.push_back(activations[numLayersNeuron-1]);
  }
  
  //calculating variational energy
  calcEnergy(); 
  backPropagate(numDets); 
  //update weights and biases
  //0th layer is the input layer,
  //the biases should not change.
  //weights have only 0 -- numLayers-2 layers.
  for (int layer=0; layer < numLayersBaisesWeights; ++layer){
    biases[layer] -= eta / numDets * nabla_biases[layer];
    weights[layer] -= eta / numDets * nabla_weights[layer];
    }
  } 
    
}

//void NeuralNetwork::feedForward();
void backPropagate(int numDets){
  int numLayersNeuron(sizes.size());
  int numLayersBaisesWeights(sizes.size()-1);
  std::vector<double> dEdC(NablaE_C());
  for (int epoch=0; epoch < numDets; ++i){
    VectorXd deltaTheLastLayer;
    deltaTheLastLayer = 
      dEdC[epoch]*inputSignal_Epochs[epoch][numLayersNeuron-1].unaryExpr(&Tanh_prime);
    //adding up all nabla_biases and nabla_weights, in the end use the average
    //of them to determine the final change of the weights and biases.
    //Treat them as Vectors and Matrices.
    //Starting from the last layer of neuron and backpropagate the error.
    //nabla_weights have the same structure as weights, only numLayers-2 layers.
    nabla_biases[numLayersBaisesWeights-1] += deltaTheLastLayer;
    nabla_weights[numLayersBaisesWeights-1] += activations[numLayersNeuron-2]
      * deltaTheLastLayer.transpose();
    for (int layer=numLayersBaisesWeights-2; layer >= 0; --layer){
      //Calculating the error from the second last layer to
      //the 1st layer (0th layer is the input neurons, have no biases.)
      //Remember weights have only 0th -- numLayers-2 layers.
      VectorXd deltaPreviousLayer(delta_LastLayer);
      VectorXd deltaThisLayer;
      MatrixXd deltaAsDiagonal;
      //\delta_l = ((weights_{l+1, l}^T\delta^{l+1})) .* tanh'(z^l);
      //where ^T means transpose, l means lth layer, tanh' means the derivative
      //of tanh and z^l= tanh(weights_{l,l-1}activations^{l-1}) is the input signal 
      //on the lth layer neurons, where weights_{l, l-1} corresponds to the 
      //connections between the lth and l-1 th layers of neurons. Notice
      // .* means elementwise multiply.
      deltaThisLayer = weights[layer+1].transpose() * deltaPreviousLayer; 
      //To achive elementwise multiply, form a diagonal vector.
      deltaAsDiagonal = deltaThisLayer.asDiagonal();
      deltaAsDiagonal = deltaAsDiagonal 
        * inputSignal_Epochs[epoch][layer+1].unaryExpr(&Tanh_prime);
      deltaThisLayer = deltaAsDiagonal * 
        VectorXd::Ones(deltaThisLayer.size());
      nabla_biases[layer] += deltaThisLayer;
      //get a weight matrix. \partial C/\partial w^{l}_{jk} = a^{l-1}_k \delta_j^l 
      //the layer here refers to the lth layer of Biases and weights, so for
      //activation layer refers to the l-1th layer.
      nabla_weights[layer] += activations[layer] * deltaThisLayer.transpose() * ;
    }
  }
}

std::vector<double> NeuralNetwork::NablaE_C(){
  std::vector<double> dEdC;
  for (int i=0; i < list.size(); ++i){
    double dEdC_i(0.);
    for (int j=0; j < list.size(); ++j){
      dEdC_i += 2* output_Cs[i] * output_Cs[j] * H(list[i], list[j]);
    }
    dEdC.push_back(dEdC_i);
  }
  return dEdC;
}
