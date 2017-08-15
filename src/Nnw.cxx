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
  for (int i=0; i< sizes.size(); ++i){
    activations.push_back(VectorXd::Zero(sizes[i]));
    if (i==0) {
      biases.push_back(VectorXd::Zero(sizes[0]));
      delta_biases.push_back(VectorXd::Zero(sizes[0]));
    }
    else{ 
      biases.push_back(VectorXd::Random(sizes[i]));
      delta_biases.push_back(VectorXd::Zero(sizes[i]));
    }
    if (i < sizes.size()-1) {
      weights.push_back(MatrixXd::Random(sizes[i+1], sizes[i]));
      delta_weights.push_back(MatrixXd::Zero(sizes[i+1], sizes[i]));
    }
  }
}

double NeuralNetwork::calcEnergy(std::vector<detType> &list){
  energy = 0.;
  for (int i=0; i < list.size(); ++i){
    for (int j=i; j < list.size(); ++j){
      Hij = 2 * H(list[i], list[j]);
      energy += output_Cs[i] * output_Cs[j] * Hij;
    }
  }
}

void() NeuralNetwork::train(std::vector<detType> &list, double eta){
  int num_Dets = list.size();
  std::vector<vector<VectorXd>> inputSignal_Epochs;
  std::vector<VectorXd> inputSignal; 
  for (int epoch=0; epoch < num_Dets; ++epoch){
    //initial input layer
    Map<VectorxXd> activations[0](*list[epoch][0], list[epoch][0].size());
    for (int layer=1; layer < sizes.size(); ++layer){
      activations[layer] = weights[layer-1]*activations[layer-1]+biases[layer];
      inputSignal.push_back(activations[layer]);
      activations[layer] = activations[layer].unaryExpr(&Tanh);
    }
    inputSignal_Epochs.push_back(input_Tanh);
    output_Cs.push_back(activations[sizes.size()-1]);
  }
  
}

//void NeuralNetwork::feedForward();
void backPropagate(int numDets){
  int numLayers(sizes.size());
  std::vector<double> dEdC(NablaE_C());
  std::vector<vector<double>> delta_Epoch;
  for (int epoch=0; epoch < numDets; ++i){
    std::vector<VectorXd> delta_Layer;
    delta_Layer.push_back(
      dEdC[epoch] * inputSignal_Epochs[epoch][numLayers-1].unaryExpr(&Tanh_prime)
    );
    for (int layer=numLayers-1; layer > 0; --layer){
      
    }
  }
  delta = 
    
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
