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
using namespace Eigen;
NeuralNetwork::NeuralNetwork(std::vector<int> sizes_){
  sizes = sizes_;
  for (int i=0; i< sizes.size(); ++i){
    activations.push_back(VectorXd::Zero(sizes[i]));
    if (i==0 || i==sizes.size()-1) biases.push_back(VectorXd::Zero(sizes[0]));
    biases.push_back(VectorXd::Random(sizes[i]));
    if (i < sizes.size()-1) weights.push_back(MatrixXd::Random(sizes[i], sizes[i+1]));
  }
}

double calcEnergy();
