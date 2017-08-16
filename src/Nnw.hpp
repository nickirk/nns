/*
 * Nnw.hpp
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */

#ifndef NeuralNetwork_DEFINED
#define NeuralNetwork_DEFINED

#include <vector>
#include <stdio.h>
#include <math.h>
#include <Eigen/Dense>
#include "Determinant.hpp"
#include "Basis.hpp"
#include "Hamiltonian.hpp"
#include "Sampler.hpp"
using namespace Eigen;
class NeuralNetwork{
  public:
    NeuralNetwork(std::vector<int> sizes_, Hamiltonian const&H_,
      Basis const &fullBasis_);
    void train();

  private:
    std::vector<int> sizes;
    std::vector<VectorXd> activations;
    std::vector<VectorXd> biases;
    std::vector<MatrixXd> weights;
    std::vector<double> output_Cs;
    std::vector<VectorXd> nabla_biases;
    std::vector<MatrixXd> nabla_weights;
    Hamiltonian &H;
    Basis &fullBasis;
    double energy;
    void calcEnergy();
    void train(std::vector<detType> &list, double eta);
    //double feedForward(double activation_);
    void backPropagate(int numDets);
    std::vector<double> NablaE_C();
};

double Tanh_prime(double in){return 1-tanh(in)*tanh(in);}
double Tanh(double in) {return tanh(in);}
#endif
