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
//#include "Sampler.hpp"
class NeuralNetwork{
  public:
    NeuralNetwork(std::vector<int> const &sizes_, Hamiltonian const&H_,
      Basis const &fullBasis_);
    std::vector<detType> train(std::vector<detType> const&listDetsToTrain, double eta);
    double getEnergy(){return energy;}
    std::vector<double> getEnergyDerivative(std::vector<detType> const &list){return NablaE_C(list);}
    int getSign(){return sign;}
    std::vector<double> getCs() const {return output_Cs;}
    Eigen::VectorXd feedForward(detType const& det);
  private:
    int sign;
    std::vector<int> sizes;
    std::vector<Eigen::VectorXd> biases;
    std::vector<Eigen::MatrixXd> weights;
    std::vector<double> output_Cs;
    std::vector<Eigen::VectorXd> inputSignal; 
    std::vector<Eigen::VectorXd> activations; 
    std::vector<Eigen::VectorXd> nabla_biases;
    std::vector<Eigen::MatrixXd> nabla_weights;
    //std::vector<detType> &listDetsToTrain;
    Hamiltonian const&H;
    Basis const&fullBasis;
    double energy;
    void calcEnergy(std::vector<detType> const &listDetsToTrain);
    //std::vector<detType> train(std::vector<detType> &listDetsToTrain_, double eta);
    //double feedForward(double activation_);
    void backPropagate(
           std::vector<detType> const &listDetsToTrain,
           std::vector<std::vector<Eigen::VectorXd>> const &inputSignal_Epochs,
           std::vector<std::vector<Eigen::VectorXd>> const &activations
         );
    std::vector<double> NablaE_C(std::vector<detType> const &listDetsToTrain);
};

double Tanh_prime(double in);
double Tanh(double in); 
double Linear(double in);
double Linear_prime(double in);
#endif
