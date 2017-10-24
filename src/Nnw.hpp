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
#include <complex>
#include <Eigen/Dense>
#include "Determinant.hpp"
#include "Basis.hpp"
#include "Hamiltonian.hpp"
//#include "Sampler.hpp"
const std::complex<double> ii(0.,1.);
class NeuralNetwork{
  public:
    NeuralNetwork(std::vector<int> const &sizes_, Hamiltonian const&H_,
      Basis const &fullBasis_);
    std::vector<detType> train(std::vector<detType> const&listDetsToTrain, double eta);
    double getEnergy(){return energy;}
    std::vector<Eigen::VectorXd> getCs() const {return outputCs;}
    double calcEnergy(std::vector<detType> const &listDetsToTrain) const;
    Eigen::VectorXd feedForward(detType const& det);
  private:
    double momentumDamping;
    bool momentum;
    std::vector<int> sizes;
    std::vector<Eigen::VectorXd> biases;
    std::vector<Eigen::MatrixXd> weights;
    std::vector<Eigen::VectorXd> outputCs;
    std::vector<Eigen::VectorXd> inputSignal; 
    std::vector<Eigen::VectorXd> activations; 
    std::vector<Eigen::VectorXd> nablaBiases;
    std::vector<Eigen::MatrixXd> nablaWeights;
    std::vector<Eigen::VectorXd> nablaBiasesPrev;
    std::vector<Eigen::MatrixXd> nablaWeightsPrev;
    std::vector<Eigen::VectorXd> gFactorBiases;
    std::vector<Eigen::MatrixXd> gFactorWeights;
    std::vector<Eigen::VectorXd> gFactorBiasesPrev;
    std::vector<Eigen::MatrixXd> gFactorWeightsPrev;
    //std::vector<detType> &listDetsToTrain;
    Hamiltonian const&H;
    Basis const&fullBasis;
    double energy;
    mutable double normalizerCoeff;
    //void calcLocalEnergy(std::vector<detType> const &listDetsToTrain);
    //std::vector<detType> train(std::vector<detType> &listDetsToTrain_, double eta);
    //double feedForward(double activation_);
    void backPropagate(
           std::vector<detType> const &listDetsToTrain,
           std::vector<std::vector<Eigen::VectorXd>> const &inputSignalEpochs,
           std::vector<std::vector<Eigen::VectorXd>> const &activations
         );
    std::vector<Eigen::VectorXd> NablaE_C(std::vector<detType> const &listDetsToTrain);
};

double Tanh_prime(double in);
double Tanh(double in); 
double Linear(double in);
double Linear_prime(double in);
double Gaussian(double in);
double Gaussian_prime(double in);
double Sigmoid(double in);
double Sigmoid_prime(double in);
#endif
