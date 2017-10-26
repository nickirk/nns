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
#include "CostFunction.hpp"
#include "Determinant.hpp"
#include "Basis.hpp"
#include "Hamiltonian.hpp"
#include "Evaluator.hpp"
#include "State.hpp"
//#include "Sampler.hpp"
const std::complex<double> ii(0.,1.);
class NeuralNetwork{
  public:
    NeuralNetwork(std::vector<int> const &sizes_, CostFunction const &externalCF);
    std::vector<detType> train(std::vector<detType> const&listDetsToTrain, double eta);
    double getEnergy(){return energy;}
    State getState() const {return outputState;}
    double calcEnergy(std::vector<detType> const &listDetsToTrain) const;
    Eigen::VectorXd feedForward(detType const& det);
    void setCostFunction(CostFunction const &externalCF) {cf = &externalCF;}
  private:
    double momentumDamping;
    bool momentum;
    State outputState;
    std::vector<int> sizes;
    std::vector<Eigen::VectorXd> biases;
    std::vector<Eigen::MatrixXd> weights;
    std::vector<Eigen::VectorXd> inputSignal; 
    std::vector<Eigen::VectorXd> activations; 
    std::vector<Eigen::VectorXd> nablaBiases;
    std::vector<Eigen::MatrixXd> nablaWeights;
    std::vector<Eigen::VectorXd> nablaBiasesPrev;
    std::vector<Eigen::MatrixXd> nablaWeightsPrev;
    std::vector<Eigen::VectorXd> gBiasesPrev;
    std::vector<Eigen::MatrixXd> gWeightsPrev;
    std::vector<Eigen::VectorXd> gFactorBiases;
    std::vector<Eigen::MatrixXd> gFactorWeights;
    std::vector<Eigen::VectorXd> gFactorBiasesPrev;
    std::vector<Eigen::MatrixXd> gFactorWeightsPrev;
    CostFunction const *cf;
    double energy;
    void backPropagate(
           std::vector<std::vector<Eigen::VectorXd>> const &inputSignal_Epochs,
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
