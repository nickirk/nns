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
#include "State.hpp"
#include "Solver.hpp"
//#include "Sampler.hpp"
const std::complex<double> ii(0.,1.);
class NeuralNetwork{
  public:
    NeuralNetwork(std::vector<int> const &sizes_, CostFunction const &externalCF);
    void train(std::vector<detType> const&listDetsToTrain, double eta, int iteration);
    double getEnergy(){return cf->calc(outputState);}
    State getState() const {return outputState;}
    double calcEnergy(std::vector<detType> const &listDetsToTrain) const;
    Eigen::VectorXd feedForward(detType const& det) const;
    void setCostFunction(CostFunction const &externalCF) {cf = &externalCF;}
    CostFunction const* getCostFunction() const {return cf;}
  private:
    double momentumDamping;
    bool momentum;
    int iteration;
    State outputState;
    std::vector<int> sizes;
    int numNNP;
    Eigen::VectorXd NNP;
    double *adNNP;
    Eigen::VectorXd nablaNNP;
    double *adNablaNNP;
    std::vector<std::vector<Eigen::VectorXd>> inputSignalEpochs;
    std::vector<std::vector<Eigen::VectorXd>> activationsEpochs;
    mutable std::vector<Eigen::VectorXd> inputSignal; 
    mutable std::vector<Eigen::VectorXd> activations; 
    std::vector<Eigen::Map<Eigen::VectorXd>> biases;
    std::vector<Eigen::Map<Eigen::MatrixXd>> weights;
    std::vector<Eigen::Map<Eigen::VectorXd>> nablaBiases;
    std::vector<Eigen::Map<Eigen::MatrixXd>> nablaWeights;
    CostFunction const *cf;
    Solver sl;
    Eigen::VectorXd backPropagate(
           std::vector<std::vector<Eigen::VectorXd>> const &inputSignalEpochs,
           std::vector<std::vector<Eigen::VectorXd>> const &activationsEpochs
         );
    Eigen::MatrixXcd backPropagateSR(
           std::vector<std::vector<Eigen::VectorXd>> const &inputSignalEpochs,
           std::vector<std::vector<Eigen::VectorXd>> const &activationsEpochs
         );
    void updateParameters(int method);
    std::vector<Eigen::VectorXd> NablaE_C(std::vector<detType> const &listDetsToTrain);
};

void preTrain(NeuralNetwork &network, State const &target, int iteration);

double NormalDistribution(double dummy);
double Tanh_prime(double in);
double Tanh(double in); 
double Linear(double in);
double Linear_prime(double in);
double Rectifier(double in);
double Rectifier_prime(double in);
double Arcsinh(double in);
double Arcsinh_prime(double in);
double Gaussian(double in);
double Gaussian_prime(double in);
double GaussianAntiSym(double in);
double GaussianAntiSym_prime(double in);
double Sigmoid(double in);
double Sigmoid_prime(double in);
#endif
