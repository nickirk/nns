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

const std::complex<double> ii(0.,1.);
class NeuralNetwork{
  public:
    NeuralNetwork(std::vector<int> const &sizes_, CostFunction const &externalCF);
    coeffType getCoeff(detType const &det) const;
    void train(std::vector<detType> const&listDetsToTrain, double eta, int iteration_);
    //double getEnergy(){return cf->calc(outputState);}
    //State getState() const {return outputState;}
    double calcEnergy(std::vector<detType> const &listDetsToTrain) const;
    Eigen::VectorXd feedForward(detType const& det) const;
    void updateParameters(int method, std::vector<State> const &outputState, double learningRate, int iteration);
    void setCostFunction(CostFunction const &externalCF) {cf = &externalCF;}
    CostFunction const* getCostFunction() const {return cf;}
    Eigen::Map<Eigen::MatrixXd> getWeights(int layer) const {return weights[layer];}
    Eigen::Map<Eigen::VectorXd> getBiases(int layer) const {return biases[layer];}
  private:
    double momentumDamping;
    bool momentum;
    double lamdaS1;
    double lamdaS;
    double gammaS;
    double gammaS1;
    Eigen::VectorXd yS;
    Eigen::VectorXd yS1;
    Eigen::VectorXd Egz2;
    std::vector<int> sizes;
   //----------------------
   //variables for ADAM
   double beta1;
   double beta2;
   Eigen::VectorXd m;
   Eigen::VectorXd v;
   Eigen::VectorXd m1;
   Eigen::VectorXd v1;
   //---------------------
    int numNNP;
    Eigen::VectorXd NNP;
    Eigen::VectorXd generlisedForcePrev;
    double *adNNP;
    Eigen::VectorXd nablaNNP;
    double *adNablaNNP;
    //These four are cache variables that are used to store the state of the network
    //They are potentially very memory expensive, we want to get rid of them first
    //mutable std::vector<std::vector<Eigen::VectorXd>> inputSignalsEpochs;
    //mutable std::vector<std::vector<Eigen::VectorXd>> activationsEpochs;
    mutable std::vector<Eigen::VectorXd> inputSignals; 
    mutable std::vector<Eigen::VectorXd> activations; 
    std::vector<Eigen::Map<Eigen::VectorXd>> biases;
    std::vector<Eigen::Map<Eigen::MatrixXd>> weights;
    std::vector<Eigen::Map<Eigen::VectorXd>> nablaBiases;
    std::vector<Eigen::Map<Eigen::MatrixXd>> nablaWeights;
    CostFunction const *cf;
    Solver sl;
    Eigen::VectorXd backPropagate(
      Eigen::VectorXd const &lastLayerFeedBack
         );
    Eigen::VectorXd calcNablaNNP(
      std::vector<State> const &outputState
         );
    Eigen::VectorXd calcNablaNNPMk(
	   std::vector<State> const &outputState
     );
    Eigen::MatrixXcd calcdCdwSR(
      std::vector<State> const &outputState
         );
    coeffType outputLayer() const {	size_t numLayers{sizes.size()};
    	return coeffType(activations[numLayers-1][0],activations[numLayers-1][1]);}
};

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
