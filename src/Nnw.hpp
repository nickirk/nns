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
#include <complex>
#include <string>
#include <Eigen/Dense>
#include "math/MathFunctions.hpp"
#include "CostFunction.hpp"
#include "Determinant.hpp"
#include "Basis.hpp"
#include "Hamiltonian.hpp"
#include "State.hpp"
#include "Solver.hpp"
#include "Layer.hpp"
#include "InputLayer.hpp"
#include "DenseLayer.hpp"
#include "ConvLayer.hpp"
#include "TypeDefine.hpp"

const std::complex<double> ii(0.,1.);
class NeuralNetwork{
public:
  NeuralNetwork();
  ~NeuralNetwork();
  NeuralNetwork(Hamiltonian const &H_,
                CostFunction const &externalCF);
  //construction functions for NNW
  void constrDenseLayer(
      std::vector<Eigen::VectorXd> const &inputs_,
      std::string actFunc_,
      int size_);

  void constrConvLayer(
    std::vector<Eigen::VectorXd> const &inputs_,
    std::string actFunc_,
    int numFilters_,
    int lengthFilter_,
    int stride_
    );

  void constrInputLayer(int numNrn);
 /* 
  void constrOutputLayer(
      std::vector<Eigen::VectorXd> const &inputs_, 
      double &(actFunc_)(double), 
      int size_
      );
 */
  void initialiseNetwork();

  void setCostFunction(CostFunction const &externalCF) {cf = &externalCF;};
  //functionalities of NNW
  Eigen::VectorXd feedForward(detType const& det) const;
  void updateParameters(int method, std::vector<State> const &outputState, 
                        double learningRate, int iteration);

  //interface API
  coeffType getCoeff(detType const &det) const;
  CostFunction const* getCostFunction() const {return cf;};
  Layer* getLayer(int layer){return Layers[layer];};
  //Eigen::Map<Eigen::MatrixXd> getWeights(int layer) const {return weights[layer];};
  //Eigen::Map<Eigen::VectorXd> getBiases(int layer) const {return biases[layer];};

private:
  // Hamiltonian
  Hamiltonian const &H;
  //Structure of the NNW
  std::vector<Layer*> Layers;
  //variables for RMSprop
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
 //variable for the network
  int numNNP;
  int numLayers;
  Eigen::VectorXd NNP;
  Eigen::VectorXd generlisedForcePrev;
  double *adNNP;
  Eigen::VectorXd nablaNNP;
  double *adNablaNNP;
  std::vector<Eigen::VectorXd> feedIns;

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
  //coeffType outputLayer() const {
  //	return coeffType(Layers[numLayers-1]->getActs()[0][0],Layers[numLayers-1]->getActs()[0][1]);
  //};
};
#endif
