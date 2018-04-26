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
#include "../HilbertSpace/Basis.hpp"
#include "../HilbertSpace/Determinant.hpp"
#include "../CostFunctions/CostFunction.hpp"
#include "../utilities/State.hpp"
#include "../utilities/TypeDefine.hpp"
#include "../math/MathFunctions.hpp"
#include "Layers/ConvLayer.hpp"
#include "Layers/DenseLayer.hpp"
#include "Layers/InputLayer.hpp"
#include "Layers/Layer.hpp"
#include "../Solver.hpp"
#include "../Hamiltonian/Hamiltonian.hpp"
#include "LayerStructure.hpp"

namespace networkVMC{

const std::complex<double> ii(0.,1.);
class NeuralNetwork{
public:
  NeuralNetwork();
  ~NeuralNetwork();
  NeuralNetwork(NeuralNetwork const &source);
  NeuralNetwork(CostFunction const &externalCF);
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
  Layer* getLayer(int layer){
	  if(layer<Layers.size() && layer >= 0)
	  return Layers[layer];
	  throw OutOfRangeError(layer);
  };
  //Eigen::Map<Eigen::MatrixXd> getWeights(int layer) const {return weights[layer];};
  //Eigen::Map<Eigen::VectorXd> getBiases(int layer) const {return biases[layer];};

private:
  //Structure of the NNW
  LayerStructure Layers;
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
  Eigen::VectorXd nablaNNP;
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
  void copyNetwork(NeuralNetwork const &source);
  //coeffType outputLayer() const {
  //	return coeffType(Layers[numLayers-1]->getActs()[0][0],Layers[numLayers-1]->getActs()[0][1]);
  //};
};

}
#endif
