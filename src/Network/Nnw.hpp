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
#include "../CostFunctions/CostFunction.hpp"
#include "../Hamiltonian/TwoBodyHamiltonian.hpp"
#include "../HilbertSpace/Basis.hpp"
#include "../HilbertSpace/Determinant.hpp"
#include "../utilities/State.hpp"
#include "../utilities/TypeDefine.hpp"
#include "../math/MathFunctions.hpp"
#include "Layers/ConvLayer.hpp"
#include "Layers/DenseLayer.hpp"
#include "Layers/InputLayer.hpp"
#include "Layers/Layer.hpp"
#include "../Solvers/Solver.hpp"
#include "LayerStructure.hpp"
#include "Parametrization.hpp"

namespace networkVMC{

const std::complex<double> ii(0.,1.);

// Neural network parametrization of a wavefunction
template<typename T=VecType>
class NeuralNetwork: public ClonableParametrization<T,NeuralNetwork<T> >{
public:
  NeuralNetwork();
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

  // implementation of updateParameters
  void updateParameters(State const &outputState,
                        double learningRate, int iteration, int method = 3);

  //implementation of getCoeff()
  coeffType getCoeff(detType const &det) const; // can throw an EmptyNetworkError if no layers exist

  //implementation of pars()
  VecType const& pars() const {return NNP;}

  //feed forward: iterative calculate the activations
  VecType feedForward(detType const& det) const; // can throw an EmptyNetworkError if no layers exist
  Layer* getLayer(int layer){
	  if(static_cast<unsigned int>(layer)<Layers.size() && layer >= 0)
	  return Layers[layer];
	  throw OutOfRangeError(layer);
  };
  //Eigen::Map<Eigen::MatrixXd> getWeights(int layer) const {return weights[layer];};
  //Eigen::Map<Eigen::VectorXd> getBiases(int layer) const {return biases[layer];};

  VecType calcNablaPars(
    State const &inputState,
	nablaType const& dEdC
  );
  // derivative taking into account connected dets
  VecType calcNablaParsConnected(State const &inputState,nablaType const&dEdC);
  // derivative of higher order(?)
  //Eigen::MatrixXcd calcdCdwSR(State const &outputState);

private:
  //Structure of the NNW
  LayerStructure Layers;
  //variables for RMSprop
  std::vector<int> sizes;
 //---------------------
 //variable for the network
  int numNNP;
  int numLayers;
  Eigen::VectorXd NNP;
  Eigen::VectorXd generlisedForcePrev;
  Eigen::VectorXd nablaNNP;
  std::vector<Eigen::VectorXd> feedIns;

  VecType backPropagate(
    coeffType const &lastLayerFeedBack
       );
  void copyNetwork(NeuralNetwork const &source);
  //coeffType outputLayer() const {
  //	return coeffType(Layers[numLayers-1]->getActs()[0][0],Layers[numLayers-1]->getActs()[0][1]);
  //};
};

}
#endif
