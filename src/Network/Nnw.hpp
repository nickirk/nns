/*
 * Nnw.hpp
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */

#ifndef NeuralNetwork_DEFINED
#define NeuralNetwork_DEFINED

#include <vector>
#include <complex>
#include <Eigen/Dense>
#include "../HilbertSpace/Basis.hpp"
#include "../utilities/StateForward.hpp"
#include "../utilities/TypeDefine.hpp"
#include "Layers/Layer.hpp"
#include "LayerStructure.hpp"
#include "Parametrization.hpp"

namespace networkVMC{

const std::complex<double> ii(0.,1.);

// Neural network parametrization of a wavefunction
template <typename F=std::complex<double>, typename coeffType=std::complex<double>>
class NeuralNetwork: public ClonableParametrization<F, coeffType, NeuralNetwork<F, coeffType> >{
  public:
    using T=Eigen::Matrix<F, Eigen::Dynamic, 1>;
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
  void updateParameters(State<coeffType> const &outputState,
                        double learningRate, int iteration, int method = 3);

  //implementation of getCoeff()
  coeffType getCoeff(detType const &det) const; // can throw an EmptyNetworkError if no layers exist

  //implementation of pars()
  T const& pars() const {return NNP;}

  //feed forward: iterative calculate the activations
  T feedForward(detType const& det) const; // can throw an EmptyNetworkError if no layers exist
  Layer<F, coeffType>* getLayer(int layer){
	  if(static_cast<unsigned int>(layer)<Layers.size() && layer >= 0)
	  return Layers[layer];
	  throw errors::OutOfRangeError(layer);
  };
  //Eigen::Map<Eigen::MatrixXd> getWeights(int layer) const {return weights[layer];};
  //Eigen::Map<Eigen::VectorXd> getBiases(int layer) const {return biases[layer];};

  T calcNablaPars(
    State<coeffType> const &inputState,
	T const& dEdC
  );
  // derivative taking into account connected dets
  T calcNablaParsConnected(State<coeffType> const &inputState,T const&dEdC);
  // derivative of higher order(?)
  //Eigen::MatrixXcd calcdCdwSR(State const &outputState);

private:
  //Structure of the NNW
  LayerStructure<F, coeffType> Layers;
  //variables for RMSprop
  std::vector<int> sizes;
 //---------------------
 //variable for the network
  int numLayers;
  int numNNP;
  T NNP;
  T generlisedForcePrev;
  T nablaNNP;
  std::vector<Eigen::VectorXd> feedIns;

  T backPropagate(
    F const &lastLayerFeedBack
       );
  void copyNetwork(NeuralNetwork<F, coeffType> const &source);
  //F outputLayer() const {
  //	return F(Layers[numLayers-1]->getActs()[0][0],Layers[numLayers-1]->getActs()[0][1]);
  //};
};

}
#endif
