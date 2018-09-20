/*
 * Layer.hpp
 * Created on 06.3.2018
 * Author: Ke Liao 
 */

#ifndef Layer_DEFINED
#define Layer_DEFINED

#include <vector>
#include <Eigen/Dense>
#include <string>

#include "../../utilities/TypeDefine.hpp"
#include "../../utilities/MatrixTypeDefine.hpp"
#include "../../math/MathFunctions.hpp"
#include "../../utilities/Errors.hpp"

namespace networkVMC{

// single layer of a neural network (abstract base class)
class Layer{
public:
  using weightType = std::vector<std::vector<Eigen::Map<Eigen::Matrix<paraType, Eigen::Dynamic, Eigen::Dynamic>>>>;
  using biasType = std::vector<Eigen::Map<Eigen::Matrix<paraType, Eigen::Dynamic, 1>>>;
  Layer(std::vector<Eigen::VectorXd> const &inputs_, std::string actFunc_);
  virtual ~Layer();
  // virtual function for copying polymorphic layers
  virtual Layer* clone() const=0;
  //functional functions
  virtual void processSignal() const=0;
  // just for convenience, we want to be able to call process signal with an argument
  // (to prevent unfortunate bugs with the inputlayer, which overrides this method)
  virtual void processSignal(detType const &det) const;
  // update of parameters will be done on network level, not on layer level.
  virtual void backProp(
                        std::vector<paraVector> const &prevDelta,
                        weightType const &prevWeights
                        )=0;
  virtual void backProp(coeffType const &prevDelta_){};
  virtual void mapPara(paraType *adNNP, paraType *adNablaNNP, int &startPoint)=0;
  //virtual void initialise(T const &NNP);
  //interface for accessing data
  virtual const std::vector<paraVector>& getInputs() const{return inputs;};
  virtual const std::vector<paraVector>& getActs() const{return activations;};
  //getWeights in base function has no return value and it is not defined for
  //inputLayer, this may cause problem. Default return a empty vector.
  //or move weights to base class
  const weightType & getWeights() const{return weights;};
  const biasType & getBiases() const{return biases;};
  virtual int getNumPara() const{return numPara;};
  const std::vector<paraVector>& getDeltas(){return deltas;}
protected:
  double (*actFunc)(double);
  double (*actFuncPrime)(double);
  int numPara;
  std::vector<paraVector> const &inputs;
  // The activations are auxiliary internal variables and not states of the layer
  // (same for z, they are just functions of the state)
  mutable std::vector<paraVector> activations;
  //z=w*input+bias, which is an intermediate variable but is needed
  //during the backpropagation step
  mutable std::vector<paraVector> z;
  std::vector<paraVector> deltas;
  biasType biases;
  // the weights have one-to-one correspondance to neurons. 
  // in convNet, each filter has the same depth as that of the 
  // previous output signal. So we define it as a Lxlxsize
  weightType weights;
  biasType nablaBiases;
  weightType nablaWeights;
};

}

#endif
