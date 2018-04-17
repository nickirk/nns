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
#include "TypeDefine.hpp"
#include "math/MathFunctions.hpp"
#include "utilities/Errors.hpp"

class Layer{
public:
  Layer(std::vector<Eigen::VectorXd> const &inputs_, std::string actFunc_);
  virtual ~Layer();
  //functional functions
  virtual void processSignal(){};
  virtual void processSignal(detType const det){};
  // update of parameters will be done on network level, not on layer level.
  //virtual void updatePara();
  virtual void backProp(
                        std::vector<Eigen::VectorXd> const &prevDelta, 
                        weightType const &prevWeights
                        ){};
  virtual void backProp(
                        Eigen::VectorXd const &prevDelta
                        ){};
  //virtual void mapPara(double *adNNP, int &startPoint);
  virtual void mapPara(double *adNNP, double *adNablaNNP, int &startPoint){};
  //virtual void initialise(Eigen::VectorXd const &NNP);
  //interface for accessing data
  virtual const std::vector<Eigen::VectorXd>& getInputs(){return inputs;};
  virtual const std::vector<Eigen::VectorXd>& getActs(){return activations;};
  //getWeights in base function has no return value and it is not defined for
  //inputLayer, this may cause problem. Default return a empty vector.
  //or move weights to base class
  //virtual const weightType & getWeights(){};
  //virtual const biasType & getBiases(){};
  const weightType & getWeights(){return weights;};
  const biasType & getbiases(){return biases;};
  virtual int getNumPara(){return numPara;};
  const std::vector<Eigen::VectorXd>& getDeltas(){return deltas;}
protected:
  double (*actFunc)(double);
  double (*actFuncPrime)(double);
  int numPara;
  std::vector<Eigen::VectorXd> const &inputs;
  std::vector<Eigen::VectorXd> activations;
  //z=w*input+bias, which is an intermediate variable but is needed
  //during the backpropagation step
  std::vector<Eigen::VectorXd> z;
  std::vector<Eigen::VectorXd> deltas;
  biasType biases;
  // the weights have one-to-one correspondance to neurons. 
  // in convNet, each filter has the same depth as that of the 
  // previous output signal. So we define it as a Lxlxsize
  weightType weights;
  biasType nablaBiases;
  weightType nablaWeights;
};

#endif
