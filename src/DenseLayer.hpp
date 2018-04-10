/*
 * DenseLayer.hpp
 * Created on 06.3.2018
 * Author: Ke Liao 
 */

#ifndef DenseLayer_DEFINED
#define DenseLayer_DEFINED

#include <vector>
#include <Eigen/Dense>
#include "Layer.hpp"

class DenseLayer: public Layer{
public:
  DenseLayer(std::vector<Eigen::VectorXd> const &inputs_, 
             double &(actFunc_)(double), int size_) 
  virtual ~DenseLayer();
  virtual void processSignal();
  virtual void backProp();
  virtual void mapPara(double *adNNP, int &startPoint);
  virtual int getNumPara(return numPara;);
  //want to achieve the effect of returning reference instead of making a 
  //copy. But Eigen::Map is already a referece object..
  weightType & getWeights(return &weights) const;
  biasType & getbiases(return &biases) const;
 
private:
  int numPara;
  int numNrn;
  //z=w*input+bias, which is an intermediate variable but is needed
  //during the backpropagation step
  std::vector<Eigen::VectorXd> z;
  biasType biases;
  weightType weights;
  biasType nablaBiases;
  weightType nablaWeights;
}

#endif
