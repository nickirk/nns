/*
 * DenseLayer.hpp
 * Created on 06.3.2018
 * Author: Ke Liao 
 */

#ifndef DenseLayer_DEFINED
#define DenseLayer_DEFINED

#include <vector>
#include <Eigen/Dense>
#include <string>
#include "Layer.hpp"
#include "math/MathFunctions.hpp"

class DenseLayer: public Layer{
public:
  DenseLayer(std::vector<Eigen::VectorXd> const &inputs_, 
             std::string actFunc_, int size_);
  virtual ~DenseLayer();
  virtual void processSignal();
  virtual void backProp(
                        std::vector<Eigen::VectorXd> const &prevDelta,
                        weightType const &prevWeights
                        );
  virtual void backProp(
                        Eigen::VectorXd dCostdC
                        );
  virtual void mapPara(double *adNNP, double *adNablaNNP, int &startPoint);
  virtual int getNumPara(){return numPara;};
  //want to achieve the effect of returning reference instead of making a 
  //copy. But Eigen::Map is already a referece object..
  virtual const weightType & getWeights(){return weights;};
  const biasType & getbiases(){return biases;};
 
protected:
  int numPara;
  int numNrn;
  biasType biases;
  weightType weights;
  biasType nablaBiases;
  weightType nablaWeights;
};

#endif
