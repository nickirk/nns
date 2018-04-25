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
  virtual void processSignal() const;
  virtual void backProp(
                        std::vector<Eigen::VectorXd> const &prevDelta,
                        weightType const &prevWeights
                        );
  virtual void backProp(
                        Eigen::VectorXd const &prevDelta
                        );
  virtual void mapPara(double *adNNP, double *adNablaNNP, int &startPoint);
  virtual int getNumPara() const{return numPara;};
  //want to achieve the effect of returning reference instead of making a 
  //copy. But Eigen::Map is already a referece object..
  //const weightType & getWeights(){return weights;};
  //const biasType & getBiases(){return biases;};
 
protected:
  int numNrn;
  //biasType biases;
  //weightType weights;
  //biasType nablaBiases;
  //weightType nablaWeights;
};

#endif
