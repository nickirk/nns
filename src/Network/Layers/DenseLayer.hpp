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

#include "../../math/MathFunctions.hpp"
#include "Layer.hpp"

namespace networkVMC{

class DenseLayer: public Layer{
public:
  DenseLayer(std::vector<Eigen::VectorXd> const &inputs_, 
             std::string actFunc_, int size_);
  virtual ~DenseLayer();
  // implement the clone() functionality
  virtual DenseLayer* clone() const{return new DenseLayer(*this);}
  virtual void processSignal() const;
  virtual void backProp(
                        std::vector<Eigen::VectorXd> const &prevDelta,
                        weightType const &prevWeights
                        );
  virtual void backProp(coeffType const &prevDelta);
  virtual void mapPara(double *adNNP, double *adNablaNNP, int &startPoint);
  virtual int getNumPara() const{return numPara;};
  //want to achieve the effect of returning reference instead of making a 
  //copy. But Eigen::Map is already a referece object..
  //const weightType & getWeights(){return weights;};
  //const biasType & getBiases(){return biases;};
 
private:
  int numNrn;
};

}

#endif