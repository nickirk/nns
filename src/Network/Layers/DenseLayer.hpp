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

template <typename F=std::complex<double>, typename coeffType=std::complex<double>>
class DenseLayer: public Layer<F, coeffType>{
public:
  using T=Eigen::Matrix<F, Eigen::Dynamic, 1>;
  using weightType = std::vector<std::vector<Eigen::Map<Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic>>>>;
  using biasType = std::vector<Eigen::Map<Eigen::Matrix<F, Eigen::Dynamic, 1>>>;
  DenseLayer(std::vector<T> const &inputs_,
             std::string actFunc_, int size_);
  virtual ~DenseLayer();
  // implement the clone() functionality
  virtual DenseLayer* clone() const{return new DenseLayer(*this);}
  virtual void processSignal() const;
  virtual void backProp(
                        std::vector<T> const &prevDelta,
                        weightType const &prevWeights
                        );
  virtual void backProp(coeffType const &prevDelta_);
  virtual void mapPara(double *adNNP, double *adNablaNNP, int &startPoint);
  virtual int getNumPara() const{return numPara;};
  //want to achieve the effect of returning reference instead of making a 
  //copy. But Eigen::Map is already a referece object..
  //const weightType & getWeights(){return weights;};
  //const biasType & getBiases(){return biases;};
 
private:
  int numNrn;
  using Layer<F, coeffType>::numPara;
  using Layer<F, coeffType>::activations;
  using Layer<F, coeffType>::inputs;
  using Layer<F, coeffType>::deltas;
  using Layer<F, coeffType>::z;
  using Layer<F, coeffType>::weights;
  using Layer<F, coeffType>::nablaWeights;
  using Layer<F, coeffType>::biases;
  using Layer<F, coeffType>::nablaBiases;
  using Layer<F, coeffType>::actFunc;
  using Layer<F, coeffType>::actFuncPrime;
};

}

#endif
