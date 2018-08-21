/*
 * ConvLayer.hpp
 * Created on 27.3.2018
 * Author: Ke Liao 
 */
#ifndef ConvLayer_DEFINED
#define ConvLayer_DEFINED

#include <vector>
#include <Eigen/Dense>
#include "../../math/MathFunctions.hpp"
#include "Layer.hpp"

namespace networkVMC{

template <typename F=std::complex<double>, typename coeffType=std::complex<double>>
class ConvLayer: public Layer<F, coeffType> {
public:
  using T=Eigen::Matrix<F, Eigen::Dynamic, 1>;
  using weightType = std::vector<std::vector<Eigen::Map<Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic>>>>;
  using biasType = std::vector<Eigen::Map<Eigen::Matrix<F, Eigen::Dynamic, 1>>>;
  ConvLayer(
          std::vector<T> const &inputs_,
          std::string actFunc_, 
          int numFilters_, 
          int lengthFilter_, 
          int stride_
            ); 
  virtual ~ConvLayer();
  // implement the clone() functionality
  virtual ConvLayer* clone() const{return new ConvLayer(*this);}
  virtual void backProp(
    std::vector<T> const &prevDelta,
    weightType const &prevWeights
    );
  virtual void processSignal() const;
  virtual void mapPara(double *adNNP, double *adNablaNNP, int &startPoint);
  int getStride() const{return stride;};
  //interface
private:
  int numFilters;
  int depthFilter;
  int lengthFilter;
  int stride;
  //sizeAct is the size of the activations of each filter.
  int sizeAct;
  // z=w*input+bias, which is an intermediate variable but is needed
  // during the backpropagation step
  //std::vector<T> z;
  // each filter (consists of several chanels) shares the same 
  // bias. So the dimension will be L vector of double
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
