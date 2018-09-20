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

class ConvLayer: public Layer {
public:
  using T=Eigen::Matrix<paraType, Eigen::Dynamic, 1>;
  using weightType = std::vector<std::vector<Eigen::Map<Eigen::Matrix<paraType, Eigen::Dynamic, Eigen::Dynamic>>>>;
  using biasType = std::vector<Eigen::Map<Eigen::Matrix<paraType, Eigen::Dynamic, 1>>>;
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
};

}

#endif
