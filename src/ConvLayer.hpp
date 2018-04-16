/*
 * ConvLayer.hpp
 * Created on 27.3.2018
 * Author: Ke Liao 
 */
#ifndef ConvLayer_DEFINED
#define ConvLayer_DEFINED

#include <vector>
#include <Eigen/Dense>
#include "Layer.hpp"

class ConvLayer: public Layer {
public:
  ConvLayer(
          std::vector<Eigen::VectorXd> const &inputs_, 
          double &(actFunc_)(double), 
          int lengthFilter_, int numFilters_, int stride_
            ) 
  virtual ~ConvLayer();
  virtual void backProp(
    std::vector<Eigen::VectorXd> const &prevDelta,
    weightType const &prevWeights
    );
  virtual void mapPara(double *adNNP, int &startPoint);
  virtual void backProp();
  virtual int getNumPara(){return numPara;};
  //interface
  weightType & getWeights(return &weights) const;
  biasType & getbiases(return &biases) const;
private:
  int numPara;
  int numFilters;
  int depthFilter;
  int lengthFilter;
  int stride;
  //sizeAct is the size of the activations of each filter.
  int sizeAct;
  // z=w*input+bias, which is an intermediate variable but is needed
  // during the backpropagation step
  std::vector<Eigen::VectorXd> z;
  // each filter (consists of several chanels) shares the same 
  // bias. So the dimension will be L vector of double
  biasType biases;
  // the weights have one-to-one correspondance to neurons. 
  // in convNet, each filter has the same depth as that of the 
  // previous output signal. So we define it as a Lxlxsize
  weightType weights;
}


#endif
