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
  DenseLayer(std::vector<Eigen::MatrixXd> const &inputs_, 
             double &(actFunc_)(double), int size_) 
  virtual ~DenseLayer();
  void processSignal();
  void mapPara(double *adNNP, int &startPoint);
  void backProp();
private:
  int numNrn;
  //z=w*input+bias, which is an intermediate variable but is needed
  //during the backpropagation step
  std::vector<Eigen::VectorXd> z;
  std::vector<Eigen::Map<Eigen::VectorXd>> biases;
  std::vector<Eigen::Map<Eigen::MatrixXd>> weights;
}

#endif
