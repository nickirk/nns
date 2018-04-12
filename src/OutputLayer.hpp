/*
 * OutputLayer.hpp
 * Created on 10.04.2018
 * Author: Ke Liao 
 */

#ifndef OutputLayer_DEFINED
#define OutputLayer_DEFINED

#include <vector>
#include <Eigen/Dense>
#include "DenseLayer.hpp"

class OutputLayer: public DenseLayer{
public:
  OutputLayer(std::vector<Eigen::VectorXd> const &inputs_, 
             double &(actFunc_)(double), int size_) 
  virtual void backProp(Eigen::VectorXd dCostdC);
}

#endif
