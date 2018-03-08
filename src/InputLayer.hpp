/*
 * InputLayer.hpp
 * Created on 08.3.2018
 * Author: Ke Liao 
 */

#ifndef InputLayer_DEFINED
#define InputLayer_DEFINED

#include <vector>
#include <Eigen/Dense>
#include "Layer.hpp"

class InputLayer: public Layer{
public:
  InputLayer(std::vector<Eigen::VectorXd> const &inputs_, int size_) 
  virtual ~InputLayer();
  void processSignal();
private:
  int numNrn;
}

#endif
