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

// First layer of a nerual network
class InputLayer: public Layer{
public:
  InputLayer(std::vector<Eigen::VectorXd> const &inputs_, int size_);
  virtual ~InputLayer();
  // implement the clone() functionality
  virtual InputLayer* clone() const {return new InputLayer(*this);}
  void processSignal(detType const det) const;
private:
  int numNrn;
};

#endif
