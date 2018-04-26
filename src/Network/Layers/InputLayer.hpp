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

namespace networkVMC{

// First layer of a nerual network
class InputLayer: public Layer{
public:
  InputLayer(std::vector<Eigen::VectorXd> const &inputs_, int size_);
  virtual ~InputLayer();
  // implement the clone() functionality
  virtual InputLayer* clone() const {return new InputLayer(*this);}
  // back-propagating on an input layer does not do anything
  virtual void backProp(std::vector<Eigen::VectorXd> const &prevDelta,
                        weightType const &prevWeights){};
  // same for mapping parameters - nothing will happen
  virtual void mapPara(double *adNNP, double *adNablaNNP, int &startPoint){};
  void processSignal(detType const &det) const;
  // input layers don't do anything when processing an empty signal
  void processSignal()const{};
private:
  int numNrn;
};

}
#endif
