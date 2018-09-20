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
  using weightType = std::vector<std::vector<Eigen::Map<Eigen::Matrix<paraType, Eigen::Dynamic, Eigen::Dynamic>>>>;
  using biasType = std::vector<Eigen::Map<Eigen::Matrix<paraType, Eigen::Dynamic, 1>>>;
  InputLayer(std::vector<Eigen::VectorXd> const &inputs_, int size_);
  virtual ~InputLayer();
  // implement the clone() functionality
  virtual InputLayer* clone() const {return new InputLayer(*this);}
private:
  // back-propagating on an input layer does not do anything
  virtual void backProp(std::vector<paraVector> const &prevDelta,
                        weightType const &prevWeights){};
  // same for mapping parameters - nothing will happen
  virtual void mapPara(double *adNNP, double *adNablaNNP, int &startPoint){};
  // input layers don't do anything when processing an empty signal
  void processSignal()const{};
  // number of neurons (better be the number of orbitals
  int numNrn;
};

}
#endif
