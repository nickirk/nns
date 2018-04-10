/*
 * Layer.hpp
 * Created on 06.3.2018
 * Author: Ke Liao 
 */

#ifndef Layer_DEFINED
#define Layer_DEFINED

#include <vector>
#include <Eigen/Dense>
#include <string>
#include "TypeDefine.hpp"
#include "math/MathFunctions.hpp"

class Layer{
public:
	Layer(std::vector<Eigen::VectorXd> const &inputs_, string actFunc_);
	virtual ~Layer();
  //functional functions
  virtual void processSignal();
  // update of parameters will be done on network level, not on layer level.
  //virtual void updatePara();
  virtual void backProp();
  virtual void mapPara(double *adNNP, int &startPoint);
  virtual void initialise(Eigen::VectorXd const &NNP);
  //interface for accessing data
  virtual std::vector<Eigen::VectorXd> getInputs(){return inputs;};
  virtual std::vector<Eigen::VectorXd> getActs(){return activations;};
  virtual int getNumPara();
protected:
  doulbe &actFunc; 
  doulbe &actFuncPrime; 
  int numPara;
	std::vector<Eigen::VectorXd> inputs;
	std::vector<Eigen::VectorXd> activations;
  std::vector<Eigen::VectorXd> deltas;
}

#endif
