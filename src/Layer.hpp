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
#include "math/MathFunctions.hpp"

class Layer{
public:
	Layer(std::vector<Eigen::VectorXd> const &inputs_, string actFunc_);
	virtual ~Layer();
  virtual void processSignal();
  virtual void updatePara();
  virtual std::vector<Eigen::VectorXd> getInputs(){return inputs;};
  virtual std::vector<Eigen::VectorXd> getActs(){return activations;};
  virtual void initialise(Eigen::VectorXd const &NNP);
  doulbe &actFunc; 
  doulbe &actFuncPrime; 
protected:
	std::vector<Eigen::VectorXd> inputs;
	std::vector<Eigen::VectorXd> activations;
}

#endif
