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
  Layer(std::vector<Eigen::VectorXd> const &inputs_, int actFunc_);
  virtual ~Layer();
  //functional functions
  virtual void processSignal();
  // update of parameters will be done on network level, not on layer level.
  //virtual void updatePara();
  virtual void backProp(
                        std::vector<Eigen::VectorXd> prevDelta, 
                        weightType &prevWeights
                        );
  virtual void backProp(
                        Eigen::VectorXd prevDelta
                        );
  //virtual void mapPara(double *adNNP, int &startPoint);
  virtual void mapPara(double *adNNP, double *adNablaNNP, int &startPoint){};
  virtual void initialise(Eigen::VectorXd const &NNP);
  //interface for accessing data
  virtual const std::vector<Eigen::VectorXd>& getInputs(){return inputs;};
  virtual const std::vector<Eigen::VectorXd>& getActs(){return activations;};
  virtual int getNumPara();
  const std::vector<Eigen::VectorXd>& getDeltas(){return deltas;}
protected:
  double &actFunc;
  double &actFuncPrime;
  int numPara;
	std::vector<Eigen::VectorXd> inputs;
	std::vector<Eigen::VectorXd> activations;
  //z=w*input+bias, which is an intermediate variable but is needed
  //during the backpropagation step
	std::vector<Eigen::VectorXd> z;
  std::vector<Eigen::VectorXd> deltas;
};

#endif
