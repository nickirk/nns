/*
 * Layer.hpp
 * Created on 06.3.2018
 * Author: Ke Liao 
 */

#ifndef Layer_DEFINED
#define Layer_DEFINED

#include <vector>
#include <Eigen/Dense>

class Layer{
public:
	Layer(std::vector<Eigen::MatrixXd> const &inputs_, double &(actFunc_)(double))
	inputs(inputs_), actFunc(actFunc_);
	virtual ~Layer();
	double &actFunc();
  virtual void processSignal();
  virtual void updatePara();
  virtual std::vector<Eigen::MatrixXd> getInputs(){return inputs;};
  virtual std::vector<Eigen::MatrixXd> getActs(){return activations;};
  virtual void initialise(Eigen::VectorXd const &NNP);
protected:
	std::vector<Eigen::MatrixXd> inputs;
	std::vector<Eigen::MatrixXd> activations;
}

#endif
