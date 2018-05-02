/*
 * ADAM.cxx
 *
 *  Created on: May 2, 2018
 *      Author: guther
 */

#include "ADAM.hpp"

namespace networkVMC {

ADAM::ADAM(double learningRate_):Solver(learningRate_),numPars(0),
		iteration(0),beta1(0.9),beta2(0.999) {
	m  = Eigen::VectorXd::Zero(numPars);
	m1  = Eigen::VectorXd::Zero(numPars);
	v  = Eigen::VectorXd::Zero(numPars);
	v1  = Eigen::VectorXd::Zero(numPars);
}

ADAM::~ADAM() {
}

void ADAM::update(VecType &w, VecType const &force, State const &input){
	// in the first iteration, we learn about the number of parameters
	if(iteration==0){
		numPars = w.size();
		m  = Eigen::VectorXd::Zero(numPars);
		m1  = Eigen::VectorXd::Zero(numPars);
		v  = Eigen::VectorXd::Zero(numPars);
		v1  = Eigen::VectorXd::Zero(numPars);
	}
	iteration+=1;
    m1 = beta1*m + (1-beta1)*force;
    v1 = beta2*v + (1-beta2)*force.array().square().matrix();
    m = m1;
    v = v1;
    m1 = m1/(1-std::pow(beta1,iteration));
    v1 = v1/(1-std::pow(beta2,iteration));
    w -= learningRate * (m1.array()/(v1.array().sqrt()+1e-8)).matrix();
}

} /* namespace networkVMC */
