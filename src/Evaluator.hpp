/*
 * Evaluator.hpp
 *
 *  Created on: Oct 25, 2017
 *      Author: guther
 */

#ifndef SRC_EVALUATOR_HPP_
#define SRC_EVALUATOR_HPP_

#include "CostFunction.hpp"
#include "State.hpp"

//Wrapper class to evaluate cost functions. Only has the purpose to prevent redundant calls
//to calc() by storing a flag if the function value has already been computed for this input

class Evaluator{
public:
	Evaluator(CostFunction const &cf_, State const &input_):cf(&cf_), input(&input_), unEvaluated(true){};
	double operator()() const{if(unEvaluated) cf->calc(*input); return cf->getValue();}
	std::vector<coeffType > nabla() const{return cf->nabla(*input); unEvaluated = false;}
private:
	// We use pointers here because assignment shall be possible
	CostFunction const *cf;
	State const *input;
	mutable bool unEvaluated;
};



#endif /* SRC_EVALUATOR_HPP_ */
