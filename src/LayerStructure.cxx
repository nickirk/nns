/*
 * LayerStructure.cxx
 *
 *  Created on: Apr 25, 2018
 *      Author: guther
 */

#include "LayerStructure.hpp"
#include "InputLayer.hpp"
#include "DenseLayer.hpp"
#include "ConvLayer.hpp"

LayerStructure::LayerStructure(){
	layers.clear();
}

//---------------------------------------------------------------------------//

void LayerStructure::addInputLayer(std::vector<Eigen::VectorXd> const &feedIns, int numStates){
	layers.push_back(InputLayer(feedIns,numStates));
}

//---------------------------------------------------------------------------//
void LayerStructure::addDenseLayer(std::vector<Eigen::VectorXd> const &input,
		std::string const &actFunc, int size){
	layers.push_back(DenseLayer(input,actFunc,size));
}

//---------------------------------------------------------------------------//
void LayerStructure::addConvLayer(std::vector<Eigen::VectorXd> const &inputs, std::string const &actFunc,
		int numFilters, int lengthFilter,int stride){
	layers.push_back(ConvLayer(inputs,actFunc,numFilters,lengthFilter,stride));
}



