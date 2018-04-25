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

// also free the layers upon destruction
LayerStructure::~LayerStructure(){
	for(size_t i{0};i<layers.size();++i){
		delete layers[i];
	}
}

//---------------------------------------------------------------------------//

void LayerStructure::addInputLayer(std::vector<Eigen::VectorXd> const &feedIns, int numStates){
	layers.push_back(new InputLayer(feedIns,numStates));
}

//---------------------------------------------------------------------------//
void LayerStructure::addDenseLayer(std::vector<Eigen::VectorXd> const &input,
		std::string const &actFunc, int size){
	layers.push_back(new DenseLayer(input,actFunc,size));
}

//---------------------------------------------------------------------------//
void LayerStructure::addConvLayer(std::vector<Eigen::VectorXd> const &inputs, std::string const &actFunc,
		int numFilters, int lengthFilter,int stride){
	layers.push_back(new ConvLayer(inputs,actFunc,numFilters,lengthFilter,stride));
}



