/*
 * LayerStructure.cxx
 *
 *  Created on: Apr 25, 2018
 *      Author: guther
 */

#include "LayerStructure.hpp"

#include "Layers/ConvLayer.hpp"
#include "Layers/DenseLayer.hpp"
#include "Layers/InputLayer.hpp"

LayerStructure::LayerStructure(){
	layers.clear();
}

//---------------------------------------------------------------------------//

// copy constructor
LayerStructure::LayerStructure(LayerStructure const &source){
	copyLayerStructure(source);
}

// move constructor
LayerStructure::LayerStructure(LayerStructure &&source){
	// take over the pointers from the source - we now own source's layers
	layers=source.layers;
	// and clear the source layers
	source.layers.clear();
}

// assignment operator
LayerStructure& LayerStructure::operator=(LayerStructure const &source){
	// check for self-assignment: we dont need to do anything in that case
	if(this!=&source){
		// else, call the copy functionality
		release();
		copyLayerStructure(source);
	}
	return *this;

}

// move assignment operator
LayerStructure& LayerStructure::operator=(LayerStructure &&source){
	// again, check for self-move
	if(this!=&source){
		release();
		// take the layers of source
		layers=source.layers;
		// and clear the original ones
		source.layers.clear();
	}
	return *this;
}

//---------------------------------------------------------------------------//

// also free the layers upon destruction
LayerStructure::~LayerStructure(){
	// use the release() function, which frees all memory of the layers
	release();
}

//---------------------------------------------------------------------------//

void LayerStructure::release(){
	for(size_t i{0};i<layers.size();++i){
		delete layers[i];
	}
}

//---------------------------------------------------------------------------//

// This is the body of copy construction/assignment
void LayerStructure::copyLayerStructure(LayerStructure const &source){
	layers.resize(source.layers.size());
	for(size_t i{0}; i<layers.size();++i){
		// the type is taken from the source via the virtual clone()
		// (it returns a pointer to a new object)
		layers[i] = source.layers[i]->clone();
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



