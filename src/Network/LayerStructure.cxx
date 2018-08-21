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

namespace networkVMC{

template <typename F, typename coeffType>
LayerStructure<F, coeffType>::LayerStructure(){
	layers.clear();
}

//---------------------------------------------------------------------------//

// copy constructor
template <typename F, typename coeffType>
LayerStructure<F, coeffType>::LayerStructure(LayerStructure<F, coeffType> const &source){
	copyLayerStructure(source);
}

// move constructor
template <typename F, typename coeffType>
LayerStructure<F, coeffType>::LayerStructure(LayerStructure<F, coeffType> &&source){
	// take over the pointers from the source - we now own source's layers
	layers=source.layers;
	// and clear the source layers
	source.layers.clear();
}

// assignment operator
template <typename F, typename coeffType>
LayerStructure<F, coeffType>& LayerStructure<F, coeffType>::operator=(LayerStructure<F, coeffType> const &source){
	// check for self-assignment: we dont need to do anything in that case
	if(this!=&source){
		// else, call the copy functionality
		release();
		copyLayerStructure(source);
	}
	return *this;

}

// move assignment operator
template <typename F, typename coeffType>
LayerStructure<F, coeffType>& LayerStructure<F, coeffType>::operator=(LayerStructure<F, coeffType> &&source){
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
template <typename F, typename coeffType>
LayerStructure<F, coeffType>::~LayerStructure(){
	// use the release() function, which frees all memory of the layers
	release();
}

//---------------------------------------------------------------------------//

template <typename F, typename coeffType>
void LayerStructure<F, coeffType>::release(){
	// deallocate all layers
	for(size_t i{0};i<layers.size();++i){
		delete layers[i];
	}
}

//---------------------------------------------------------------------------//

// This is the body of copy construction/assignment
template <typename F, typename coeffType>
void LayerStructure<F, coeffType>::copyLayerStructure(LayerStructure<F, coeffType> const &source){
	layers.resize(source.layers.size());
	for(size_t i{0}; i<layers.size();++i){
		// the type is taken from the source via the virtual clone()
		// (it returns a pointer to a new object)
		layers[i] = source.layers[i]->clone();
	}
}

//---------------------------------------------------------------------------//

// add the different kinds of layers to the structure

template <typename F, typename coeffType>
void LayerStructure<F, coeffType>::addInputLayer(std::vector<Eigen::VectorXd> const &feedIns, int numStates){
	layers.push_back(new InputLayer(feedIns,numStates));
}

//---------------------------------------------------------------------------//
template <typename F, typename coeffType>
void LayerStructure<F, coeffType>::addDenseLayer(std::vector<T> const &input,
		std::string const &actFunc, int size){
	layers.push_back(new DenseLayer(input,actFunc,size));
}

//---------------------------------------------------------------------------//
template <typename F, typename coeffType>
void LayerStructure<F, coeffType>::addConvLayer(std::vector<T> const &inputs, std::string const &actFunc,
		int numFilters, int lengthFilter,int stride){
	layers.push_back(new ConvLayer(inputs,actFunc,numFilters,lengthFilter,stride));
}

}


