/*
 * LayerStructure.hpp
 *
 *  Created on: Apr 25, 2018
 *      Author: guther
 */

#ifndef SRC_HEADERS_LAYERSTRUCTURE_HPP_
#define SRC_HEADERS_LAYERSTRUCTURE_HPP_

#include <vector>
#include "../utilities/TypeDefine.hpp"
#include "../utilities/MatrixTypeDefine.hpp"
#include "Layers/Layer.hpp"

namespace networkVMC{

// The array containing the layers
// I externalized this to circumvent the need of RAII overhead
class LayerStructure {
  public:
	LayerStructure();
	// do all the nasty RAII stuff here instead of the NNW
	LayerStructure(LayerStructure const &source);
	LayerStructure(LayerStructure &&source);
	LayerStructure& operator=(LayerStructure const &source);
	LayerStructure& operator=(LayerStructure &&source);
	virtual ~LayerStructure();
	// Give it the same functionality as std::vector, so it feels like a std::vector
	size_t size() const {return layers.size();}
	// Except it does not return the objects themselves, but their addresses
	// (in fact, it is a vector of pointers and the assignment/construction handles the
	//  resource ownership, so this is a vector of polymorphic objects)
	Layer* operator[](int i) {return const_cast<Layer*>(static_cast<LayerStructure const&>(*this)[i]);}
	Layer const* operator[](int i) const{return layers[i];}
	// The push_back method might not be so commonly required, as we already have
	// explicit addXYLayer functions
	void push_back(Layer* &input){layers.push_back(input);}
	// Even though it also has a push_back method, I though it
	// easier to maintain with these functions
	void addInputLayer(std::vector<Eigen::VectorXd> const &feedIns, int numStates);
	void addDenseLayer(std::vector<paraVector> const &input, std::string const &actFunc, int size);
	void addConvLayer(std::vector<paraVector> const &inputs, std::string const &actFunc, int numFilters,
    int lengthFilter,int stride);
private:
	std::vector<Layer*> layers;
	// Internal directives for copying the vector
	void copyLayerStructure(LayerStructure const &source);
	// free all memory
	void release();
};

}

#endif /* SRC_HEADERS_LAYERSTRUCTURE_HPP_ */
