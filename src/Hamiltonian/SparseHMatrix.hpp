/*
 * SparseHMatrix.hpp
 *
 *  Created on: Apr 10, 2018
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_SPARSEHMATRIX_HPP_
#define SRC_HAMILTONIAN_SPARSEHMATRIX_HPP_

#include <vector>

#include "../utilities/TypeDefine.hpp"
#include "Hamiltonian.hpp"

namespace networkVMC{

class SparseHMatrix {
// Class for translating a second quantized Hamiltonian into a sparse matrix
public:
	SparseHMatrix();
	SparseHMatrix(Hamiltonian const &H){load(H);};
	virtual ~SparseHMatrix();
// Multiplication function for the Arpack interface
	void MatMul(std::vector<coeffType> const &in, std::vector<coeffType> const &out);
// Load a Hamiltonian and construct the sparse matrix
	void load(Hamiltonian const &H);
private:
// use a standard sparse storage scheme
	std::vector<coeffType> entries;
	std::vector<int> rows, cols;
};

}

#endif /* SRC_HAMILTONIAN_SPARSEHMATRIX_HPP_ */
