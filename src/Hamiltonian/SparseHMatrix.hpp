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
#include "../utilities/State.hpp"

namespace networkVMC{

class SparseHMatrix {
// Class for translating a second quantized Hamiltonian into a sparse matrix
// This is not really storing the Hamiltonian, but applying it to a vector
public:
	SparseHMatrix();
	SparseHMatrix(Hamiltonian const &H,State const &subspace){load(H,subspace);};
	virtual ~SparseHMatrix();
// Multiplication function for the Arpack interface
	void MatMul(coeffType *in, coeffType *out);

// create a sparse representation of H in subspace
	void load(Hamiltonian const &H, State const &subspace);
private:
// We consider the Hamiltonian H in the space spanned by a the determinants of a given State subspace
	size_t dim;
// sparse storage scheme (CRS)
// entries are the values of the nonzero elements of H
	std::vector<double> entries;
// cols are the columns belonging to the values, and rowPos contains the indices of those entries
// that start a new row
	std::vector<int> rowPos,cols;
};



}

#endif /* SRC_HAMILTONIAN_SPARSEHMATRIX_HPP_ */
