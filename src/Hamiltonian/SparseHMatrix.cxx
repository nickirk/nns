/*
 * SparseHMatrix.cxx
 *
 *  Created on: Apr 10, 2018
 *      Author: guther
 */

#include <cmath>
#include "../math/constants.hpp"
#include "SparseHMatrix.hpp"

namespace networkVMC{
SparseHMatrix::~SparseHMatrix() {}

// This is what the arpack++ interface takes
// Here, we just have to rely on in/out being of size dim
void SparseHMatrix::MatMul(coeffType *in, coeffType *out){
// Do the Matrix-Vector multiplication
	// first, set the output vector to 0
	for(size_t i{0};i<dim;++i){
		out[i] = 0;
	}
	// then, do the sparse matrix-vector mutltiplication
	for(size_t i{0};i<rowPos.size()-1;++i){
		// for each row (starting with 0), get all the corresponding contributions
		for(size_t j{rowPos[j]};j<rowPos[i+1];++j){
			// entries[j] is the contribution to row i from the entry of in at cols[j]
			out[i] += entries[j]*in[cols[j]];
		}
	}
}

// This sets up the sparse matrix representing H in subspace
void SparseHMatrix::load(Hamiltonian const &H, State const &subspace){
	dim = subspace.size();
	entries.clear();
	rowPos.clear();
	cols.clear();
// scan the Hamiltonian for nonzero entries and store them
	for(size_t i{0};i<dim;++i){
		rowPos.push_back(entries.size());
		for(size_t j{0};j<dim;++j){
			double tmp{H(subspace.det(i),subspace.det(j))};
			if(std::abs(tmp) > epsilon){
				entries.push_back(tmp);
				cols.push_back(j);
			}
		}
	}
// by convention and for convenience, add the number of nonzero entries as a last entry to rowPos
	rowPos.push_back(entries.size());
}

}

