/*
 * SubspaceCF.cxx
 *
 *  Created on: Apr 10, 2018
 *      Author: guther
 */

#include "SubspaceCF.hpp"

#include "../../lib/arpackpp/include/arssym.h"
#include "../../lib/arpackpp/include/arscomp.h"
#include "NormCF.hpp"
#include "../Hamiltonian/SparseHMatrix.hpp"

namespace networkVMC{

SubspaceCF::~SubspaceCF() {
}

nablaType SubspaceCF::nabla(State const &input) const{
  auto dist = NormCF(diagonalizeSubspace(input));
  distance = dist.calc(input);
  return dist.nabla(input);
}

//---------------------------------------------------------------------------------------------------//

State SubspaceCF::diagonalizeSubspace(State const & input) const{
	SparseHMatrix HMatrix(H, input);
	double tol = 1e-10;
	int maxIter = 40;
	// the output State looks the same as the input state
	State output{input};
	// use the input state's copy as initial state
	coeffType *outputVec{&(output.coeff(0))};
	// also store the subspace energy: Even though we cannot use it as a cost function
	// it gives a useful estimate of energy
	coeffType *energyPtr{&subspaceEnergy};
	// set up the eigenvalue problem
	ARCompStdEig<double, SparseHMatrix> eigProblem(HMatrix.dimension(),1,
			&HMatrix,&SparseHMatrix::MatMul,"SR",0,tol,maxIter,outputVec);
	// solve it and write the eigenvector to output
	int nconv = eigProblem.EigenValVectors(outputVec,energyPtr);
	if(nconv == 0) throw UnconvergedEigenproblem();
	return output;
}

}
