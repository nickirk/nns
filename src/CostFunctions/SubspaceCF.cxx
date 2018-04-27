/*
 * SubspaceCF.cxx
 *
 *  Created on: Apr 10, 2018
 *      Author: guther
 */

#include "SubspaceCF.hpp"

#include "../../lib/arpackpp/include/arssym.h"
#include "NormCF.hpp"

namespace networkVMC{

SubspaceCF::~SubspaceCF() {
}

std::vector<Eigen::VectorXd > SubspaceCF::nabla(State const &input) const{
	NormCF dist(diagonalizeSubspace(input));
	return dist.nabla(input);
}

State SubspaceCF::diagonalizeSubspace(State const & input) const{
	return input;
}

}
