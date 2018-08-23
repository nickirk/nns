/*
 * TypeDefine.hpp
 *
 *  Created on: Oct 25, 2017
 *      Author: guther
 */

#ifndef SRC_UTILITIES_TYPEDEFINE_HPP_
#define SRC_UTILITIES_TYPEDEFINE_HPP_


#include <complex>
#include <vector>

namespace networkVMC{

// We eventually want to use complex, but the current version uses VectorXd
// Only keep it as long as nececcary in this stage
// Typedef for coefficients
using cType = std::complex<double>;
using detType = std::vector<bool>;

// Supertypes for samplers and Hamiltonians to set defaults
// and avoid invalid combinations
enum SamplerType {Markov, PreFetched};
enum HType {Hubbard, AbInitio, Constant, Heisenberg};

}
#endif /* SRC_UTILITIES_TYPEDEFINE_HPP_ */
