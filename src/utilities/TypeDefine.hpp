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
#include <Eigen/Dense>

namespace networkVMC{

// We eventually want to use complex, but the current version uses VectorXd
// Only keep it as long as nececcary in this stage
// Typedef for coefficients
//using coeffType = std::complex<double>;
using detType = std::vector<bool>;
// the weightType and biasType will be eventually removed to be consistent with
// the template implementation.
//using weightType = std::vector<std::vector<Eigen::Map<Eigen::MatrixXd>>>;
//numFilter<lengthFilter>
//using biasType = std::vector<Eigen::Map<Eigen::VectorXd>>;
// return type for cost function nabla (derivative with respect to a vector)
//using nablaType = std::vector<coeffType>;

// Supertypes for samplers and Hamiltonians to set defaults
// and avoid invalid combinations
enum SamplerType {Markov, PreFetched};
enum HType {Hubbard, AbInitio, Constant, Heisenberg};

}
#endif /* SRC_UTILITIES_TYPEDEFINE_HPP_ */
