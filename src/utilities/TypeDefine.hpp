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
//using coeffType = std::complex<double>;
// Only keep it as long as nececcary in this stage
// Typedef for coefficients
using coeffType = std::complex<double>;
using detType = std::vector<bool>;
//numFilter<depthFilter<lengthFilter>>
// general eigen vector type for complex numbers, used for parameters etc
using VecType = Eigen::VectorXd;
using VecCType = Eigen::VectorXcd;
using weightType = std::vector<std::vector<Eigen::Map<Eigen::MatrixXd>>>;
//numFilter<lengthFilter>
using biasType = std::vector<Eigen::Map<VecType>>;
// return type for cost function nabla (derivative with respect to a vector)
using nablaType = std::vector<coeffType>;

// Supertypes for samplers and Hamiltonians to set defaults
// and avoid invalid combinations
enum SamplerType {Markov, PreFetched};
enum HType {Hubbard, AbInitio, Constant};

}
#endif /* SRC_UTILITIES_TYPEDEFINE_HPP_ */
