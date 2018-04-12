/*
 * TypeDefine.hpp
 *
 *  Created on: Oct 25, 2017
 *      Author: guther
 */

#ifndef SRC_TYPEDEFINE_HPP_
#define SRC_TYPEDEFINE_HPP_


#include <complex>
#include <vector>
#include "Eigen/Dense"

// We eventually want to use complex, but the current version uses VectorXd
//using coeffType = std::complex<double>;
// Only keep it as long as nececcary in this stage
// Typedef for coefficients
using coeffType = std::complex<double>;
using weightType = std::vector<std::vector<Eigen::Map<Eigen::MatrixXd>>>;
using biasType = std::vector<Eigen::Map<Eigen::VectorXd>>;
#endif /* SRC_TYPEDEFINE_HPP_ */
