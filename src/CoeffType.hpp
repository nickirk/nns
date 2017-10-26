/*
 * CoeffType.hpp
 *
 *  Created on: Oct 25, 2017
 *      Author: guther
 */

#ifndef SRC_COEFFTYPE_HPP_
#define SRC_COEFFTYPE_HPP_


#include <complex>
#include "Eigen/Dense"

// We eventually want to use complex, but the current version uses VectorXd
//using coeffType = std::complex<double>;
// Only keep it as long as nececcary in this stage
// Typedef for coefficients
using coeffType = Eigen::VectorXd;

#endif /* SRC_COEFFTYPE_HPP_ */
