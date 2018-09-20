/*
 * MatrixTypeDefine.hpp
 *
 *  Created on: Sep 20, 2018
 *      Author: guther
 */

#ifndef SRC_UTILITIES_MATRIXTYPEDEFINE_HPP_
#define SRC_UTILITIES_MATRIXTYPEDEFINE_HPP_

#include <Eigen/Dense>
#include "TypeDefine.hpp"

namespace networkVMC{
	using paraVector = Eigen::Matrix<paraType, Eigen::Dynamic, 1>;
}



#endif /* SRC_UTILITIES_MATRIXTYPEDEFINE_HPP_ */
