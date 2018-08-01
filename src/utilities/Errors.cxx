/*
 * Errors.cxx
 *
 *  Created on: Jul 27, 2018
 *      Author: guther
 */


#include "Errors.hpp"
#include "../HilbertSpace/Determinant.hpp"
#include <iostream>

namespace networkVMC{

InvalidDeterminantError::InvalidDeterminantError(detType const &a)
{
		std::cout << "Error: invalid determinant encountered:\n";
		printDet(a);
};

}

