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

FileNotFound::FileNotFound(std::string file_):file(file_){
    std::cout << "File " << file << " not found" << std::endl;
}

}

