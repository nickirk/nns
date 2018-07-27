/*
 * FermionBasis.hpp
 *
 *  Created on: Jul 26, 2018
 *      Author: guther
 */

#ifndef SRC_HILBERTSPACE_FERMIONBASIS_HPP_
#define SRC_HILBERTSPACE_FERMIONBASIS_HPP_

#include "Basis.hpp"
#include "../utilities/SpinConfig.hpp"

namespace networkVMC {

class FermionBasis: public Basis {
public:
	FermionBasis(SpinConfig const &spinConfig_);
	// Return the determinant with index 'index'
	detType getDetByIndex(int index) const; // can throw an OutOfRangeError
	// Return the index of the determinant 'det_'
	int getIndexByDet(detType const & det_) const; // can throw an InvalidDeterminantError
	// Return the alpha/beta spin distribution
	SpinConfig const& getSpinConfig() const {return spinConfig;};
private:
  int numEle;
  SpinConfig spinConfig;
  int numOrb;
  int indexOfDet;
  std::vector<int> listOfOrbNum;
  std::vector<int> combination;
// internal methods for handling the map
  int calcSize(int numOrb_, int numEle_);
  void createFermionBasisDet(int offset, int numEle_);
  // determinants need to be accessed often, no need to add some
  // overhead by using an extra class here, better use an alias
  std::vector<detType > basis;
  std::vector<int> indexBasis;
};

} /* namespace networkVMC */

#endif /* SRC_HILBERTSPACE_FERMIONBASIS_HPP_ */
