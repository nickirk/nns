/*
 * ExcitmatType.hpp
 *
 *  Created on: Jun 20, 2018
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_EXCITATIONGENERATORS_EXCITMATTYPE_HPP_
#define SRC_HAMILTONIAN_EXCITATIONGENERATORS_EXCITMATTYPE_HPP_

namespace networkVMC {

/**
 * \class ExcitmatType
 * \brief type of excitation matrices: deep-copyable 2x2 array
 */
class ExcitmatType{
public:
	// we store the entries contiguously
	ExcitmatType():
		entries(std::vector<int>(4,-1)){};
	// access operators (const/non-const version)
	// do some profiling to see if this slows us down

	/**
	 * \brief Get an entry of the excitation matrix
	 * \param[in] i,j indices of the excitation matrix
	 * \return entry of the (2,2) excitation matrix at (i,j)
	 */
	int const& operator()(int i, int j) const{return entries[j+2*i];}
	/// \overload
	int& operator()(int i, int j) {return const_cast<int&>(
			static_cast<ExcitmatType const&>(*this)(i,j));}
	// static cast is compile-time, so we do not expect a performance
	// effect
private:
	/// contiguous array of the entries
	std::vector<int> entries;
};

} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_EXCITATIONGENERATORS_EXCITMATTYPE_HPP_ */
