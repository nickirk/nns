/*
 * ExcitmatType.hpp
 *
 *  Created on: Jun 20, 2018
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_EXCITATIONGENERATORS_EXCITMATTYPE_HPP_
#define SRC_HAMILTONIAN_EXCITATIONGENERATORS_EXCITMATTYPE_HPP_

namespace networkVMC {

// type of excitation matrices: copyable 2x2 array
class ExcitmatType{
public:
	// we store the entries contiguously
	ExcitmatType():
		entries(std::vector<int>(4,-1)){};
	// access operators (const/non-const version)
	// do some profiling to see if this slows us down
	int const& operator()(int i, int j) const{return entries[j+2*i];}
	int& operator()(int i, int j) {return const_cast<int&>(
			static_cast<ExcitmatType const&>(*this)(i,j));}
	// static cast is compile-time, so we do not expect a performance
	// effect
private:
	// array of the entries
	std::vector<int> entries;
};

} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_EXCITATIONGENERATORS_EXCITMATTYPE_HPP_ */
