/*
 * AbInitioHamiltonian.h
 *
 *  Created on: Feb 13, 2018
 *      Author: Lauretta Schwarz
 */

#ifndef SRC_HAMILTONIAN_ABINITIOHAMILTONIAN_HPP_
#define SRC_HAMILTONIAN_ABINITIOHAMILTONIAN_HPP_

#include <string>
#include "TwoBodyHamiltonian.hpp"
#include "../HilbertSpace/Determinant.hpp"

namespace networkVMC{

/**
 * \class AbInitoHamiltonian
 * \brief Fermionic generic two-body Hamiltonian
 *
 * Implements the framework for the ab-initio Hamiltonian in a Slater-Determinant basis, that is, we parametrize
 * the Hamiltonian as a generic fermionic two-body Hamiltonian. The 1-e and 2-e integrals are stored
 * and the Slater-Condon rules are used for evaluating matrix elements. This implies the dimension is not specified.
 */
class AbInitioHamiltonian: public TwoBodyHamiltonian {
  public:
	/// \param[in] dimension the number of orbitals (i.e. dimension of the single-particle Hilbert space
    AbInitioHamiltonian(int dimension):TwoBodyHamiltonian(dimension){};
	virtual ~AbInitioHamiltonian();
    /**
     * \brief Counts the number of states coupling to a given basis vector
     * \param[in] source basis vector to count from
     * \param exflag indicates the type of excitation:
     *               1 - single
     *               2 - double
     *               3 - both, first singles, then doubles
     * \param[out] nsingleexcit number of single excitations from source
     * \param[out] ndoubleexcit number of double excitations from source
     * \return number of possible excitations from source
     */
    int countNumberCoupledStates(detType const &source, int exflag, int &nsingleexcit, int &ndoubleexcit);

    // deterministic excitgen
    // TODO: move to its own class
    // deterministic excitation generator: generates all connected determinants
    std::vector<detType> getCoupledStates(detType const &source) const;
    // calculate the generation probability biased according to the connecting
    // Hamiltonian matrix element

    /// \return AbInitio
    virtual HType type() const {return AbInitio;}

    // The fermi sign is the main difference to the plain TwoBodyHamiltonian
    /// for AbInitioHamiltonian, the sign created by an excitation operator is given by the fermi sign
    int getFermiSign(detType const &alpha, int annihilatorIndex, int creatorIndex) const{
    	return JWStringLength(alpha,annihilatorIndex,creatorIndex);
    }
  private:
    /**
     * \brief Deterministically create a single excitation from source
     * \param[in] source determinant from which to excite
     * \param[in] source_orbs orbitals occupied in source
     * \param holes holes of the last excitation
     * \param particles electrons of the last excitation
     * \param exflag indicates the type of excitation:
     *               1 - single
     *               2 - double
     *               3 - both, first singles, then doubles
     * \param[out] ballexcitfound flag indicating if all excitations have been generated
     * \return basis vector representing the next excitation of source
     *
     * This (and getDoubleExcitation) are used internally by getCoupledStates to generate a list
     * of all excitations from source. They are called repeatedly using the last output to supply the input arguments
     * holes and particles for the next call. This way, all excitations are generated.
     */
    detType getSingleExcitation(detType const &source,
    		std::vector<int> const &source_orbs, std::vector<int> &holes,
			std::vector<int> &particles, int &exflag, bool &ballexcitfound) const;

    /**
     * \brief Deterministically generate a double excitation from source
     * \param[in] source determinant from which to excite
     * \param[in] source_orbs orbitals occupied in source
     * \param holes holes of the last excitation
     * \param particles electrons of the last excitation
     * \param exflag indicates the type of excitation:
     *               1 - single
     *               2 - double
     *               3 - both, first singles, then doubles
     * \param[out] ballexcitfound flag indicating if all excitations have been generated
     * \return basis vector representing the next excitation of source
     *
     * See getSingleExcitation
     */
    detType getDoubleExcitation(detType const &source, std::vector<int>
    const &source_orbs, std::vector<int> &holes, std::vector<int> &particles,
	int &exflag, bool &ballexcitfound) const;

    /**
     * \brief Pick a pair of electrons
     * \param[in] source_orbs list of orbitals occupied in source
     * \param[out] el1 index of the first electron
     * \param[out] el2 index of the second electron
     * \param[out] spinpair indicates the spin of the two electrons
     *           spinpair = 1 -> alpha + alpha
     *			 spinpair = 2 -> alpha + beta
     * 			 spinpair = 3 -> beta + alpha
     *			 spinpair = 4 -> beta + beta
     * \param[in] ind index of the electron pair in a list of all electron pairs
     *
     * This is called internally in the deterministic excitation generation to loop over all
     * pairs of electrons
     */
    void getElecPair(std::vector<int> const &source_orbs,
            int &el1, int &el2, int &spinpair, int ind) const;
};

AbInitioHamiltonian readAbInitioHamiltonian(std::string file_name, bool molpro_fcidump=false);

}

#endif /* SRC_HAMILTONIAN_ABINITIOHAMILTONIAN_HPP_ */
