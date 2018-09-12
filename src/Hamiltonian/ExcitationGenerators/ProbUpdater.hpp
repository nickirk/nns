/*
 * probUpdater.hpp
 *
 *  Created on: Jun 19, 2018
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_EXCITATIONGENERATORS_PROBUPDATER_HPP_
#define SRC_HAMILTONIAN_EXCITATIONGENERATORS_PROBUPDATER_HPP_

#include "../../utilities/TypeDefine.hpp"
#include "ExcitationGenerator.hpp"

namespace networkVMC {

class ExcitmatType;

/**
 * \class ExcitationGenerator::ProbUpdater
 * \brief Manages internal biases of ExcitationGenerators
 *
 * This is used to keep statistics of already generated excitations and re-adjust internal biases
 */
class ExcitationGenerator::ProbUpdater {
public:
	ProbUpdater();
	/// \param[in] exampleDet examplary basis vector to estimate initial values for the biases
	explicit ProbUpdater(detType const &exampleDet){
		setProbabilities(exampleDet);}
	virtual ~ProbUpdater();

	/**
	 * \brief getter for the probability to pick a spin-parallel excitation
	 * \return probability to pick a spin-parallel excitation
	 */
	virtual double pParallel() const {return pParallelInternal;}
	/**
	 * \brief getter for the probability to pick a double excitation
	 * \return probability to pick a double excitation
	 */
	virtual double pDoubles() const {return pDoublesInternal;}

    // set values for the generation probabilities based on the number
    // of electrons and holes
	void setProbabilities(detType example_det);
    // initialise these parameters
    void setProbabilitiesParallelBias(detType example_det,
    		int exflag);

    // for updating these parameters
    void checkProbabilities(ExcitmatType const &excitmat, double hel, int nspawns, double pGen);
    void adjustProbabilities();
private:
    /// This is what is accumulated accumulate during sampling
    struct gatheredData{
    	// collection of POD
        int nsingle;
        int ndouble;
        int nparadouble;
        int noppdouble;
        double max_hel_single_ratio;
        double max_hel_double_ratio;
        double max_hel_para_double_ratio;
        double max_hel_opp_double_ratio;

        /// default initialization: all data members 0
        gatheredData();
        /**
         * \brief accumulation: add up the integers, take the max of the floats
         * \param[in] rhs gatheredData to add
         * \return *this with accumulated members
         * Adds the integer members to those of *this and takes the max of the float members of rhs and *this
         */
        gatheredData& operator+=(gatheredData const &rhs);
    };

    /// the relevant probabilities: Biases for double and spin-parallel excitations
    double pDoublesInternal, pParallelInternal;
    // for the dynamic adjustment of the probabilities
    int minsingle;
    int mindouble;
    int minparadouble;
    int minoppdouble;
    /// statistics accumulated in this instance
    gatheredData localStatistics;
    bool bbiasSd;
    bool bbiasPo;

    /// for paralellization: thread-global containers for gathered data
    static gatheredData globalStatistics;
    /// accumulate across threads
    void globalizeMax();
};

void setNewBiases(ExcitationGenerator::ProbUpdater &pBiasGen, double &pDoubles, double &pParallel);

} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_EXCITATIONGENERATORS_PROBUPDATER_HPP_ */
