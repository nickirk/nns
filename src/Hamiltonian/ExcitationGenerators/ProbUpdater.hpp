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

class ExcitationGenerator::ProbUpdater {
public:
	ProbUpdater();
	explicit ProbUpdater(detType const &exampleDet){
		setProbabilities(exampleDet);}
	virtual ~ProbUpdater();

	virtual double pParallel() const {return pParallelInternal;}
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
    // This is what we accumulate during sampling
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

        // default initialization: all data members 0
        gatheredData();
        // accumulation: add up the integers, take the max of the floats
        gatheredData& operator+=(gatheredData const &rhs);
    };

    // the relevant probabilities
    double pDoublesInternal, pParallelInternal;
    // for the dynamic adjustment of the probabilities
    int minsingle;
    int mindouble;
    int minparadouble;
    int minoppdouble;
    gatheredData localStatistics;
    bool bbiasSd;
    bool bbiasPo;

    // for paralellization: thread-global containers for gathered data
    static gatheredData globalStatistics;
    void globalizeMax();
};

void setNewBiases(ExcitationGenerator::ProbUpdater &pBiasGen, double &pDoubles, double &pParallel);

} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_EXCITATIONGENERATORS_PROBUPDATER_HPP_ */
