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
    // the relevant probabilities
    double pDoublesInternal, pParallelInternal;
    // for the dynamic adjustment of the probabilities
    int minsingle;
    int mindouble;
    int minparadouble;
    int minoppdouble;
    double max_hel_single_ratio;
    double max_hel_double_ratio;
    double max_hel_para_double_ratio;
    double max_hel_opp_double_ratio;
    int nsingle;
    int ndouble;
    int nparadouble;
    int noppdouble;
    bool bbiasSd;
    bool bbiasPo;
};

} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_EXCITATIONGENERATORS_PROBUPDATER_HPP_ */
