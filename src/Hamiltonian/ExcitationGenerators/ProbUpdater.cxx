/*
 * probUpdater.cxx
 *
 *  Created on: Jun 19, 2018
 *      Author: Lauretta Schwarz, guther
 */

#include "ProbUpdater.hpp"
#include "ExcitmatType.hpp"
#include <cmath>
#include <iostream>
#include <omp.h>

namespace networkVMC {

// here we store the probabilities obtained on a thread globally
std::vector<double> ExcitationGenerator::ProbUpdater::threadedPParallel = std::vector<double>(omp_get_num_threads(),0.0);
std::vector<double> ExcitationGenerator::ProbUpdater::threadedPDoubles = std::vector<double>(omp_get_num_threads(),0.0);
// this marks if we already set the threaded probabilities on this thread
std::vector<bool> ExcitationGenerator::ProbUpdater::threadedPSet = std::vector<bool>(omp_get_num_threads(),false);

ExcitationGenerator::ProbUpdater::ProbUpdater() {
}

ExcitationGenerator::ProbUpdater::~ProbUpdater() {
}

void ExcitationGenerator::ProbUpdater::setThreadP(){
	int threadID = omp_get_thread_num();
	threadedPParallel[threadID] = pParallelInternal;
	threadedPParallel[threadID] = pDoublesInternal;
	if(not threadedPSet[threadID]) threadedPSet[threadID] = true;
};

void ExcitationGenerator::ProbUpdater::readThreadP(){
	int threadID = omp_get_thread_num();
	if(threadedPSet[threadID]){
		pParallelInternal = threadedPParallel[omp_get_thread_num()];
		pDoublesInternal = threadedPDoubles[omp_get_thread_num()];
	}
}

//---------------------------------------------------------------------------------------------------//

void ExcitationGenerator::ProbUpdater::setProbabilities(
		detType example_det){
    // set appropriate initial values for the probabilities based
    // on an example determinant

	// number of alpha, beta spin electrons
    int nel{0},nalpha{0},nbeta{0},norbs{0};

    // number of spin orbitals
    norbs = example_det.size();

    for (size_t i=0; i<example_det.size(); ++i){
        if (example_det[i]){
            // electron
            nel += 1;
            if ((i%2)==0){
                // alpha spin
                nalpha += 1;
            }
            else{
                // beta spin
                nbeta += 1;
            }
        }
    }

    // number of parallel spin pairs
    int parallel_pairs = (nalpha*(nalpha-1)/2) + (nbeta*(nbeta-1)/2);
    // number of opposite spin pairs
    int opp_pairs = nalpha*nbeta;

    // p_parallel + p_opposite = 1
    pParallelInternal = static_cast<double>(parallel_pairs)/static_cast<double>(parallel_pairs+opp_pairs);

    // number of single excitations
    int nsingleexcit = (nalpha*((norbs/2)-nalpha)) + (nbeta*((norbs/2)-nbeta));

    // number of double excitations
    // aa-pairs + bb-pairs + ab-pairs
    int ndoubleexcit = ((nalpha*(nalpha-1)/2)*(((norbs/2)-nalpha)*((norbs/2)-nalpha-1)/2)) +\
                       ((nbeta*(nbeta-1)/2)*(((norbs/2)-nbeta)*((norbs/2)-nbeta-1)/2)) +\
                       (nalpha*nbeta*((norbs/2)-nalpha)*((norbs/2)-nbeta));

    // p_singles + p_doubles = 1
    double pSingles = static_cast<double>(nsingleexcit)/
    		static_cast<double>(nsingleexcit+ndoubleexcit);
    pDoublesInternal = 1.0 - pSingles;

    std::cout << "\n Setting the probabilities to their initial values: \n" << std::endl;
    std::cout << "p_singles: " << pSingles << std::endl;
    std::cout << "p_doubles: " << pDoublesInternal << std::endl;
    std::cout << "p_parallel: " << pParallelInternal << std::endl;

}

//---------------------------------------------------------------------------------------------------//

void ExcitationGenerator::ProbUpdater::setProbabilitiesParallelBias(
		detType example_det, int exflag){
    // set appropriate initial values for the probabilities based
    // on an example determinant

    int nel,nalpha,nbeta,norbs;

    // number of spin orbitals
    norbs = example_det.size();

    bbiasSd = false;
    bbiasPo = false;

    if (exflag==1){
        // no adjustment needed since only single excitations a generated
        minsingle = 0;
        mindouble = 0;
        minparadouble = 0;
        minoppdouble = 0;
        bbiasSd = false;
        bbiasPo = false;
    }
    else{
        // adjust probabilities

        // number of alpha, beta spin electrons
        nel = 0;
        nalpha = 0;
        nbeta = 0;
        for (size_t i=0; i<example_det.size(); ++i){
            if (example_det[i]){
                // electron
                nel += 1;
                if ((i%2)==0){
                    // alpha spin
                    nalpha += 1;
                }
                else{
                    // beta spin
                    nbeta += 1;
                }
            }
        }

        // number of parallel spin pairs
        int parallel_pairs = (nalpha*(nalpha-1)/2) + (nbeta*(nbeta-1)/2);
        // number of opposite spin pairs
        int opp_pairs = nalpha*nbeta;

        // p_parallel + p_opposite = 1
        pParallelInternal = static_cast<double>(parallel_pairs)/static_cast<double>(parallel_pairs+opp_pairs);

        // number of single excitations
        int nsingleexcit = (nalpha*((norbs/2)-nalpha)) + (nbeta*((norbs/2)-nbeta));

        // number of double excitations
        // aa-pairs + bb-pairs + ab-pairs
        int ndoubleexcit = ((nalpha*(nalpha-1)/2)*(((norbs/2)-nalpha)*((norbs/2)-nalpha-1)/2)) +\
                           ((nbeta*(nbeta-1)/2)*(((norbs/2)-nbeta)*((norbs/2)-nbeta-1)/2)) +\
                           (nalpha*nbeta*((norbs/2)-nalpha)*((norbs/2)-nbeta));

        // at least 10 percent of the total number of excitation should be generated
        minsingle = static_cast<int>(0.1*static_cast<double>(nsingleexcit));
        mindouble = static_cast<int>(0.1*static_cast<double>(ndoubleexcit));
        minparadouble = static_cast<int>(0.1*static_cast<double>(parallel_pairs));
        minoppdouble = static_cast<int>(0.1*static_cast<double>(opp_pairs));

        if (exflag == 2){
            bbiasSd = false;
            bbiasPo = true;
            if ((parallel_pairs==0) or (opp_pairs==0)){
                bbiasPo = false;
                bbiasSd = false;
            }
        }
        else if (exflag == 3){
            bbiasPo = true;
            bbiasSd = true;
            if ((parallel_pairs==0) or (opp_pairs==0)){
                bbiasPo = false;
            }
            if ((nsingleexcit==0) or (ndoubleexcit==0)){
                bbiasSd = false;
            }

        }

        if (bbiasSd){
            std::cout << "\n Biasing p_singles / p_doubles" << std::endl;
            std::cout << "Minimum values: " << minsingle << " " << mindouble << "\n" << std::endl;
        }
        if (bbiasPo){
            std::cout << "\n Biasing p_parallel / p_oppposite" << std::endl;
            std::cout << "Minimum values: " << mindouble << " " << minparadouble << " " << minoppdouble << "\n" << std::endl;
        }
    }
}

//---------------------------------------------------------------------------------------------------//

void ExcitationGenerator::ProbUpdater::checkProbabilities(
		ExcitmatType const &excitmat, double hel, int nspawns, double pGen){
    // for updating these parameters

    double prob_ratio = 0.0;
    double hel_ratio = 0.0;
    double pSingles = 1.0 - pDoublesInternal;

    if (excitmat(1,0)<0){
        // single excitation
        prob_ratio = pGen/pSingles;
        hel_ratio = std::fabs(hel)/(prob_ratio*static_cast<double>(nspawns));
        max_hel_single_ratio = std::max(max_hel_single_ratio,hel_ratio);
        nsingle += 1;
    }
    else{
        // double excitation
        prob_ratio = pGen/pDoublesInternal;
        hel_ratio = std::fabs(hel)/(prob_ratio*static_cast<double>(nspawns));;
        max_hel_double_ratio = std::max(max_hel_double_ratio,hel_ratio);
        ndouble += 1;

        // distinguish between opposite and parallel spin
        if ((excitmat(0,0)%2)==(excitmat(1,0)%2)){
            // parallel spin
            prob_ratio /= pParallelInternal;
            hel_ratio = std::fabs(hel)/prob_ratio;

            max_hel_para_double_ratio = std::max(max_hel_para_double_ratio,hel_ratio);
            nparadouble += 1;
        }
        else{
            // opposite spin
            prob_ratio /= (1.0-pParallelInternal);
            hel_ratio = std::fabs(hel)/prob_ratio;

            max_hel_opp_double_ratio = std::max(max_hel_opp_double_ratio,hel_ratio);
            noppdouble += 1;
        }
    }
}

void ExcitationGenerator::ProbUpdater::adjustProbabilities(){
    // adjust the probabilities

    double pSingles_new,pDoubles_new,pParallel_new;

    double pSingles = 1.0 - pDoublesInternal;
    pSingles_new = pSingles;
    pDoubles_new = pDoublesInternal;
    pParallel_new = pParallelInternal;

    // check that there have been enough excitations
    // of each type
    bool enough_single = false;
    bool enough_double = false;
    bool enough_paradouble = false;
    bool enough_oppdouble = false;
    if (nsingle > minsingle){
        enough_single = true;
    }
    else if (bbiasPo and (not bbiasSd)){
        enough_single = true;
    }
    if (ndouble > mindouble){
        enough_double = true;
    }
    if (nparadouble > minparadouble){
        enough_paradouble = true;
    }
    if (noppdouble > minoppdouble){
        enough_oppdouble = true;
    }

    if (bbiasPo){
        if (enough_single and enough_double){
            pParallel_new = max_hel_para_double_ratio / \
                            (max_hel_para_double_ratio + max_hel_opp_double_ratio);
            pSingles_new = max_hel_single_ratio*pParallel_new / \
                           (max_hel_single_ratio*pParallel_new + max_hel_para_double_ratio);
            pDoubles_new = 1.0 - pSingles_new;
        }
        else{
            pParallel_new = pParallelInternal;
            pSingles_new = pSingles;
            pDoubles_new = 1.0 - pSingles_new;
        }

        if (enough_single and enough_double and bbiasSd){
            if ((std::fabs(pSingles_new-pSingles)/pSingles) > 1e-4){
                pSingles = pSingles_new;
                pDoublesInternal = pDoubles_new;
            }
        }

        if (enough_paradouble and enough_oppdouble and bbiasPo){
            if ((std::fabs(pParallel_new - pParallelInternal)/pParallelInternal) > 1e-4){
                pParallelInternal = pParallel_new;
            }
        }
    }
    else if (bbiasSd and (not bbiasPo)){
        if (enough_single and enough_double){
            pSingles_new = max_hel_single_ratio / \
                           (max_hel_single_ratio + max_hel_double_ratio);
            pDoubles_new = 1.0 - pSingles_new;
        }
        else{
            pSingles_new = pSingles;
            pDoubles_new = 1.0 - pSingles_new;
        }

        if (enough_single and enough_double and bbiasSd){
            if ((std::fabs(pSingles_new-pSingles)/pSingles) > 1e-4){
                pSingles = pSingles_new;
                pDoublesInternal = pDoubles_new;
            }
        }
    }


    max_hel_single_ratio = 0.0;
    max_hel_double_ratio = 0.0;
    max_hel_para_double_ratio = 0.0;
    max_hel_opp_double_ratio = 0.0;
    nsingle = 0;
    ndouble = 0;
    nparadouble = 0;
    noppdouble = 0;

    return;
}


} /* namespace networkVMC */
