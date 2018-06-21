/*
 * WeightedExcitgen.cxx
 *
 *  Created on: Jun 15, 2018
 *      Author: Lauretta Schwarz, guther
 */

#include "WeightedExcitgen.hpp"
#include "WeightedSelector.hpp"
#include "../../HilbertSpace/Determinant.hpp"
#include "../../utilities/TypeDefine.hpp"
#include "../../utilities/Errors.hpp"
#include "../Hamiltonian.hpp"
#include <iostream>

namespace networkVMC {

WeightedExcitgen::WeightedExcitgen(Hamiltonian const &H_, detType const &HF):
	clonableExcitgen<WeightedExcitgen>(),pBiasGen(ProbUpdater(HF)),H(&H_){
	pParallel = pBiasGen.pParallel();
	pDoubles = pBiasGen.pDoubles();
}

//---------------------------------------------------------------------------------------------------//

WeightedExcitgen::~WeightedExcitgen() {
}

//---------------------------------------------------------------------------------------------------//

void WeightedExcitgen::constructClassCount(detType const &source) {
    // fill in information necessary for random excitation generation

    // get the occupied spin orbitals
    sourceOrbs = getOccupiedPositions(source);

    excitmat(0,0) = -1;
    excitmat(0,1) = -1;
    excitmat(1,0) = -1;
    excitmat(1,1) = -1;
    pgen = 0.0;

    nel = sourceOrbs.size();
    norbs = source.size();
    nalphaels = 0;
    nbetaels = 0;
    // number of occupied and unoccupied alpha and beta spin electrons
    for (size_t i=0; i<sourceOrbs.size(); ++i){
        if ((sourceOrbs[i]%2) == 0){
            // alpha spin
            nalphaels += 1;
        }
        else{
            // beta spin
            nbetaels += 1;
        }
    }
    nalphaholes = (source.size()/2) - nalphaels;
    nbetaholes = (source.size()/2) - nbetaels;

    // number of electron pairs
    aa_elec_pairs = nalphaels*(nalphaels-1)/2;
    bb_elec_pairs = nbetaels*(nbetaels-1)/2;
    ab_elec_pairs = nalphaels*nbetaels;

    bfilled = true;

}


//---------------------------------------------------------------------------------------------------//

detType WeightedExcitgen::generateExcitation(
			detType const &source, double &pGen){
	// generate a random excitation

	// set up the random number generator
	// initialise the seed engine
	std::random_device rd;
	// use Mersenne-Twister random number generator
	std::mt19937 rng(rd());
	// uniform distribution from 0.0 to 1.0
	std::uniform_real_distribution<double> uniform_dist(0.0,1.0);


	// get information on the number of alpha and beta soin
	// electrons and holes for the determinant
	constructClassCount(source);

	// return variable
	detType target = source;
	// decide whether a single or double excitation is done
	try{
	if (uniform_dist(rng) < 1 - pDoubles){
		// single excitation
		target = generateSingleExcit(source);
		pgen *= (1- pDoubles);
	}
	else{
		// double excitation
		target = generateDoubleExcit(source);
		pgen *= pDoubles;
	}
	}
	catch(NoExcitationFound const&){
		// is this legit? do not make a move if there is none
		target = source;
	}

	// assing the probability
	pGen = pgen;

	return target;
}

//---------------------------------------------------------------------------------------------------//

detType WeightedExcitgen::generateSingleExcit(detType const &source) {
    // generate a single excitation

    // set up the random number generator
    // initialise the seed engine
    std::random_device rd;
    // use Mersenne-Twister random number generator
    std::mt19937 rng(rd());
    // uniform distribution from 0.0 to 1.0
    std::uniform_real_distribution<double> uniform_dist(0.0,1.0);

    // pick an electron at random with uniform probability
    // floor: greatest integer <= x
    //int elec = 1 + std::floor((uniform_dist(rng)*nel));
    int elec = std::floor((uniform_dist(rng)*nel));
    int src = sourceOrbs[elec];
    // return variable: initialize it with the input
    auto target = source;
    // auxiliary selector
    WeightedSelector selector(*H,source);
    // select the target hole by connection strength
    int tgt = selector.selectSingleHole(src,pgen);
    // if we did not find a hole, throw an exception
    if (tgt < 0){
    	throw NoExcitationFound();
    }

    // generate the new determinant
    target = source;
    annihilate(target,src);
    create(target,tgt);
    excitmat(1,0) = -1;
    excitmat(1,1) = -1;
    excitmat(0,0) = src;
    excitmat(0,1) = tgt;

    pgen /= static_cast<double>(nel);

    return target;

}

//---------------------------------------------------------------------------------------------------//

detType WeightedExcitgen::generateDoubleExcit(detType const &source){
    // generate a double excitation

    std::vector<double> cpt_pair,int_cpt,sum_pair,cum_sum;
    std::vector<int> elecs;
    std::vector<int> tgt;
    tgt.assign(2,-1);

    // auxiliary selector for picking holes/electron
    WeightedSelector doublesSelector(*H,source);

    // the return value, initially just the input, we will then
    // do some manipulations with that and return it
    auto target = source;

    // select a pair of electrons
    std::vector<int> src = pickBiasedElecs(elecs, source);


    // select the two target holes by approximate connection
    // strength
    if ((not H->partExact()) and (not H->linExact()) and ((src[0]%2)==(src[1]%2))){
        // in this case the CDFs are the same
        std::vector<int> tgt =
        		doublesSelector.selectDoubleHoles(src,cum_sum,int_cpt);
    }
    else{
        //std::vector<int> tgt;
        //tgt.assign(2,-1);
        int_cpt.assign(2,0.0);
        cum_sum.assign(2,0.0);

        tgt[0] = doublesSelector.selectDoubleHole(
        		src,-1,cum_sum[0],int_cpt[0]);

        if (tgt[0] > -1){
            tgt[1] =doublesSelector.selectDoubleHole(
            				src,tgt[0],cum_sum[1],int_cpt[1]);
        }
    }

    // check that the two electrons have not been excited into the
    // same orbitals
    if (std::any_of(tgt.begin(),tgt.end(), [](const int i) {return i<0;})){
    	throw NoExcitationFound();
    }

    // adjust the probabilities: if the two holes are selected from
    // the same list, they could have been selected in any order
    cpt_pair.assign(2,0.0);
    sum_pair.assign(2,0.0);
    if ((tgt[0]%2)==(tgt[1]%2)){
        if (H->partExact() or H->linExact()){
        	cpt_pair[0] = pGenSelectDoubleHole(source,src,-1,sum_pair[0],tgt[1]);
            cpt_pair[1] = pGenSelectDoubleHole(source,src,tgt[1],sum_pair[1],tgt[0]);
        }
        else{
            cpt_pair[0] = int_cpt[1];
            cpt_pair[1] = int_cpt[0];
            sum_pair[0] = cum_sum[0];
            sum_pair[1] = cum_sum[0] - int_cpt[1];
        }
        //double prod_int_cpt = std::accumulate(int_cpt.begin(),int_cpt.end(),1,std::multiplies<double>());
        double prod_int_cpt = int_cpt[0]*int_cpt[1];
        //double prod_cum_sum = std::accumulate(cum_sum.begin(),cum_sum.end(),1,std::multiplies<double>());
        double prod_cum_sum = cum_sum[0]*cum_sum[1];
        //double prod_cpt_pair = std::accumulate(cpt_pair.begin(),cpt_pair.end(),1,std::multiplies<double>());
        double prod_cpt_pair = cpt_pair[0]*cpt_pair[1];
        //double prod_sum_pair = std::accumulate(sum_pair.begin(),sum_pair.end(),1,std::multiplies<double>());
        double prod_sum_pair = sum_pair[0]*sum_pair[1];
        pgen *= (prod_int_cpt/prod_cum_sum + prod_cpt_pair/prod_sum_pair);
    }
    else{
        //double prod_int_cpt = std::accumulate(int_cpt.begin(),int_cpt.end(),1,std::multiplies<double>());
        double prod_int_cpt = int_cpt[0]*int_cpt[1];
        //double prod_cum_sum = std::accumulate(cum_sum.begin(),cum_sum.end(),1,std::multiplies<double>());
        double prod_cum_sum = cum_sum[0]*cum_sum[1];
        pgen *= (prod_int_cpt/prod_cum_sum);
    }

    // generate the new determinant
    annihilate(target,src[0]);
    annihilate(target,src[1]);
    create(target,tgt[0]);
    create(target,tgt[1]);
    excitmat(0,0) = src[0];
    excitmat(1,0) = src[0];
    if ((src[0]%2)==(tgt[0]%2)){
        excitmat(0,1) = tgt[0];
        excitmat(1,1) = tgt[1];
    }
    else{
        excitmat(0,1) = tgt[1];
        excitmat(1,1) = tgt[0];
    }

    return target;

}

//---------------------------------------------------------------------------------------------------//

std::vector<int> WeightedExcitgen::pickBiasedElecs(std::vector<int> &elecs, detType const &source){
    // pick a pair of electrons

    // set up the random number generator
    // initialise the seed engine
    std::random_device rd;
    // use Mersenne-Twister random number generator
    std::mt19937 rng(rd());
    // uniform distribution from 0.0 to 1.0
    std::uniform_real_distribution<double> uniform_dist(0.0,1.0);
    std::vector<int> alpha_num;
    std::vector<int> beta_num;
    int alpha_req=0;
    int beta_req=0;
    std::vector<int> src;

    elecs.clear();
    src.clear();

    // pick the n-th alpha or beta electrons from the determinant
    // according to the availability of pairs and the weighting of
    // opposite-spin pairs relative to same-spin ones
    double rand = uniform_dist(rng);
    if (rand < pParallel){
        // same spin pair
        pgen = pParallel / static_cast<double>(aa_elec_pairs+bb_elec_pairs);
        rand = (rand/pParallel)*static_cast<double>(aa_elec_pairs+bb_elec_pairs);
        int id = std::floor(rand);
        if (id < aa_elec_pairs){
            // alpha alpha spin pair
            alpha_req = 2;
            beta_req = 0;
            // which alpha spin electrons ?
            // triangular indexing system
            alpha_num.push_back(std::ceil((1+std::sqrt(9 + 8*static_cast<double>(id)))/2));
            alpha_num.push_back(id+1-((alpha_num[0]-1)*(alpha_num[0]-2))/2);
        }
        else{
            // beta beta spin pair
            alpha_req = 0;
            beta_req = 2;
            // which beta spin electrons ?
            // triangular indexing system
            id -= aa_elec_pairs;
            beta_num.push_back(std::ceil((1+std::sqrt(9+8*static_cast<double>(id)))/2));
            beta_num.push_back(id+1-((beta_num[0]-1)*(beta_num[0]-2))/2);
        }
    }
    else{
        // opposite spin pair
        alpha_req = 1;
        beta_req = 1;
        pgen = (1.0-pParallel)/static_cast<double>(ab_elec_pairs);
        rand = ((rand - pParallel)/(1.0 - pParallel))*ab_elec_pairs;
        int id = std::floor(rand);
        alpha_num.push_back(1+(id % nalphaels));
        beta_num.push_back(1+std::floor(id/static_cast<double>(nalphaels)));

    }


    // loop through the source determinant and choose the relevant electrons
    int alpha_count = 0;
    int beta_count = 0;
    int elecs_found = 0;

    for (size_t i=0; i<sourceOrbs.size(); ++i){
        if ((sourceOrbs[i]%2)==0){
            // alpha spin
            alpha_count += 1;
            if (alpha_req > 0){
                if (alpha_count == alpha_num[alpha_req-1]){
                    // found an alpha electron
                    elecs_found += 1;
                    elecs.push_back(i);
                    alpha_req -= 1;
                }
            }
        }
        else{
            // beta spin
            beta_count += 1;
            if (beta_req > 0){
                if (beta_count == beta_num[beta_req-1]){
                    // found a beta electron
                    elecs_found += 1;
                    elecs.push_back(i);
                    beta_req -= 1;
                }
            }
        }
        if ((alpha_req==0) and (beta_req==0)){
            break;
        }
    }

    // the spin orbitals
    src.push_back(sourceOrbs[elecs[0]]);
    src.push_back(sourceOrbs[elecs[1]]);

    return src;

}

//---------------------------------------------------------------------------------------------------//

double WeightedExcitgen::pGenSelectDoubleHole(detType const &source, std::vector<int> const &src,
		int const &orb_pair, double &cum_sum, int const &tgt){
    // for a double excitation pick a pair of holes using <e1e1|h1h1> and
    // <e2e2|h2h2> for the probability distribution


    double cpt=0.0;
    double tmp=0.0;
    // selector
    WeightedSelector selector(*H,source);
    // for a parallel spin pair without the partExact() and linExact()
    // schemes the CDFs need not be recomputed

    if ((src[0]%2) == (src[1]%2)){
        // same spin pair

        cum_sum = 0.0;
        tmp = 0.0;
        cpt = 0.0;

        if ((src[0]%2)==1){
            // beta spin pair
            for (int i=1; i<norbs; i+=2){
                if ((not source[i]) and (i!=orb_pair)){
                    // spin orbital is vacant
                    tmp = selector.sameSpinPairContribution(src[0],src[1],i,orb_pair);
                    cum_sum += tmp;
                    if (i==tgt){
                        cpt = tmp;
                    }
                }
            }
        }
        else{
            // alpha spin pair
            for (int i=0; i<norbs; i+=2){
                if ((not source[i]) and (i!=orb_pair)){
                    // spin orbital is vacant
                    tmp = selector.sameSpinPairContribution(
                    		src[0],src[1],i,orb_pair);
                    cum_sum += tmp;
                    if (i==tgt) cpt = tmp;
                }
            }
        }
    }
    else{
        // opposite spin pair
        // only the electron-hole interactions with the same spin is required

        cum_sum = 0.0;
        tmp = 0.0;
        cpt = 0.0;

        if ((tgt%2)==(src[0]%2)){
            if ((tgt%2)==1){
                // beta spin pair
                for (int i=1; i<norbs; i+=2){
                    if ((not source[i]) and (i!=orb_pair)){
                        // spin orbital is vacant
                        tmp = selector.oppSpinPairContribution(src[0],src[1],i,orb_pair);
                        cum_sum += tmp;
                        if (i==tgt){
                            cpt = tmp;
                        }
                    }
                }
            }
            else{
                // alpha spin pair
                for (int i=0; i<norbs; i+=2){
                    if ((not source[i]) and (i!=orb_pair)){
                        // spin orbital is vacant
                        tmp = selector.oppSpinPairContribution(
                        		src[0],src[1],i,orb_pair);
                        cum_sum += tmp;
                        if (i==tgt){
                            cpt = tmp;
                        }
                    }
                }
            }
        }
        else{
            if ((tgt%2)==1){
                // beta spin pair
                for (int i=1; i<norbs; i+=2){
                    if ((not source[i]) and (i!=orb_pair)){
                        // spin orbital is vacant
                        tmp = selector.oppSpinPairContribution(
                        		src[1],src[0],i,orb_pair);
                        cum_sum += tmp;
                        if (i==tgt){
                            cpt = tmp;
                        }
                    }
                }
            }
            else{
                // alpha spin pair
                for (int i=0; i<norbs; i+=2){
                    if ((not source[i]) and (i!=orb_pair)){
                        // spin orbital is vacant
                        tmp = selector.oppSpinPairContribution(
                        		src[1],src[0],i,orb_pair);
                        cum_sum += tmp;
                        if (i==tgt){
                            cpt = tmp;
                        }
                    }
                }
            }

        }

    }


    return cpt;

}

//---------------------------------------------------------------------------------------------------//

double WeightedExcitgen::pGenSingleExcitCs(detType const &source, int srcOrb, int tgt) {
    // evaluate the generation probability for a single excitation

    double cum_sum=0.0;
    double hel=0.0;
    double cpt_tgt=0.0;
    int idi,ida;
    double pgen_single;

    pgen_single = 0.0;

    // the electron involved in the excitation is picked uniformly
    int nel = sourceOrbs.size();
    pgen_single = 1.0 / static_cast<double>(nel);

    idi = 0;
    ida = 0;
    idi = H->getId(srcOrb+1);

    cpt_tgt = 0.0;
    cum_sum = 0.0;

    // construct the cumulative list of connection strengths
    if ((srcOrb%2)==0){
        // alpha spin orbital
        cum_sum = 0.0;
        for (int i=0; i<norbs; i+= 2){
            if (not source[i]){
                // orbital is not in the source determinant
                // get Hamiltonian matrix element
                hel = 0.0;
                ida = H->getId(i+1);
                std::vector<int> idsrc;
                for (size_t j=0; j<sourceOrbs.size(); ++j){
                    idsrc.push_back(H->getId(sourceOrbs[j]+1));
                    if (sourceOrbs[j] == srcOrb){
                        continue;
                    }
                    hel += H->getMatrixElement(idi,idsrc[j],ida,idsrc[j]);
                    if ((srcOrb%2) == (sourceOrbs[j]%2)){
                        hel -= H->getMatrixElement(idi,idsrc[j],idsrc[j],ida);
                    }
                }
                hel += H->getMatrixElement(idi,ida);
                cum_sum += std::fabs(hel);
                if (i == tgt){
                    cpt_tgt = std::fabs(hel);
                }
            }
        }
    }
    else{
        // beta spin orbital
        cum_sum = 0.0;
        for (int i=1; i<norbs; i+=2){
            if (not source[i]){
                // orbital is not in the source determinant
                // get Hamiltonian matrix element
                hel = 0.0;
                ida = H->getId(i+1);
                std::vector<int> idsrc;
                for (size_t j=0; j<sourceOrbs.size(); ++j){
                    idsrc.push_back(H->getId(sourceOrbs[j]+1));
                    if (sourceOrbs[j] == srcOrb){
                        continue;
                    }
                    hel += H->getMatrixElement(idi,idsrc[j],ida,idsrc[j]);
                    if ((srcOrb%2) == (sourceOrbs[j]%2)){
                        hel -= H->getMatrixElement(idi,idsrc[j],idsrc[j],ida);
                    }
                }
                hel += H->getMatrixElement(idi,ida);
                cum_sum += std::fabs(hel);
                if (i == tgt){
                    cpt_tgt = std::fabs(hel);
                }
            }
        }
    }

    // the generation probability
    if (std::fabs(cum_sum) < 1e-12){
        pgen_single = 0.0;
    }
    else{
        pgen_single *= (cpt_tgt/cum_sum);
    }

    return pgen_single;
}

//---------------------------------------------------------------------------------------------------//

double WeightedExcitgen::getExcitationProb(
		detType const &source, detType const &target){
    // evaluate the generation probability biased according to
    // the strength of the connecting Hamiltonian matrix elements

    double diff = 0.0;

    std::vector<int> src,tgt;
    //src.assign(2,-1);
    //tgt.assign(2,-1);

    // get the source and target orbitals of the excitation
    for (size_t i=0; i<source.size(); ++i){
        diff = static_cast<int>(source[i]) - static_cast<int>(target[i]);
        if (diff > 0){
            // holes in source
            src.push_back(i);
        }
        if (diff < 0){
            // particles in tgt
            tgt.push_back(i);
        }
    }
    if (src.size() < 2){
        src.push_back(-1);
    }
    if (tgt.size() < 2){
        tgt.push_back(-1);
    }
    if (src.size() != tgt.size()){
        throw SizeMismatchError(src.size(),tgt.size());
    }
// TODO WARNING the ordering of tgt is not necessarily the same as in
// the generation. In some cases, the probability depends on the order
// of tgt, then, this will return another probability
   return calcPgen(source,src,tgt);
}

//---------------------------------------------------------------------------------------------------//

double WeightedExcitgen::calcPgen(detType const &source, std::vector<int> const &src,
		std::vector<int> const &tgt) {
    // calculate the generation probability

    std::vector<double> cpt_pair,int_cpt,sum_pair,cum_sum;
    std::vector<int> tmptgt = tgt;
    cpt_pair.assign(2,0.0);
    int_cpt.assign(2,0.0);
    cum_sum.assign(2,0.0);
    sum_pair.assign(2,0.0);
    double calcedP{0.0};

    if ((src[0] >= 0) and (src[1] < 0)){
        // single excitation
        calcedP = (1-pDoubles)* pGenSingleExcitCs(source,src[0],tmptgt[0]);
    }
    else if ((src[0] >= 0) and (src[1] >= 0)){
        // double eccitation
    	calcedP = pDoubles;

        // biasing of selection of electron pair
        if ((src[0]%2)==(src[1]%2)){
        	calcedP *= pParallel / static_cast<double>(aa_elec_pairs + \
                    bb_elec_pairs);
        }
        else{
        	calcedP *= (1.0-pParallel) / static_cast<double>(ab_elec_pairs);
        }

        // ensure that the target orbitals match the source orbitals
        if (((src[0]%2)!=(tgt[0]%2)) or ((src[1]%2)!=(tgt[1]%2))){
            tmptgt[0] = tgt[1];
            tmptgt[1] = tgt[0];
        }

        int_cpt[0] = pGenSelectDoubleHole(source,src,-1,cum_sum[0],tmptgt[0]);
        int_cpt[1] = pGenSelectDoubleHole(source,src,tmptgt[0],cum_sum[1],tmptgt[1]);

        // deal with the cases when there are no available excitation
        if (std::any_of(cum_sum.begin(),cum_sum.end(), [](const double i) \
                    {return std::fabs(i)<1e-12;})){
            cum_sum[0] = 0.0;
            cum_sum[1] = 0.0;
            int_cpt[0] = 0.0;
            int_cpt[1] = 0.0;
        }

        // adjust the probabilities: if the two holes are selected from the
        // same list, they could have been selected either way round
        cpt_pair.assign(2,0.0);
        sum_pair.assign(2,0.0);
        if ((tmptgt[0]%2)==(tmptgt[1]%2)){
            if (H->partExact() or H->linExact()){
            	cpt_pair[0] = pGenSelectDoubleHole(source,src,-1,sum_pair[0],tmptgt[1]);
                cpt_pair[1] = pGenSelectDoubleHole(source,src,tmptgt[1],sum_pair[1],tmptgt[0]);
            }
            else{
                cpt_pair[0] = int_cpt[1];
                cpt_pair[1] = int_cpt[0];
                sum_pair[0] = cum_sum[0];
                sum_pair[1] = cum_sum[0] - int_cpt[1];
            }
            // deal with the cases when there are no available excitation
            if (std::any_of(sum_pair.begin(),sum_pair.end(), [](const double i) \
                        {return std::fabs(i)<1e-12;})){
                sum_pair[0] = 0.0;
                sum_pair[1] = 0.0;
                cpt_pair[0] = 0.0;
                cpt_pair[1] = 0.0;
            }
            double prod_int_cpt = std::accumulate(int_cpt.begin(),int_cpt.end(),1,std::multiplies<double>());
            double prod_cum_sum = std::accumulate(cum_sum.begin(),cum_sum.end(),1,std::multiplies<double>());
            double prod_cpt_pair = std::accumulate(cpt_pair.begin(),cpt_pair.end(),1,std::multiplies<double>());
            double prod_sum_pair = std::accumulate(sum_pair.begin(),sum_pair.end(),1,std::multiplies<double>());
            // it seems the std::accumulate does not work
            // I use explicit reduction until I can check this thoroughly
            prod_int_cpt = int_cpt[0]*int_cpt[1];
            prod_cpt_pair = cpt_pair[0]*cpt_pair[1];
            prod_cum_sum = cum_sum[0]*cum_sum[1];
            prod_sum_pair = sum_pair[0]*sum_pair[1];
            calcedP *= (prod_int_cpt/prod_cum_sum + prod_cpt_pair/prod_sum_pair);
        }
        else{
            double prod_int_cpt = std::accumulate(int_cpt.begin(),int_cpt.end(),1,std::multiplies<double>());
            double prod_cum_sum = std::accumulate(cum_sum.begin(),cum_sum.end(),1,std::multiplies<double>());
            // same as above
            prod_cum_sum = cum_sum[0]*cum_sum[1];
            prod_int_cpt = int_cpt[0]*int_cpt[1];
            calcedP *= (prod_int_cpt/prod_cum_sum);
        }

    }
    else{
        // neither single nor double excitation
    	calcedP = 0.0;
    }
    return calcedP;
}

} /* namespace networkVMC */
