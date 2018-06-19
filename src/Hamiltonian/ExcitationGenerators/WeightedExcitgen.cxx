/*
 * WeightedExcitgen.cxx
 *
 *  Created on: Jun 15, 2018
 *      Author: guther
 */

#include "WeightedExcitgen.hpp"
#include "../../HilbertSpace/Determinant.hpp"

namespace networkVMC {

WeightedExcitgen::WeightedExcitgen() {
}

//---------------------------------------------------------------------------------------------------//

WeightedExcitgen::~WeightedExcitgen() {
}

//---------------------------------------------------------------------------------------------------//

detType WeightedExcitgen::generateDoubleExcit() {
    // generate a double excitation

    std::vector<double> cpt_pair,int_cpt,sum_pair,cum_sum;
    std::vector<int> elecs;
    std::vector<int> tgt;
    tgt.assign(2,-1);

    // select a pair of electrons
    std::vector<int> src = pickBiasedElecs(elecs);


    // select the two target holes by approximate connection
    // strength
    if ((not H.part_exact) and (not H.lin_exact) and ((src[0]%2)==(src[1]%2))){
        // in this case the CDFs are the same
        std::vector<int> tgt = this->select_double_holes(H,src,cum_sum,int_cpt);
    }
    else{
        //std::vector<int> tgt;
        //tgt.assign(2,-1);
        int_cpt.assign(2,0.0);
        cum_sum.assign(2,0.0);

        tgt[0] = this->select_double_hole(H,src,-1,cum_sum[0],int_cpt[0]);

        if (tgt[0] > -1){
            tgt[1] = this->select_double_hole(H,src,tgt[0],cum_sum[1],int_cpt[1]);
        }
    }

    // check that the two electrons have not been excited into the
    // same orbitals
    if (std::any_of(tgt.begin(),tgt.end(), [](const int i) {return i<0;})){
        pgen = 0.0;
        excitmat[0][0] = -1;
        excitmat[0][1] = -1;
        excitmat[1][0] = -1;
        excitmat[1][1] = -1;
        target.clear();
    }

    // adjust the probabilities: if the two holes are selected from
    // the same list, they could have been selected in any order
    cpt_pair.assign(2,0.0);
    sum_pair.assign(2,0.0);
    if ((tgt[0]%2)==(tgt[1]%2)){
        if (H.part_exact or H.lin_exact){
            cpt_pair[0] = this->pgen_select_double_hole(H,src,-1,sum_pair[0],tgt[1]);
            cpt_pair[1] = this->pgen_select_double_hole(H,src,tgt[1],sum_pair[1],tgt[0]);
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
    target = source;
    annihilate(target,src[0]);
    annihilate(target,src[1]);
    create(target,tgt[0]);
    create(target,tgt[1]);
    excitmat[0][0] = src[0];
    excitmat[1][0] = src[0];
    if ((src[0]%2)==(tgt[0]%2)){
        excitmat[0][1] = tgt[0];
        excitmat[1][1] = tgt[1];
    }
    else{
        excitmat[0][1] = tgt[1];
        excitmat[1][1] = tgt[0];
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
    auto sourceOrbs = getOccupiedPositions(source);

    elecs.clear();
    src.clear();

    // pick the n-th alpha or beta electrons from the determinant
    // according to the availability of pairs and the weighting of
    // opposite-spin pairs relative to same-spin ones
    double rand = uniform_dist(rng);
    if (rand < pParallel){
        // same spin pair
        pGen = pParallel / static_cast<double>(aa_elec_pairs+bb_elec_pairs);
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
        pGen = (1.0-pParallel)/static_cast<double>(ab_elec_pairs);
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


} /* namespace networkVMC */
