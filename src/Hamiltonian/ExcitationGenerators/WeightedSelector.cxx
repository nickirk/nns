/*
 * WeightedSelector.cxx
 *
 *  Created on: Jun 19, 2018
 *      Author: guther
 */

#include "WeightedSelector.hpp"
#include "../Hamiltonian.hpp"
#include "../../HilbertSpace/Determinant.hpp"
#include <vector>
#include <cmath>
#include <iterator>

namespace networkVMC{

WeightedSelector::~WeightedSelector() {
}

int WeightedSelector::selectSingleHole(int src, double &pgen) {
    // select a hole for a single excitation

    // set up the random number generator
    // initialise the seed engine
    std::random_device rd;
    // use Mersenne-Twister random number generator
    std::mt19937 rng(rd());
    // uniform distribution from 0.0 to 1.0
    std::uniform_real_distribution<double> uniform_dist(0.0,1.0);
    std::vector<double> cdf_single_a;
    std::vector<int> holes_a;
    double cum_sum=0.0;
    double hel=0.0;
    int nexcit=0;
    int idi,ida;
    int tgt = -1;

    // orbitals that are occupied in source
    auto sourceOrbs = getOccupiedPositions(source);
    // number of electrons as given by number of occupied orbitals
    int nel = sourceOrbs.size();

    pgen = 0.0;

    idi = 0;
    ida = 0;
    idi = H.getId(src+1);

    // maximum number of available holes if no symmetry and
    // spin information are considered
    int nelem = norbs-nel;

    // construct the cumulative list of connection strengths
    if ((src%2)==0){
        // alpha spin orbital
        cum_sum = 0.0;
        for (size_t i=0; i<norbs; i+= 2){
            if (not source[i]){
                // orbital is not in the source determinant
                // get Hamiltonian matrix element
                nexcit += 1;
                hel = 0.0;
                ida = H.getId(i+1);
                std::vector<int> idsrc;
                for (size_t j=0; j<sourceOrbs.size(); ++j){
                    idsrc.push_back(H.getId(sourceOrbs[j]+1));
                    if (sourceOrbs[j] == src){
                        continue;
                    }
                    hel += H.getMatrixElement(idi,idsrc[j],ida,idsrc[j]);
                    if ((src%2) == (sourceOrbs[j]%2)){
                        hel -= H.getMatrixElement(idi,idsrc[j],idsrc[j],ida);
                    }
                }
                hel += H.getMatrixElement(idi,ida);
                cum_sum += std::fabs(hel);
                holes_a.push_back(i);
                cdf_single_a.push_back(cum_sum);
            }
        }
    }
    else{
        // beta spin orbital
        cum_sum = 0.0;
        for (size_t i=1; i<norbs; i+=2){
            if (not source[i]){
                // orbital is not in the source determinant
                // get Hamiltonian matrix element
                hel = 0.0;
                nexcit += 1;
                ida = H.getId(i+1);
                std::vector<int> idsrc;
                for (size_t j=0; j<sourceOrbs.size(); ++j){
                    idsrc.push_back(H.getId(sourceOrbs[j]+1));
                    if (sourceOrbs[j] == src){
                        continue;
                    }
                    hel += H.getMatrixElement(idi,idsrc[j],ida,idsrc[j]);
                    if ((src%2) == (sourceOrbs[j]%2)){
                        hel -= H.getMatrixElement(idi,idsrc[j],idsrc[j],ida);
                    }
                }
                hel += H.getMatrixElement(idi,ida);
                cum_sum += std::fabs(hel);
                holes_a.push_back(i);
                cdf_single_a.push_back(cum_sum);
            }
        }
    }

    // pick a particular hole a for use
    if ((nexcit == 0) or (std::fabs(cum_sum) < 1e-12)){
        // no available excitation
        tgt = -1;
        pgen = 0.0;

        return tgt;
    }
    else if (nexcit > nelem){
        throw OutOfRangeError(nexcit);
    }
    else{
        double hel_picked = uniform_dist(rng)*cum_sum;
        std::vector<double>::iterator up;
        up = std::upper_bound(cdf_single_a.begin(),cdf_single_a.end(),hel_picked);
        int orbind = up - cdf_single_a.begin();
        tgt = holes_a[orbind];

        // the generation probability
        if (orbind > 0){
            pgen = (cdf_single_a[orbind]-cdf_single_a[orbind-1])/cum_sum;
        }
        else{
            pgen = cdf_single_a[orbind]/cum_sum;
        }

        return tgt;
    }
}


//---------------------------------------------------------------------------------------------------//

int WeightedSelector::selectDoubleHole(std::vector<int> const &src,
		std::size_t orb_pair, double &cum_sum, double &cpt){
    // for a double excitation pick a pair of holes using <e1e1|h1h1> and
    // <e2e2|h2h2> for the probability distribution


    // set up the random number generator
    // initialise the seed engine
    std::random_device rd;
    // use Mersenne-Twister random number generator
    std::mt19937 rng(rd());
    // uniform distribution from 0.0 to 1.0
    std::uniform_real_distribution<double> uniform_dist(0.0,1.0);
    std::vector<int> alpha_num;
    int nexcit;
    std::vector<double> cdf;
    std::vector<int> holes;
    int tgt=0;

    // for a parallel spin pair without the partExact() and linExact
    // schemes the CDFs need not be recomputed

    if ((src[0]%2) == (src[1]%2)){
        // same spin pair

        cum_sum = 0.0;
        nexcit = 0;

        if ((src[0]%2)==1){
            // beta spin pair
            for (size_t i=1; i<norbs; i+=2){
                if ((not source[i]) and (i!=orb_pair)){
                    // spin orbital is vacant
                    cum_sum += this->sameSpinPairContribution(src[0],src[1],i,orb_pair);
                    nexcit += 1;
                    cdf.push_back(cum_sum);
                    holes.push_back(i);
                }
            }
        }
        else{
            // alpha spin pair
            for (size_t i=0; i<norbs; i+=2){
                if ((not source[i]) and (i!=orb_pair)){
                    // spin orbital is vacant
                    cum_sum += this->sameSpinPairContribution(src[0],src[1],i,orb_pair);
                    nexcit += 1;
                    cdf.push_back(cum_sum);
                    holes.push_back(i);
                }
            }
        }
    }
    else{
        // opposite spin pair
        // only the electron-hole interactions with the same spin is required

        cum_sum = 0.0;
        nexcit = 0;

        if (orb_pair < 0){
            if ((src[0]%2)==1){
                // beta spin pair
                for (size_t i=1; i<norbs; i+=2){
                    if ((not source[i]) and (i!=orb_pair)){
                        // spin orbital is vacant
                        cum_sum += this->oppSpinPairContribution(src[0],src[1],i,orb_pair);
                        nexcit += 1;
                        cdf.push_back(cum_sum);
                        holes.push_back(i);
                    }
                }
            }
            else{
                // alpha spin pair
                for (size_t i=0; i<norbs; i+=2){
                    if ((not source[i]) and (i!=orb_pair)){
                        // spin orbital is vacant
                        cum_sum += this->oppSpinPairContribution(src[0],src[1],i,orb_pair);
                        nexcit += 1;
                        cdf.push_back(cum_sum);
                        holes.push_back(i);
                    }
                }
            }
        }
        else{
            if ((orb_pair%2)==0){
                // beta spin since first orital is alpha spin
                for (size_t i=1; i<norbs; i+=2){
                    if ((not source[i]) and (i!=orb_pair)){
                        // spin orbital is vacant
                        cum_sum += this->oppSpinPairContribution(src[0],src[1],i,orb_pair);
                        nexcit += 1;
                        cdf.push_back(cum_sum);
                        holes.push_back(i);
                    }
                }
            }
            else{
                // alpha spin since first orbital is beta spin
                for (size_t i=0; i<norbs; i+=2){
                    if ((not source[i]) and (i!=orb_pair)){
                        // spin orbital is vacant
                        cum_sum += oppSpinPairContribution(src[0],src[1],i,orb_pair);
                        nexcit += 1;
                        cdf.push_back(cum_sum);
                        holes.push_back(i);
                    }
                }
            }

        }

    }

    // if there are no available hole pairs
    if ((nexcit==0) or (std::fabs(cum_sum)<1e-12)){
        tgt = -1;

        return tgt;
    }

    // choose a value
    double hel_picked = uniform_dist(rng)*cum_sum;
    std::vector<double>::iterator up;
    up = std::upper_bound(cdf.begin(),cdf.end(),hel_picked);
    int orbind = up - cdf.begin();
    tgt = holes[orbind];

    // the generation probability
    if (orbind > 0){
        cpt = cdf[orbind]-cdf[orbind-1];
    }
    else{
        cpt = cdf[orbind];
    }


    return tgt;

}

//---------------------------------------------------------------------------------------------------//

double WeightedSelector::sameSpinPairContribution(int i, int j, int a, int b){
    // the contribution for a pair of same spin electrons
    int idi = H.getId(i+1);
    int idj = H.getId(j+1);
    int ida = H.getId(a+1);
    int idb = H.getId(b+1);
    double contrib = 0.0;

    if (H.partExact() and (b > 0)){
        // include sqrt(abs(<ij|ab>-<ij|ba>))
        contrib = std::max(std::sqrt(std::fabs(H.getMatrixElement(idi,idj,ida,idb)-H.getMatrixElement(idi,idj,idb,ida))),1e-5);
    }
    else if (H.linExact()){
        if (b > 0){
            // include abs(<ij|ab>-<ij|ba>)
            contrib = std::fabs(H.getMatrixElement(idi,idj,ida,idb)-H.getMatrixElement(idi,idj,idb,ida));
        }
        else{
            // select the first orbital linearly
            contrib = 1.0;
        }
    }
    else{
        // include the contribution of this terms sqrt(<ii|aa>+<jj|aa>)
        contrib = std::sqrt(std::fabs(H.getMatrixElement(idi,idi,ida,ida)+H.getMatrixElement(idj,idj,ida,ida)));
    }

    return contrib;

}

//---------------------------------------------------------------------------------------------------//

double WeightedSelector::oppSpinPairContribution(int i, int j, int a, int b){
    // the contribution for a pair of opposite spin electrons

    int idi = H.getId(i+1);
    int idj = H.getId(j+1);
    int ida = H.getId(a+1);
    int idb = H.getId(b+1);
    double contrib = 0.0;

    if (H.partExact() and (b > 0)){
        // include sqrt(abs(<ij|ab>))
        // this can only be used if <ij|ba> = 0
        contrib = std::max(std::sqrt(std::fabs(H.getMatrixElement(idi,idj,ida,idb))),1e-5);
    }
    else if (H.linExact()){
        if (b > 0){
            // include abs(<ij|ab>)
            // this can only be used if <ij|ba> = 0
            contrib = std::fabs(H.getMatrixElement(idi,idj,ida,idb));
        }
        else{
            // select the first orbital linearly
            contrib = 1.0;
        }
    }
    else{
        // include the contribution of this terms sqrt(<ii|aa>)
        contrib = std::sqrt(std::fabs(H.getMatrixElement(idi,idi,ida,ida)));
    }

    return contrib;

}

//---------------------------------------------------------------------------------------------------//

std::vector<int> WeightedSelector::selectDoubleHoles(std::vector<int> const &src, std::vector<double> &cum_sum, std::vector<double> &cpt){
    // for a double excitation pick a pair of holes using <e1e1|h1h1> and
    // <e2e2|h2h2> for the probability distribution


    // set up the random number generator
    // initialise the seed engine
    std::random_device rd;
    // use Mersenne-Twister random number generator
    std::mt19937 rng(rd());
    // uniform distribution from 0.0 to 1.0
    std::uniform_real_distribution<double> uniform_dist(0.0,1.0);
    std::vector<int> alpha_num;
    int nexcita,nexcitb;
    std::vector<double> cdf_a,cdf_b;
    std::vector<int> holes_a,holes_b;
    std::vector<int> tgt;

    // for a parallel spin pair without the partExact() and linExact()
    // schemes the CDFs need not be recomputed

    if ((src[0]%2) == (src[1]%2)){
        // same spin pair

        cum_sum.clear();
        cum_sum.assign(2,0.0);
        //cum_sum.push_back(0.0);
        //cum_sum.push_back(0.0);
        nexcita = 0;
        nexcitb = 0;

        if ((src[0]%2)==1){
            // beta spin pair
            for (size_t i=1; i<norbs; i+=2){
                if (not source[i]){
                    // spin orbital is vacant
                    cum_sum[0] += sameSpinPairContribution(src[0],src[1],i,-1);
                        nexcita += 1;
                        cdf_a.push_back(cum_sum[0]);
                        holes_a.push_back(i);
                }
            }
        }
        else{
            // alpha spin pair
            for (size_t i=0; i<norbs; i+=2){
                if (not source[i]){
                    // spin orbital is vacant
                    cum_sum[0] += sameSpinPairContribution(src[0],src[1],i,-1);
                        nexcita += 1;
                        cdf_a.push_back(cum_sum[0]);
                        holes_a.push_back(i);
                }
            }
        }

        // if there are no available hole pairs
        if ((nexcita==0) or (std::fabs(cum_sum[0])<1e-12)){
            tgt.clear();
            tgt.assign(2,-1);
            //tgt.push_back(-1);
            //tgt.push_back(-1);

            return tgt;
        }

        // choose a value to pick orbital a
        double hel_picked_a = uniform_dist(rng)*cum_sum[0];
        std::vector<double>::iterator up;
        up = std::upper_bound(cdf_a.begin(),cdf_a.end(),hel_picked_a);
        int orbind_a = up - cdf_a.begin();
        tgt.push_back(holes_a[orbind_a]);

        // the generation probability
        if (orbind_a > 0){
            cpt.push_back(cdf_a[orbind_a]-cdf_a[orbind_a-1]);
        }
        else{
            cpt.push_back(cdf_a[orbind_a]);
        }

        // the CDF for the orbital b is the same for b<a
        //vector<double> cdf_b(cdf_a.begin(),(cdf_a.begin()+orbind_a+1-1));
        //vector<int> holes_b(holes_a.begin(),(holes_a.begin()+orbind_a+1-1));
        std::vector<double> cdf_b(cdf_a.begin(),(cdf_a.begin()+orbind_a));
        std::vector<int> holes_b(holes_a.begin(),(holes_a.begin()+orbind_a));
        // for b>=a it is the CDF with the contribution of a removed
        nexcitb = nexcita -1;
        // if there are no available hole pairs
        if ((nexcitb==0) or (std::fabs(cum_sum[0]-cpt[0])<1e-12)){
            tgt.clear();
            tgt.assign(2,-1);
            //tgt.push_back(-1);
            //tgt.push_back(-1);

            return tgt;
        }

        for (size_t i=orbind_a; i<(cdf_a.size()-1); ++i){
            cdf_b.push_back(cdf_a[i+1]-cpt[0]);
            holes_b.push_back(holes_a[i+1]);
        }

        cum_sum[1] = cum_sum[0] - cpt[0];

        // choose a value to pick orbital b
        double hel_picked_b = uniform_dist(rng)*cum_sum[1];
        std::vector<double>::iterator up2;
        up2 = std::upper_bound(cdf_b.begin(),cdf_b.end(),hel_picked_b);
        int orbind_b = up2 - cdf_b.begin();
        tgt.push_back(holes_b[orbind_b]);

        // the generation probability
        if (orbind_b > 0){
            cpt.push_back(cdf_b[orbind_b]-cdf_b[orbind_b-1]);
        }
        else{
            cpt.push_back(cdf_b[orbind_b]);
        }


        return tgt;

    }
    else{
        tgt.clear();
        tgt.assign(2,-1);
        //tgt.push_back(-1);
        //tgt.push_back(-1);

        return tgt;
    }

}

}

