/*
 * UniformExcitgen.cxx
 *
 *  Created on: Jun 15, 2018
 *      Author: guther
 */

#include "UniformExcitgen.hpp"
#include <vector>
#include "../HilbertSpace/Determinant.cxx"

namespace networkVMC {

UniformExcitgen::UniformExcitgen(detType const &HF){
	setProbabilities(HF);
}

//---------------------------------------------------------------------------------------------------//

UniformExcitgen::~UniformExcitgen() {
}

//---------------------------------------------------------------------------------------------------//

detType UniformExcitgen::generateExcitation(detType const &source, double &pGet) const{

    // random excitation generator: randomly pick a connected determinant

    int exflag;
    //double pdoubles;
    detType target;
    std::vector<int> holes,particles;
    double rand;
    std::random_device rng;
    int nsingleexcit,ndoubleexcit,nexcit;
    double const normalisation=static_cast<double>(rng.max()); // minimal potentially returned values is 0
    std::vector<int> noccs,nunoccs;

    exflag = 3;
    //pdoubles = 0.0;
    target = source;


    // get the occupied spin orbitals
    std::vector<int> source_orbs = getOccupiedPositions(source);

    // number of occupied and unoccupied alpha and beta spin orbitals
    // noccs[0] : number of occupied alpha spin orbitals
    // noccs[1] : number of occupied beta spin orbitals
    // nunoccs[0] : number of unoccupied alpha spin orbitals
    // nunoccs[1] : number of unoccupied beta spin orbitals
    noccs.push_back(0);
    noccs.push_back(0);
    nunoccs.push_back(0);
    nunoccs.push_back(0);

    for (size_t i=0; i<source_orbs.size(); ++i){
        //if (source[i]){
        if ((source_orbs[i]%2)==0){
            // alpha spin
            noccs[0] += 1;
        }
        else{
            // beta spin
            noccs[1] += 1;
        }
        //}
    }
    for (size_t i=0; i<noccs.size(); ++i){
        nunoccs[i] = (source.size()/2) - noccs[i];
    }

    // for generation probabilities
    // number of excitations
    nsingleexcit = 0;
    ndoubleexcit = 0;
    if ((exflag==1)||(exflag==3)){

        // single excitations
        nsingleexcit = (noccs[0]*nunoccs[0]) + (noccs[1]*nunoccs[1]);

    }
    if ((exflag==2)||(exflag==3)){
        // double excitations

        // alpha-alpha, beta-beta, alpha-beta pairs
        ndoubleexcit = ((noccs[0]*(noccs[0]-1)/2)*(nunoccs[0]*(nunoccs[0]-1)/2)) +
            ((noccs[1]*(noccs[1]-1)/2)*(nunoccs[1]*(nunoccs[1]-1)/2)) +
            (noccs[0]*noccs[1]*nunoccs[0]*nunoccs[1]);
    }
    nexcit = nsingleexcit + ndoubleexcit;

    // use relative number of single and double excitations for probability
    //pdoubles = static_cast<double>(ndoubleexcit)/static_cast<double>(nexcit);

    // exflag = 3: single or double excitation
    // exflag = 1: single excitations only
    // exflag = 2: double excitations only
    int ic = exflag;
    double pdoubnew=0.0;

    if (ic==1){
        pdoubnew = 0.0;
    }
    else if (ic==2){
        pdoubnew = 1.0;
    }
    else if (ic==3){
        rand = static_cast<double>(rng())/normalisation;
        if (rand<pDoubles){
            ic = 2;
        }
        else{
            ic = 1;
        }
        pdoubnew = pDoubles;
    }

    // call the respective single or double excitation generators
    if (ic==2){
        // double excitation
        target = this->genDoubleExcitation(source, source_orbs, holes, particles, pGet, pdoubnew, noccs, nunoccs);
    }
    else{
        // single excitation
        target = this->genSingleExcitation(source, source_orbs, holes, particles, pGet, pdoubnew, noccs, nunoccs);
    }


    return target;

}

//---------------------------------------------------------------------------------------------------//

detType UniformExcitgen::genSingleExcitation(detType const &source, std::vector<int> const &source_orbs,
		std::vector<int> &holes, std::vector<int> &particles, double &pgen, double pdoubnew,
		std::vector<int> noccs, std::vector<int> nunoccs) const{

    // randomly pick a single excitation (random excitation generator)

    int attempts,npairs,elecwnoexcit,nexcit;
    detType target;
    double rand;
    int iorb,aorb,ielec,aelec,ispin;
    int nel,norbs;
    std::random_device rng;
    double const normalisation=static_cast<double>(rng.max()); // minimal potentially returned values is 0

    // number of electrons and spin orbitals
    nel = source_orbs.size();
    norbs = source.size();

    iorb = 0;
    aorb = 0;
    ielec = 0;
    aelec = 0;
    ispin = 0;

    target = source;


    // find the number of available pairs of particles and holes for each spin type
    npairs = 0;
    npairs = (noccs[0]*nunoccs[0]) + (noccs[1]*nunoccs[1]);
    // find the number of electrons with no available
    // excitations
    elecwnoexcit = 0;
    if (nunoccs[0]==0){
        elecwnoexcit += noccs[0];
    }
    if (nunoccs[1]==0){
        elecwnoexcit += noccs[1];
    }

    // if there are no single excitations return a null excitations
    if (npairs==0){
        particles.clear();
        holes.clear();
        target.clear();
        pgen = 0.0;

        return target;
    }

    attempts = 0;

    while (true){

        // choose a random electron
        rand = rng()/normalisation;
        // conversion into an integer in the range [a,b]
        // int(a+rand*(b-a+1))
        ielec = static_cast<int>(static_cast<double>(nel*rand));
        iorb = source_orbs[ielec];

        nexcit = 0;
        // find the spin of the electron
        // and the number of possible excitations for it
        if ((iorb%2)==0){
            // alpha spin
            ispin = 0;
            nexcit = nunoccs[0];
        }
        else{
            // beta spin
            ispin = 1;
            nexcit = nunoccs[1];
        }

        // there are allowed excitations
        if (nexcit!=0){
            break;
        }

        // to avoid infinity loops
        if (attempts>250){
            throw OutOfRangeError(attempts);
        }

        attempts += 1;

    }

    // need to pick an unoccupied orbital a
    // there are tqo poaaibilities in order to achieve this
    // possibility 1:
    // cycle through all orbitals of the same spin
    // only counting the unoccupied ones to find the
    // correct target determinant
    // possibility 2:
    // pick orbitals at random and redraw until an allowed one is chosen

    if (false){
        // run through all orbitals of the same spin
        // to find the allowed one out of nexcit

        // choose the unoccupied orbital a
        rand = rng()/normalisation;
        aelec = static_cast<int>(static_cast<double>(nexcit)*rand) + 1;

        // run through all allowed orbitals until this one is found
        int counter = 0;
        for (int i=ispin; i<norbs; ++i){
            if (!source[i]){
                // found an unoccupied orbitals of the correct spin
                counter += 1;
                if (counter==aelec){
                    // found orbital a
                    aorb = i;
                    break;
                }
            }
        }
    }
    else{
        // keep drawing orbitals from the desired spin
        // until an unoccupied one has been found
        attempts = 0;

        while (true){
            // draw randomly from the set of orbitals
            // impose the spin restriction by only
            // drawing either even or odd numbered orbitals
            rand = rng()/normalisation;
            aorb = 2*(1+static_cast<int>(rand*static_cast<double>(norbs/2))) + (ispin-2);

            // test whether it is in the source determinant
            if (!source[aorb]){
                // orbital a is unoccupied
                break;
            }

            // to avoid infinity loops
            if (attempts>250){
                throw OutOfRangeError(attempts);
            }

            attempts += 1;

        }
    }

    // generate the new determinant
    target = source;
    annihilate(target,iorb);
    create(target,aorb);
    holes.clear();
    particles.clear();
    holes.push_back(iorb);
    particles.push_back(aorb);

    // calculate the generation probability
    pgen = (1.0 - pdoubnew) / (static_cast<double>(nexcit*(nel-elecwnoexcit)));

    return target;

}

//---------------------------------------------------------------------------------------------------//

detType UniformExcitgen::genDoubleExcitation(detType const &source, std::vector<int> const &source_orbs,
		std::vector<int> &holes, std::vector<int> &particles, double &pgen, double pdoubnew,
		std::vector<int> noccs, std::vector<int> nunoccs) const{

    // randomly generate a double excitation (random excitation generator)

    detType target=source;
    int nexcita,nexcitb,nexcitotherway;
    bool baorbfail;
    int elecpairs,nel,norbs;
    int aorb,aspin,borb,bspin;
    int nforbiddenorbs;
    std::vector<int> elecs,spin;

    nexcita = 0;
    nexcitb = 0;
    nexcitotherway = 0;
    baorbfail = true;

    // number of electrons and spin orbitals
    nel = source_orbs.size();
    norbs = source.size();

    // number of electron pairs
    elecpairs = (nel*(nel-1))/2;

    // pick an unbiased and distinct electron pair
    elecs = pickElecPair(source_orbs, spin, elecpairs);

    // the number of forbidden orbitals which are not allowed by spin
    nforbiddenorbs = countForbiddenOrbs(spin, nunoccs);

    // pick orbital a
    aorb = pickOrbA(source, spin, noccs, nunoccs, norbs, nel, nforbiddenorbs, nexcita, aspin, baorbfail);

    if (baorbfail){
        // return null excitation
        aorb = -1;
        particles.clear();
        holes.clear();
        target.clear();
        pgen = 0.0;

        return target;
    }

    // pick orbital b
    borb = pickOrbB(source, spin, noccs, nunoccs, norbs, nel, nforbiddenorbs, aorb, aspin, nexcitb, bspin, nexcitotherway);

    // evaluate the probability
    // p(ij->ab) = p_doub x p(i,j) x [p(a|ij)p(b|a,ij) + p(b|ij)p(a|b,ij)]
    pgen = pdoubnew*((1.0/static_cast<double>(nexcitb)) + (1.0/static_cast<double>(nexcitotherway))) /
        (static_cast<double>((elecpairs*(nexcita-nforbiddenorbs))));

    // make the double excitation
    // generate the new determinant
    target = source;
    holes.clear();
    particles.clear();
    for (size_t i=0; i<elecs.size(); ++i){
        annihilate(target,source_orbs[elecs[i]]);
        holes.push_back(source_orbs[elecs[i]]);
    }
    create(target,aorb);
    create(target,borb);
    particles.push_back(aorb);
    particles.push_back(borb);


    return target;

}

//---------------------------------------------------------------------------------------------------//

void UniformExcitgen::setProbabilities(detType example_det){
    // set appropriate initial values for the probabilities based
    // on an example determinant

    int nel,nalpha,nbeta,norbs;

    // number of spin orbitals
    norbs = example_det.size();

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
    pParallel = static_cast<double>(parallel_pairs)/static_cast<double>(parallel_pairs+opp_pairs);

    // number of single excitations
    int nsingleexcit = (nalpha*((norbs/2)-nalpha)) + (nbeta*((norbs/2)-nbeta));

    // number of double excitations
    // aa-pairs + bb-pairs + ab-pairs
    int ndoubleexcit = ((nalpha*(nalpha-1)/2)*(((norbs/2)-nalpha)*((norbs/2)-nalpha-1)/2)) +\
                       ((nbeta*(nbeta-1)/2)*(((norbs/2)-nbeta)*((norbs/2)-nbeta-1)/2)) +\
                       (nalpha*nbeta*((norbs/2)-nalpha)*((norbs/2)-nbeta));

    // p_singles + p_doubles = 1
    double pSingles = static_cast<double>(nsingleexcit)/static_cast<double>(nsingleexcit+ndoubleexcit);
    pDoubles = 1.0 - pSingles;

    std::cout << "\n Setting the probabilities to their initial values: \n" << std::endl;
    std::cout << "p_singles: " << pSingles << std::endl;
    std::cout << "p_doubles: " << pDoubles << std::endl;
    std::cout << "p_parallel: " << pParallel << std::endl;

    return;
}

} /* namespace networkVMC */
