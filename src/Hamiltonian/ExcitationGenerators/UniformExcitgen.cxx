/*
 * UniformExcitgen.cxx
 *
 *  Created on: Jun 15, 2018
 *      Author: Lauretta Schwarz, guther
 */

#include "UniformExcitgen.hpp"
#include <vector>
#include "../../HilbertSpace/Determinant.hpp"
#include "ExcitmatType.hpp"

namespace networkVMC {

UniformExcitgen::UniformExcitgen(detType const &HF):
	clonableExcitgen<UniformExcitgen>(),pBiasGen(ProbUpdater(HF)){
	pParallel = pBiasGen.pParallel();
	pDoubles = pBiasGen.pDoubles();
}

//---------------------------------------------------------------------------------------------------//

UniformExcitgen::~UniformExcitgen() {
}

//---------------------------------------------------------------------------------------------------//

// call the ProbUpdater to get new values for pDoubles/pParallel
void UniformExcitgen::updateBiases(){
	setNewBiases(pBiasGen,pDoubles,pParallel);
}

//---------------------------------------------------------------------------------------------------//

detType UniformExcitgen::generateExcitation(detType const &source, double &pGet){

    // random excitation generator: randomly pick a connected determinant

    int exflag{0};
    //double pdoubles;
    detType target;
    std::vector<int> holes,particles;
    double rand{0.0};
    std::random_device rng;
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

    ExcitmatType exmat;
    // call the respective single or double excitation generators
    if (ic==2){
        // double excitation
        target = this->genDoubleExcitation(source, source_orbs, holes, particles, pGet, pdoubnew, noccs, nunoccs);
        // this part of exmat is only assigned for doubles (default is -1, singles keep that)
        exmat(1,0) = holes[1];
        exmat(1,1) = particles[1];
    }
    else{
        // single excitation
        target = this->genSingleExcitation(source, source_orbs, holes, particles, pGet, pdoubnew, noccs, nunoccs);
    }

    // assume particles + holes are not empty
    // TODO: Add exception safety
    exmat(0,0) = holes[0];
    exmat(0,1) = particles[0];

    // Using a uniform excitgen is like having a molecular Hamiltonian with all integrals being 1.0
    pBiasGen.checkProbabilities(exmat,1.0,1,pGet);

    return target;

}

//---------------------------------------------------------------------------------------------------//

detType UniformExcitgen::genSingleExcitation(detType const &source, std::vector<int> const &source_orbs,
		std::vector<int> &holes, std::vector<int> &particles, double &pgen, double pdoubnew,
		std::vector<int> noccs, std::vector<int> nunoccs) const{

    // randomly pick a single excitation (random excitation generator)

    int attempts{0},npairs{0},elecwnoexcit{0},nexcit{0};
    detType target;
    double rand{0.0};
    int iorb{0},aorb{0},ielec{0},aelec{0},ispin{0};
    int nel{0},norbs{0};
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
    int nexcita{0},nexcitb{0},nexcitotherway{0};
    bool baorbfail{true};
    int elecpairs{0},nel{0},norbs{0};
    int aorb{0},aspin{0},borb{0},bspin{0};
    int nforbiddenorbs{0};
    std::vector<int> elecs,spin;

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

double UniformExcitgen::getExcitationProb(detType const &source, detType const &target){

    // evaluate the generation probability of the uniform generation
	// from source to target

    double pgen;
    //double pdouble;
    double diff;
    int elecwnoexcit,ispin,nexcita,nexcitb,nexcitotherway;
    std::vector<int> holes,particles;
    std::vector<int> noccs,nunoccs;
    std::vector<int> spin;
    int elecpairs,norbs,nel;

    //pdouble = 0.0;
    pgen = 0.0;
    nexcita = 0;
    nexcitb = 0;
    nexcitotherway = 0;

    // get the particles and holes involved in this excitation
    for (size_t i=0; i<source.size(); ++i){
        diff = static_cast<int>(source[i]) - static_cast<int>(target[i]);
        if (diff > 0){
            // holes
            holes.push_back(i);
        }
        else if (diff < 0){
            // particles
            particles.push_back(i);
        }
    }

    if ((holes.size() != particles.size()) || (holes.size() > 2)){
        pgen = 0.0;

        return pgen;
    }

    if (holes.size() == 0){
        // source and target determinant are the same

        pgen = 0.0;

        return pgen;
    }

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

    // ensure that holes and particles have the same spin
    if (holes.size()==2){
        if ((holes[0]%2)!=(particles[0]%1)){
            // swap particles
            int tmp{0};
            tmp = particles[0];
            particles[0] = particles[1];
            particles[1] = tmp;
        }
    }

    norbs = source.size();
    nel = source_orbs.size();

    if (holes.size() == 1){
        // single excitaition

        // number of electrons with no available excitation
        elecwnoexcit = 0;
        for (size_t i=0; i<nunoccs.size(); ++i){
            if (nunoccs[i]==0){
                elecwnoexcit += nunoccs[i];
            }
        }

        // spin of the electron
        if ((holes[0]%2)==0){
            // alpha
            ispin = 0;
        }
        else{
            //beta spin
            ispin = 1;
        }

        // number of single excitation for electron i
        nexcita = nunoccs[ispin];

        // p_single = 1 - p_double
        // pgen = p_single * p(i) * p(a|i) * n/(n-elecwnoexcit)
        pgen = (1.0 - pDoubles) / (static_cast<double>(nexcita*(nel-elecwnoexcit)));

        return pgen;

    }
    else{
        // double excitaiton

        // number of electron pairs
        elecpairs = (nel*(nel-1))/2;

        // spin of electron
        for (size_t i=0; i< holes.size(); ++i){
            if ((holes[i]%2)==0){
                // alpha
                spin.push_back(0);
            }
            else{
                //beta spin
                spin.push_back(1);
            }
        }

        int nforbiddenorbs = countForbiddenOrbs(spin,nunoccs);

        // number of orbitals to select a from
        if (spin[0] != spin [1]){
            // alpha beta pair
            nexcita = norbs - nel;
        }
        else{
            if (spin[0]==0){
                // alpha
                nexcita = (norbs/2) - noccs[0];
            }
            else{
                //beta
                nexcita = (norbs/2) - noccs[1];
            }
        }

        // number of orbitals to select b from given that a has been selected
        if (spin[0] != spin[1]){
            // alpha beta pair
            if (spin[0]==0){
                // a is alpha
                nexcitb = nunoccs[1];
                nexcitotherway = nunoccs[0];
            }
            else{
                // a is beta
                nexcitb = nunoccs[0];
                nexcitotherway = nunoccs[1];
            }
        }
        else{
            // alpha alpha / beta beta
            nexcitb = nunoccs[spin[0]];
            nexcitotherway = nunoccs[spin[0]];

            // need to ensure that b is not picked again
            nexcitb -= 1;
            nexcitotherway -= 1;
        }

        // p_gen = p_double * p(ij) * [p(a|ij)p(b|a,ij) + p(b|ij)p(a|b,ij)]
        pgen = pDoubles*((1.0/static_cast<double>(nexcitb)) + (1.0/static_cast<double>(nexcitotherway))) /
            (static_cast<double>(elecpairs*(nexcita - nforbiddenorbs)));

        return pgen;
    }

    return pgen;

}

} /* namespace networkVMC */
