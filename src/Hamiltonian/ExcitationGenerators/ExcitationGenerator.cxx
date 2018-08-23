/*
 * ExcitationGenerator.cxx
 *
 *  Created on: Jun 15, 2018
 *      Author: Lauretta Schwarz (functions moved here by guther)
 */

#include <vector>
#include <cmath>
#include "../../utilities/Errors.hpp"
#include "../../utilities/RNGWrapper.hpp"
#include "../../utilities/TypeDefine.hpp"

// auxiliary functions for excitation generation

namespace networkVMC{

//---------------------------------------------------------------------------------------------------//

int countForbiddenOrbs(std::vector<int> const &spins, std::vector<int> nunoccs){

    // count the number of forbidden spin orbitals (required in random excitation generator)

    int nforbiddenorbs{0};

    nforbiddenorbs = 0;

    if (spins[0]!=spins[1]){
        // alpha beta pair
        if (nunoccs[0]==0){
            // there are no allowed alpha spin orbitals
            // are there any beta spin orbitals which are
            // thus forbidden ?
            nforbiddenorbs += nunoccs[1];
        }
        if (nunoccs[1]==0){
            // there are no allowed beta spin orbitals
            // are there any alpha spin orbitals which are
            // thus forbidden?
            nforbiddenorbs += nunoccs[0];
        }
    }
    else{
        if (spins[0]==0){
            // alpha alpha pair
            // if there is only one allowed alpha orbital a pair cannot be chosen
            if (nunoccs[0] == 1){
                nforbiddenorbs = 1;
            }
        }
        else{
            // beta beta pair
            // if there is only one allowed beta orbital a pair cannot be chosen
            if (nunoccs[1] == 1){
                nforbiddenorbs = 1;
            }
        }
    }


    return nforbiddenorbs;

}

//---------------------------------------------------------------------------------------------------//

std::vector<int> pickElecPair(std::vector<int> const &source_orbs, std::vector<int> &spin, int &elecpairs){

    // use a triangular indexing system in order to pick one electron pairs using just one single random
    // number for each pair

    int nel{0};
    double rand{0.0};
    std::vector<int> elecs;
    std::vector<int> orbs;
    int el1{0},el2{0};
    RNGWrapper rng;

    // number of electrons
    nel = source_orbs.size();
    // number of electron pairs
    elecpairs = (nel*(nel-1))/2;

    // pick a pair
    rand = rng();
    int ind = 1 + static_cast<int>(static_cast<double>(elecpairs)*rand);

    // generate the two indices and corresponding orbitals
    el1 = std::ceil((1+ std::sqrt(1+8*static_cast<double>(ind)))/2);
    el2 = ind - ((el1 - 1) *(el1 - 2)) / 2;
    // from indices from range 1,... to 0,...
    el1 -= 1;
    el2 -= 1;

    elecs.push_back(el1);
    elecs.push_back(el2);
    orbs.push_back(source_orbs[elecs[0]]);
    orbs.push_back(source_orbs[elecs[1]]);

    // and the spin of the electrons
    for (std::size_t i=0; i<orbs.size(); ++i){
        if ((orbs[i]%2)==0){
            // alpha spin
            spin.push_back(0);
        }
        else{
            // beta spin
            spin.push_back(1);
        }
    }

    return elecs;

}

//---------------------------------------------------------------------------------------------------//

int pickOrbA(detType const &source, std::vector<int> const &spins, std::vector<int> noccs,
		std::vector<int> nunoccs, int norbs, int nel, int nforbiddenorbs, int &nexcita,
		int &aspin, bool &baorbfail){

    // pick orbital a in an ij->ab excitation (random excitation generator)

    int attemptsa{0},attemptsoverall{0};
    int aelec{0};
    int aorb{0},counter{0};
    double rand{0.0};
    RNGWrapper rng;

    aorb = -1;

    nexcita = 0;
    // pick the orbital a, i.e. the first unoccupied orbital
    if (spins[0]!=spins[1]){
        // alpha beta pair
        // -> no restrictions on the spin or orbital a
        nexcita = norbs - nel;
    }
    else{
        if (spins[0]==0){
            // alpha alpha pair
            nexcita = (norbs/2) - noccs[0];
        }
        else{
            // beta beta pair
            nexcita = (norbs/2) - noccs[1];
        }
    }

    // if there are no allowed orbitals left
    if (nexcita==nforbiddenorbs){
        baorbfail = true;

        return aorb;

    }
    else{
        baorbfail = false;
    }

    attemptsa = 0;
    attemptsoverall = 0;

    while (true){
        // keep drawing unoccupied orbitals until one with
        // an allowed partner is found

        if (spins[0]!=spins[1]){
            // alpha beta pair
            // randomly choose spin
            if ((nexcita-nforbiddenorbs)<3){
                // run through all orbitals to find the desired one

                // unoccupied orbital a
                rand = rng();
                aelec = static_cast<int>((static_cast<double>(nexcita-nforbiddenorbs)*rand)) + 1;

                counter = 1;
                for (size_t i=0; i<source.size(); ++i){
                    // is it not occupied ?
                    if (!source[i]){
                        // orbital not in determinant and thus it is allowed
                        if ((i%2)==0){
                            // alpha spin
                            aspin = 0;
                        }
                        else{
                            // beta spin
                            aspin = 1;
                        }
                        if (counter==aelec){
                            // found the allowed orbital
                            aorb = i;
                            return aorb;
                        }
                        counter += 1;
                    }
                }
                // unable to find orbital a
                throw OutOfRangeError(aelec);
            }
            else{
                // keep drawing orbitals randomly until an
                // unoccupied one is found
                attemptsa = 0;
                while (true){
                    // draw random orbital
                    rand = rng();
                    aorb = static_cast<int>(static_cast<double>(norbs)*rand);
                    // test if the orbital is occupied
                    if (!source[aorb]){
                        break;
                    }

                    if (attemptsa > 250){
                        // unable to find orbital a
                        throw OutOfRangeError(attemptsa);
                    }
                    attemptsa += 1;
                }
            }

            if ((aorb%2)==0){
                // alpha orbital
                aspin = 0;
            }
            else{
                // beta orbital
                aspin = 1;
            }
        }
        else{
            // either alpha alpha or beta beta pair
            if ((nexcita-nforbiddenorbs) < 3){
                // run through all orbitals until the desired one
                // from the (nexcit-forbidden orbs) is found
                rand = rng();
                aelec = static_cast<int>(static_cast<double>(nexcita-nforbiddenorbs)*rand) + 1;
                aspin = spins[0];

                counter = 1;
                for (int i=1; i<((norbs/2)+1); ++i){
                    if (spins[0]==0){
                        // alpha spin
                        aorb = (2*i) - 2;
                    }
                    else{
                        // beta spin
                        aorb = (2*i) - 1;
                    }

                    // test if allowed
                    if (!source[aorb]){
                        // orbital is not occupied
                        if (counter==aelec){
                            // found orbital

                            return aorb;
                        }
                        counter += 1;
                    }
                }
                throw OutOfRangeError(aelec);
            }
            else{
                // draw orbitals randomly until an unoccupied one is found
                attemptsa = 0;

                while (true){
                    // draw orbitals randomly from set of orbital
                    rand = rng();
                    aorb = 2*(1 + static_cast<int>(static_cast<double>(norbs/2)*rand)) + (spins[0]-2);
                    aspin = spins[0];

                    // test if orbital is occupied
                    if (!source[aorb]){
                        break;
                    }

                    if (attemptsa > 250){
                        throw OutOfRangeError(attemptsa);
                    }
                    attemptsa += 1;
                }
            }
        }

        // need to test whether there are unoccupied orbitals available for the
        // second particle
        if (spins[0]!=spins[1]){
            // alpha beta spin
            if (aspin==0){
                // a has alpha spin
                if (nunoccs[1]!=0){
                    // a is allowed
                    break;
                }
            }
            else{
                // a has beta spin
                if (nunoccs[0]!=0){
                    // a is allowed
                    break;
                }
            }
        }
        else{
            if (nunoccs[spins[0]]!=0){
                // orbital a is allowed
                break;
            }
        }

        if (attemptsoverall > (norbs*20)){
            throw OutOfRangeError(attemptsoverall);
        }
        attemptsoverall += 1;
    }

    return aorb;

}

//---------------------------------------------------------------------------------------------------//

int pickOrbB(detType const &source, std::vector<int> const &spins, std::vector<int> noccs,
		std::vector<int> nunoccs, int norbs, int nel, int nforbiddenorbs, int aorb,
		int aspin, int &nexcitb, int &bspin, int &nexcitotherway){

    // pick orbital b in an ij->ab excitation (random excitation generator)

    int borb{0};
    int attemptsb{0},counter{0};
    double rand{0.0};
    int belec{0};
    RNGWrapper rng;

    borb = -1;


    // find out what spin b has to have
    if (spins[0]!=spins[1]){
        // alpha beta pair
        if (aspin==0){
            // a has alpha spin
            // b needs to pick a beta spin orbital
            nexcitb = nunoccs[1];
            nexcitotherway = nunoccs[0];
            bspin = 1;
        }
        else{
            // a has beta spin
            // need to pick an alpha spin orbital for b
            nexcitb = nunoccs[0];
            nexcitotherway = nunoccs[1];
            bspin = 0;
        }
    }
    else{
        // alpha alpha or beta beta spin pair
        nexcitb = nunoccs[spins[0]];
        nexcitotherway = nunoccs[spins[0]];
        bspin = spins[0];

        // in this case need to ensure that b is not
        // picked ti be the same as a
        // need to redraw in this case
        nexcitb -= 1;
        nexcitotherway -= 1;
    }


    if (false){
        // run though all orbials of the same spin to find an
        // allowed one

        rand = rng();
        belec = static_cast<int>(static_cast<double>(nexcitb)*rand) + 1;

        counter = 0;
        for (int i=0; i<(norbs/2); ++i){
            borb = (2*(i+1)) + (bspin-2);
            if ((!source[borb]&&(borb!=aorb))){
                // orbital is not occupied
                counter += 1;
                if (counter==belec){
                    // found b
                    break;
                }
            }
        }

        if (counter != belec){
            throw OutOfRangeError(belec);
        }
    }
    else{
        // keep drawing orbitals with the same spin until an
        // unoccupied one is found

        attemptsb = 0;

        while (true){

            // random orbitals
            rand = rng();
            borb = 2*(1 + static_cast<int>(static_cast<double>(norbs/2)*rand)) + (bspin-2);

            if ((!source[borb])&&(borb!=aorb)){
                // orbital is not occupied
                // orbital b is found
                break;
            }

            if (attemptsb>1000){
                throw OutOfRangeError(attemptsb);
            }
            attemptsb += 1;
        }
    }


    return borb;

}

}



