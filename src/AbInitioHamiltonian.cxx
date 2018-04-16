/*
 * AbInitioHamiltonian.cxx
 *
 *  Created on: Feb 13, 2018
 *      Author: Lauretta Schwarz
 */

#include "AbInitioHamiltonian.hpp"
#include "Determinant.hpp"
#include "utilities/Errors.hpp"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <random>

AbInitioHamiltonian::~AbInitioHamiltonian() {
}

int AbInitioHamiltonian::getFermiSign(detType const &alpha, int annihilatorIndex, int creatorIndex) const{
  //construct sign for conversion canonical shape (basisState 1(up),1(down),...,L(up),L(down)) to relative shape (a_exc^\dagger a_holes source)
  int fermiSign{0};
  int start{0}, end{0}, offset{0};
  fermiSign = 0;
  if(annihilatorIndex>creatorIndex){
    start=creatorIndex;
    end=annihilatorIndex;
    offset=1;
  }
  else{
    start=annihilatorIndex;
    end=creatorIndex;
    offset=1;
  }
  for(int k=start+offset;k<end;++k){
    if(alpha[k]){
      fermiSign += 1;
    }
  }
  return fermiSign;
}

AbInitioHamiltonian readAbInitioHamiltonian(int dim, std::string file_name){
    // this function reads in the 1- and 2-electron integrals of 
    // an ab-initio Hamiltonian from an FCIDUMP file

    AbInitioHamiltonian H(dim);

    // open the file in read mode
    std::ifstream inputfile;
    inputfile.open(file_name,std::ios::in);
    std::string line;
    int i,j,k,l;//,si,sj,sk,sl;
    double val = 0.0;
    bool buhf = false;
    bool end_header = false;
 
    std::cout << "Reading in 1- and 2-electron integrals from file "
        << file_name << std::endl;

    if (inputfile.is_open() & inputfile.good()){
        while(!inputfile.eof()){
            std::getline(inputfile,line);
            //std::cout << line << std::endl;
            // split the line
            std::istringstream iss(line);
            std::vector<std::string> parts{
                std::istream_iterator<std::string>(iss), {}
            };
            //// for testing
            //std::cout << "Number of strings in line: " <<
            //    parts.size() << std::endl;
            //for (std::vector<std::string>::iterator it=parts.begin(); 
            //        it != parts.end(); ++it){
            //    std::cout << *it << std::endl;
            //}

            // needed because .eof() does not work properly
            if (parts.size()==0) break;

            // FCIDUMP files are in chemical notation, (ik|jl), i.e.
            // (ik|jl) = <ij|kl>
            // but the integral storage is in physical notation <ij|kl>
            // need a case insensitive comparison of the string
            if ((parts.size() > 1)&&(!end_header)){
                // find uhf=.false. or .true.
                for (size_t a=0; a<parts.size(); ++a){
                    std::string lower_parts = parts[a]; 
                    std::transform(lower_parts.begin(), lower_parts.end(), 
                            lower_parts.begin(), ::tolower);
                    if (lower_parts.substr(0,3) == "uhf"){
                        if (lower_parts == "uhf=.false."){
                            // integrals are in spatial orbitals in FCIDUMP files
                            // but stored in spin orbital basis
                            buhf = false;
                            H.initMatrixStorage(buhf);
                        }
                        else {
                            // integrals are in spin orbitals in FCIDUMP file
                            // and stored in spin orbital basis
                            buhf = true;
                            H.initMatrixStorage(buhf);
                        }
                    }
                }
            }
            if (end_header) {
                if (buhf){
                    // integrals are in spin orbitals in FCIDUMP files
                    // and stored in spin orbital basis
                    // convert chemical to physical notation
                    i = std::stoi(parts[1]);
                    j = std::stoi(parts[3]);
                    k = std::stoi(parts[2]);
                    l = std::stoi(parts[4]);
                    val = std::stod(parts[0]);

                    if ((i!=0)&&(j!=0)&&(k!=0)&&(l!=0)){
                        // 2-electron integral
                        // <ij|kl> a+_i a+_j a_l a_k
                        // set all integrals which ought to be equal 
                        // by symmetry (real basis set: 8-fold)
                        // <ij|kl>
                        H.setMatrixElement(i,j,k,l,val);
                        //// <ij|kl> = <ji|lk> = <kl|ij> = <lk|ji> =
                        //// <kj|il> = <li|jk> = <il|kj> = <jk|li>
                    }
                    else if ((i!=0)&&(k!=0)&&(j==0)&&(l==0)){
                        // 1-electron integral
                        // <i|h|k> a+_i a_k
                        // set all integrals which ought to be equal
                        // by symmetry (real basis set: 2-fold)
                        // <i|h|k> 
                        H.setMatrixElement(i,k,val);
                        // <i|h|k> = <k|h|i>
                    }
                }
                else{
                    // integrals are in spatial orbitals in FCIDUMP file
                    // but stored in spin orbital basis
                    // convert chemical to physical notation
                    i = std::stoi(parts[1]);
                    j = std::stoi(parts[3]);
                    k = std::stoi(parts[2]);
                    l = std::stoi(parts[4]);
                    val = std::stod(parts[0]);
                    // integrals / indices refer to spatial orbitals now 

                    if ((i!=0)&&(j!=0)&&(k!=0)&&(l!=0)){
                        // 2-electron integral
                        // <ij|kl> a+_i a+_j a_l a_k
                        // <ij|kl> 
                        H.setMatrixElement(i,j,k,l,val);
                        // set all integrals which ought to be equal 
                        // by symmetry (real basis set: 8-fold)
                        // <ij|kl> = <ji|lk> = <kl|ij> = <lk|ji> =
                        // <kj|il> = <li|jk> = <il|kj> = <jk|li>
                    }
                    else if ((i!=0)&&(k!=0)&&(j==0)&&(l==0)){
                        // 1-electron integral
                        // <i|h|k> a+_i a_k
                        // <i|h|k> 
                        H.setMatrixElement(i,k,val);
                        // set all integrals which ought to be equal
                        // by symmetry (real basis set: 2-fold)
                        // <i|h|k> = <k|h|i>
                    }
                }
            }
            // in order to start reading in the integrals
            std::string end_string = parts[0];
            std::transform(end_string.begin(),end_string.end(),end_string.begin(), ::tolower);
            if (end_string.substr(1,3) == "end"){
                end_header = true;
            }
        }
    }
    else {
        std::cerr << "WARNING !" << std::endl;
        std::cerr << "Could not read in file " << file_name << std::endl;

        exit(1);
    }

    // close the file
    inputfile.close();

    return H;
}


detType AbInitioHamiltonian::getSingleExcitation(detType const &source, std::vector<int> const &source_orbs, std::vector<int> &holes, std::vector<int> &particles, int &exflag, bool &ballexcitfound) const{

    // get the next single excitation (deterministic excitation generator)

    int orbi,orba;
    int nel,norbs;
    bool binitorbsfound;
    detType target=source;
    static int orbiind = 0;
    static int orbaind = 0;

    // number of electrons and spin orbitals
    nel = source_orbs.size();
    norbs = source.size();

    binitorbsfound = false;

    // use the spin orbitals of the last excitation
    if (holes.size()>0){
        orbi = holes[0];
        orba = particles[0];
    }
    else{
        orbi = -1;
        orba = -1;
    }


    // the first excitation ?
    if ((orbi==-1)||(orba==-1)){
        // the first occupied spin orbital
        orbiind = 0;
        // orbitals start from 0
        orbi = source_orbs[orbiind];
        if ((orbi%2)==0){
            // alpha spin orbital
            orba = 0;
            orbaind = 0;
        }
        else{
            // beta spin orbital
            orba = 1;
            orbaind = 1;
        }
    }
    else{
        // start with the same orbital i of the last excitation 
        // and check whether there are any more possible excitations
        orbi = holes[0];
        orba = particles[0];
        if ((((orbi%2)==0)&&(orbaind==(norbs-2)))||(((orbi%2)==1)&&
                    (orbaind==(norbs-1)))){
            // orbital a is the last spin orbital 
            // -? need to move to next i and a
            orbiind += 1;
            if (orbiind<=(nel-1)){
                orbi = source_orbs[orbiind];
                if ((orbi%2)==0){
                    // alpha spin orbital
                    orba = 0;
                    orbaind = 0;
                }
                else{
                    // beta spin orbital
                    orba = 1;
                    orbaind = 1;
                }
            }
            else{
                // reached the end of the determinant
                if (exflag != 1){
                    // move to double excitations
                    exflag = 2;
                    holes.clear();
                    particles.clear();
                    orbiind = 0;
                    target.clear();
                    target = this->getDoubleExcitation(source,source_orbs,holes,particles,exflag,ballexcitfound);

                    return target;
                }
                else{
                    ballexcitfound = true;
                    binitorbsfound = true;
                    target.clear();

                    return target;
                }
            }
        }
        else{
            // there are more possible excitations from orbital i
            // to orbital a
            orbaind += 2;
        }
    }
    

    while (!binitorbsfound){

        int bendaorbs=false;

        if (orbiind>(nel-1)){
            // last single excitation has been reached
            // find first double excitaiton
            if (exflag != 1){
                // move to double excitations
                exflag = 2;
                holes.clear();
                particles.clear();
                orbiind = 0;
                target.clear();
                target = this->getDoubleExcitation(source,source_orbs,holes,particles,exflag,ballexcitfound);

                return target;
            }
            else{
                ballexcitfound = true;
                binitorbsfound = true;
                target.clear();

                return target;
            }
        }

        // to find orbital a take the first with the appropriate spin
        if (orbaind > (norbs-1)){
            // orbital a is the last orbital in the basis set
            bendaorbs = true;
        }
        else{
            // this is the appropriate obital a
            bendaorbs = false;
            orba = orbaind;
        }

        // check that orbital a is not occupied in the determinant
        int no_occ=0;
        if (!bendaorbs){
            while (source[orba]){
                // orbital a is occupied in the determinant
                no_occ += 2;
                if ((orbaind+no_occ) > (norbs-1)){
                    // end of all orbitals a 
                    // need to move to a new i
                    bendaorbs = true;
                    break;
                }
                else{
                    // new orbital a
                    orba = orbaind + no_occ;
                    if ((((orba%2)==0)&&(orba > (norbs-2)))||
                            (((orba%2)==1)&&(orba > (norbs-1)))){
                        break;
                    }
                }
            }
        }

        // check for spin symmetry
        if (((orbi%2)==(orba%2))&&(!bendaorbs)){
            // these are the new orbitals
            binitorbsfound = true;
            orbaind += no_occ;
        }
        else{
            // move to next i
            orbiind += 1;
            if (orbiind<=(nel-1)){
                orbi = source_orbs[orbiind];
                if ((orbi%2)==0){
                    // alpha spin orbital
                    orba = 0;
                    orbaind = 0;
                }
                else{
                    // beta spin orbital
                    orba = 1;
                    orbaind = 1;
                }
            }
        }
    }

    if (!ballexcitfound){

        // make the new single excitation
        target = source;
        annihilate(target,orbi);
        create(target,orba);
        holes.clear();
        particles.clear();
        holes.push_back(orbi);
        particles.push_back(orba);
    }


    return target;

}

detType AbInitioHamiltonian::getDoubleExcitation(detType const &source, std::vector<int> const &source_orbs, std::vector<int> &holes, std::vector<int> &particles, int &exflag, bool &ballexcitfound) const{
   
    // generate the next double excitation (deterministic excitation generator)

    bool bdoubleexcitfound=false;
    bool bfirsta,bfirstb,bnewij,bnewa,bnewb;
    int orbi,orbj,orba,orbb;
    int elecpairs,indij,spinpair;
    int el1{0},el2{0},spinb{0};
    int nel,norbs;
    static int ijind=0;
    static int orbbind=0;
    static int spina=0;
    static int orbachosen=0;
    detType target=source;

    // number of electrons and spin orbitals
    nel = source_orbs.size();
    norbs = source.size();

    bdoubleexcitfound = false;
    bfirsta = false;
    bfirstb = false;
    bnewa = false;
    bnewb = false;
    
    // orbitals of the last excitation
    if (holes.size()>0){
        orbi = holes[0];
        orbj = holes[1];
        orba = particles[0];
        orbb = particles[1];
    }
    else{
        orbi = -1;
        orbj = -1;
        orba = -1;
        orbb = -1;
    }

    // number of electron pairs
    elecpairs = nel*(nel-1)/2;

    if (orbi==-1){
        // first double excitation
        ijind = 1;
        bfirsta = true;
        bfirstb = true;
    }
    

    while (!bdoubleexcitfound){
        // use the previous ijind and save the indices for a and b
        // -> this routine pick an electron pair i,j, specified
        // by the index ijind
        // the i and j orbitals are ni(el1) and ni(el2) and the 
        // spin of the pair is spinpair
        this->getElecPair(source_orbs, el1, el2, spinpair, ijind);
        // spinpair = 1 -> alpha + alpha
        // spinpair = 2 -> alpha + beta
        // spinpair = 3 -> beta + alpha
        // spinpair = 4 -> beta + beta
        indij = (((source_orbs[el2]+1-2)*(source_orbs[el2]+1-1))/2) + source_orbs[el1]+1;
        // this becomes true when there is no longer an allowed orbital a for this i,j 
        // pair and the next one has to be considered
        bnewij = false;
        // this loop runs through allowed a orbitals until a double excitaton is found
        while ((!bnewij)&&(!bdoubleexcitfound)){
            // if this is the first a,b pair which is picked with this i,j pair, start with 
            // the alpha spin, unless i and j are both beta
            if (bfirsta){
                if (spinpair==4){
                    // beta + beta
                    spina = 1;
                    // start with first beta orbital
                    orbachosen = 1;
                }
                else if (spinpair==1){
                    // alpha + alpha
                    spina = 0;
                    // start with first alpha orbital
                    orbachosen = 0;
                }
                else{
                    spina = 0;
                    // start with first alpha orbital
                    orbachosen = 0;
                }
            }
            // if this is not the first excittion, it is checked whether there are any more excitations
            // remaingin by starting with the previous spina and orba index
            orba = orbachosen;
            // test that this spin orbital is unoccupied
            while (source[orba]){
                // orbital a is occupied
                if ((spinpair!=2)&&(spinpair!=3)){
                    // for alpha aloha or beta beta pairs need to look at the 
                    // same spin state
                    orbachosen += 2;
                }
                else{
                    // mixed spin pairs need to look at both alpha and beta
                    orbachosen += 1;
                    if (spina==0){
                        spina = 1;
                    }
                    else{
                        spina = 0;
                    }
                }

                if (orbachosen > (norbs-1)){
                    // reached the end of the basis set
                    // -> new ij pair
                    bnewij = true;
                    break;
                }

                // otherwise orbital a is the first unoccupied spin orbital
                orba = orbachosen;
            }


            // if the end of orbitals a has been reached, a new ij pair is needed
            // igind needs to be incremented and checked that it is not above the limit
            if (bnewij){
                ijind += 1;
                // new ij pair is above the limit ?
                if (ijind > elecpairs){
                    bdoubleexcitfound = true;
                    // there are no more allowed excitaitons
                    ballexcitfound = true;
                }
                break;
            }

            bnewa = false;
            // find the ab pair
            while ((!bnewa)&&(!bdoubleexcitfound)){
                // the ij and a orbitals have been found
                // -> need to pick final orbital b 
                // first find spin of b
                if (spinpair==1){
                    // alpha + alpha
                    spinb = 0;
                }
                else if (spinpair==4){
                    // beta beta
                    spinb = 1;
                }
                else{
                    // alpha beta / beta alpha pair
                    if (spina==1){
                        // alpha spin
                        spinb = 0;
                    }
                    else{
                        // beta spin
                        spinb = 1;
                    }
                }
                // if this is the first time taht an orbital b is picked for ij and a
                // begin at the start of orbitals
                // otherwise pick up where things werwe left off last time
                if (bfirstb){
                    if (spinb==0){
                        // alpha spin
                        orbbind = 0;
                    }
                    else{
                        // beta spin
                        orbbind = 1;
                    }
                }
                else{
                    // increment by 2 to retain spin symmetry
                    orbbind += 2;
                }

                // if the new b orbital is still within the basis set, check 
                // that it is unoccupied, if not move to next orbitals
                if (orbbind > (norbs-1)){
                    // choose new a orbital
                    bnewa = true;
                    bfirsta = false;
                }

                if (!bnewa){
                    // check that for the given spin orbital b is no out of range
                    if (((spinb==0)&&(orbbind > (norbs-2)))||((spinb==1)&&(orbbind > (norbs-1)))){
                        // out of range of alpha and beta spin orbitals
                        bnewa = true;
                        bfirsta = false;
                    }
                    else{
                        orbb = orbbind;
                        // check that b is unoccupied
                        while (source[orbb]||(orbb<=orba)){
                            // orbital is occupied 
                            // -> move one up
                            orbbind += 2;
                            if (((spinb==0)&&(orbbind > (norbs-2)))||((spinb==1)&&(orbbind > (norbs-1)))){
                                // b has gone out of range
                                bnewa = true;
                                bfirsta = false;
                                break;
                            }
                            // update orbital b
                            orbb = orbbind;
                        }
                    }
                }

                // if a new orbital a is required nned to check that a new ij pair is not 
                // also needed
                if (bnewa){
                    if ((spinpair==1)||(spinpair==4)){
                        // increment by 2 because the same spin is needed
                        orbachosen += 2;
                    }
                    else{
                        // only increment by 1 since both alpha and beta spin need to be looked at
                        if (spina==1){
                            spina = 0;
                        }
                        else{
                            spina = 1;
                        }
                        orbachosen += 1;
                    }
                    // need the first b orbital
                    bfirstb = true;
                    // moved out of range of basis set ?
                    if (orbachosen > (norbs-1)){
                        // need new ij pair
                        bnewij = true;
                        ijind += 1;
                        // moved out of range ?
                        if (ijind > elecpairs){
                            ballexcitfound = true;
                            bdoubleexcitfound = false;
                            // need to break out of all nested loops
                            goto found_double_excit;
                        }
                    }
                }
                else{
                    // allowed excitation
                    bdoubleexcitfound = true;
                }
            }
        }
        // this is the loop for the new ij pair 
        // if a new ij pair is needed -> automatically choose a 
        // new a and b
        bfirsta = true;
        bfirstb = true;
    }

found_double_excit:

    target = source;

    // make the new double excitation
    if (bdoubleexcitfound&&(!ballexcitfound)){
        orbi = source_orbs[el1];
        orbj = source_orbs[el2];
        target = source;
        annihilate(target,orbi);
        annihilate(target,orbj);
        create(target,orba);
        create(target,orbb);
        holes.clear();
        particles.clear();
        holes.push_back(orbi);
        holes.push_back(orbj);
        particles.push_back(orba);
        particles.push_back(orbb);
    }
    else{
        // no double excitation found
        target.clear();
        holes.clear();
        particles.clear();
    }

    return target;


}


void AbInitioHamiltonian::getElecPair(std::vector<int> const &source_orbs, 
        int &el1, int &el2, int &spinpair, int &ind) const{
    // return a pair of electrons (el1,el2) in the source determinant
    // uses a triangular indexing system in order to pick out two 
    // distinct electron 

    int elecpairs,x,k,orb1,orb2;
    int nel;

    // number of electrons and spin orbitals
    nel = source_orbs.size();

    elecpairs = nel*(nel-1)/2;
    x = elecpairs - ind;
    k = static_cast<int>((std::sqrt(8.0*static_cast<double>(x)+1.0)-1.0)/2.0);
    el1 = nel - 1 - k;
    el2 = nel - x + ((k*(k+1))/2);

    // in order to move indices in range 0,..,m-1
    el1 -= 1;
    el2 -= 1;

    // the orbitals
    orb1 = source_orbs[el1];
    orb2 = source_orbs[el2];

    // which type of spin pairs
    // spinpair = 1 -> alpha + alpha
    // spinpair = 2 -> alpha + beta
    // spinpair = 3 -> beta + alpha
    // spinpair = 4 -> beta + beta
    spinpair = 0;
    if (((orb1%2)==0)&&((orb2%2)==0)){
        // alpha + alpha
        spinpair = 1;
    }
    else if (((orb1%2)==0)&&((orb2%2)==1)){
        // alpha + beta
        spinpair = 2;
    }
    else if (((orb1%2)==1)&&((orb2%2)==0)){
        // beta + alpha
        spinpair = 3;
    }
    else if (((orb1%2)==1)&&((orb2%2)==1)){
        // beta + beta
        spinpair = 4;
    }
}


std::vector<detType> AbInitioHamiltonian::getCoupledStates(detType const &source) const{
   
    // deterministic excitation generator: returns a list of all connected determinants

    int exflag;
    detType target;
    bool ballexcitfound=false;
    std::vector<int> holes,particles;
    std::vector<detType> coupledList;

    // get the occupied spin orbitals
    std::vector<int> source_orbs = getOccupiedPositions(source);
    // start with single excitations
    exflag = 3;

    while (!ballexcitfound){
        if ((exflag==1)||(exflag==3)){
            target = this->getSingleExcitation(source, source_orbs, holes, particles, exflag, ballexcitfound);
        }
        else if (exflag==2){
            target = this->getDoubleExcitation(source, source_orbs, holes, particles, exflag, ballexcitfound);
        }
        if (ballexcitfound){
            break;
        }
        coupledList.push_back(target);
    }

    return coupledList;

}

detType AbInitioHamiltonian::genSingleExcitation(detType const &source, std::vector<int> const &source_orbs, std::vector<int> &holes, std::vector<int> &particles, double &pgen, double pdoubnew, std::vector<int> noccs, std::vector<int> nunoccs) const{

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


std::vector<int> AbInitioHamiltonian::pickElecPair(std::vector<int> const &source_orbs, std::vector<int> &spin, int &elecpairs) const{

    // use a triangular indexing system in order to pick one electron pairs using just one single random 
    // number for each pair

    int nel;
    double rand;
    std::vector<int> elecs;
    std::vector<int> orbs;
    int el1,el2;
    std::random_device rng;
    double const normalisation=static_cast<double>(rng.max()); // minimal potentially returned values is 0

    // number of electrons
    nel = source_orbs.size();
    // number of electron pairs
    elecpairs = (nel*(nel-1))/2;

    // pick a pair
    rand = static_cast<double>(rng())/normalisation;
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
    for (size_t i=0; i<orbs.size(); ++i){
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


int AbInitioHamiltonian::countForbiddenOrbs(std::vector<int> const &spins, std::vector<int> nunoccs) const{

    // count the number of forbidden spin orbitals (random excitation generator)

    int nforbiddenorbs;

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


int AbInitioHamiltonian::pickOrbA(detType const &source, std::vector<int> const &spins, std::vector<int> noccs, std::vector<int> nunoccs, int norbs, int nel, int nforbiddenorbs, int &nexcita, int &aspin, bool &baorbfail) const{

    // pick orbital a in an ij->ab excitation (random excitation generator)

    int attemptsa,attemptsoverall;
    int aelec;
    int aorb,counter;
    double rand;
    std::random_device rng;
    double const normalisation=static_cast<double>(rng.max()); // minimal potentially returned values is 0

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
                rand = static_cast<double>(rng())/normalisation;
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
                    rand = static_cast<double>(rng())/normalisation;
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
                rand = static_cast<double>(rng())/normalisation;
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
                    rand = static_cast<double>(rng())/normalisation;
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


int AbInitioHamiltonian::pickOrbB(detType const &source, std::vector<int> const &spins, std::vector<int> noccs, std::vector<int> nunoccs, int norbs, int nel, int nforbiddenorbs, int aorb, int aspin, int &nexcitb, int &bspin, int &nexcitotherway) const{

    // pick orbital b in an ij->ab excitation (random excitation generator)

    int borb;
    int attemptsb,counter;
    double rand;
    int belec;
    std::random_device rng;
    double const normalisation=static_cast<double>(rng.max()); // minimal potentially returned values is 0

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

        rand = static_cast<double>(rng())/normalisation;
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
            rand = static_cast<double>(rng())/normalisation;
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

detType AbInitioHamiltonian::genDoubleExcitation(detType const &source, std::vector<int> const &source_orbs, std::vector<int> &holes, std::vector<int> &particles, double &pgen, double pdoubnew, std::vector<int> noccs, std::vector<int> nunoccs) const{

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
    elecs = this->pickElecPair(source_orbs, spin, elecpairs); 

    // the number of forbidden orbitals which are not allowed by spin
    nforbiddenorbs = this->countForbiddenOrbs(spin, nunoccs);

    // pick orbital a
    aorb = this->pickOrbA(source, spin, noccs, nunoccs, norbs, nel, nforbiddenorbs, nexcita, aspin, baorbfail);

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
    borb = this->pickOrbB(source, spin, noccs, nunoccs, norbs, nel, nforbiddenorbs, aorb, aspin, nexcitb, bspin, nexcitotherway);

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


detType AbInitioHamiltonian::getRandomCoupledState(detType const &source, double &pGet) const{

    // random excitation generator: randomly pick a connected determinant

    int exflag;
    double pdoubles;
    detType target;
    std::vector<int> holes,particles;
    double rand;
    std::random_device rng;
    int nsingleexcit,ndoubleexcit,nexcit;
    double const normalisation=static_cast<double>(rng.max()); // minimal potentially returned values is 0
    std::vector<int> noccs,nunoccs;

    exflag = 3;
    pdoubles = 0.0;
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
    pdoubles = static_cast<double>(ndoubleexcit)/static_cast<double>(nexcit);

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
        if (rand<pdoubles){
            ic = 2;
        }
        else{
            ic = 1;
        }
        pdoubnew = pdoubles;
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

int AbInitioHamiltonian::countNumberCoupledStates(detType const &source, int exflag, int &nsingleexcit, int &ndoubleexcit){

    // count the number of connected determinants

    int ncoupledstates=0;
    int nalpha,nbeta,nunoccalpha,nunoccbeta;

    if ((exflag<1)||(exflag>3)){
        throw OutOfRangeError(exflag);
    }

    // number of occupied and unoccupied alpha and beta spin orbitals
    nalpha = 0;
    nbeta = 0;
    for (size_t i=0; i<source.size(); ++i){
        if (source[i]){
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
    nunoccalpha = (source.size()/2) - nalpha;
    nunoccbeta = (source.size()/2) - nbeta;


    if ((exflag==1)||(exflag==3)){
        
        // single excitations
        nsingleexcit = (nalpha*nunoccalpha) + (nbeta*nunoccbeta);

    }
    if ((exflag==2)||(exflag==3)){
        // double excitations

        // alpha-alpha, beta-beta, alpha-beta pairs
        ndoubleexcit = ((nalpha*(nalpha-1)/2)*(nunoccalpha*(nunoccalpha-1)/2)) + 
            ((nbeta*(nbeta-1)/2)*(nunoccbeta*(nunoccbeta-1)/2)) +
            (nalpha*nbeta*nunoccalpha*nunoccbeta);
    }

    if (exflag==1){
        ndoubleexcit = 0;
    }
    if (exflag==2){
        nsingleexcit = 0;
    }

    ncoupledstates = nsingleexcit + ndoubleexcit;

    return ncoupledstates;

}

double AbInitioHamiltonian::calcGenProp(detType const &source, detType const &target){

    // evaluate the generation probability

    double pgen;
    double pdouble;
    double diff;
    int nexcit,nsingleexcit,ndoubleexcit,exflag;
    int elecwnoexcit,ispin,nexcita,nexcitb,nexcitotherway;
    std::vector<int> holes,particles;
    std::vector<int> noccs,nunoccs;
    std::vector<int> spin;
    int elecpairs,norbs,nel;

    pdouble = 0.0;
    pgen = 0.0;
    nexcita = 0;
    nexcitb = 0;
    nexcitotherway = 0;


    exflag = 3;

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
    pdouble = static_cast<double>(ndoubleexcit)/static_cast<double>(nexcit);



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
        pgen = (1.0 - pdouble) / (static_cast<double>(nexcita*(nel-elecwnoexcit)));

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

        int nforbiddenorbs = this->countForbiddenOrbs(spin,nunoccs);

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
        pgen = pdouble*((1.0/static_cast<double>(nexcitb)) + (1.0/static_cast<double>(nexcitotherway))) / 
            (static_cast<double>(elecpairs*(nexcita - nforbiddenorbs)));

        return pgen;
    }

    return pgen;

}
