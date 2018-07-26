/*
 * AbInitioHamiltonian.cxx
 *
 *  Created on: Feb 13, 2018
 *      Author: Lauretta Schwarz
 */

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
#include <numeric>
#include <functional>
#include "../HilbertSpace/Determinant.hpp"
#include "AbInitioHamiltonian.hpp"
#include "../utilities/Errors.hpp"

namespace networkVMC{

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

AbInitioHamiltonian readAbInitioHamiltonian(int dim, std::string file_name, bool molpro_fcidump){
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
    int spin_types = 0;
    bool delimiter = false;
 
    std::cout << "Reading in 1- and 2-electron integrals from file "
        << file_name << std::endl;
    if (molpro_fcidump){
        std::cout << "Reading from a Molpro FCIDUMP file " << std::endl;
    }
    else{
        std::cout << "Reading from a normal FCIDUMP file " << std::endl;
    }

    if (inputfile.is_open() & inputfile.good()){
        while(!inputfile.eof()){
            std::getline(inputfile,line);
            //std::cout << line << std::endl;
            // split the line
            std::istringstream iss(line);
            std::vector<std::string> parts{
                std::istream_iterator<std::string>(iss), {}
            };
            // for testing
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
            if ((parts.size() > 0)&&(!end_header)){
                // find uhf=.false. or .true.
                for (size_t a=0; a<parts.size(); ++a){
                    std::string lower_parts = parts[a]; 
                    std::transform(lower_parts.begin(), lower_parts.end(), 
                            lower_parts.begin(), ::tolower);
                    if ((lower_parts.substr(0,3) == "uhf") || molpro_fcidump){
                        if (lower_parts == "uhf=.false."){
                            // integrals are in spatial orbitals in FCIDUMP files
                            // but stored in spin orbital basis
                            buhf = false;
                            H.initMatrixStorage(buhf);
                        }
                        else if (molpro_fcidump){
                            if ((!buhf) && (lower_parts.substr(0,6) == "iuhf=1")){
                                // integrals are in spin orbitals in FCIDUMP file
                                // and stored in spin orbital basis
                                buhf = true;
                                H.initMatrixStorage(buhf);
                            }
                            else if ((!buhf) && (lower_parts == "/")){
                                // complete header has been read in but no indication 
                                // for spin orbitals in FCIDUMP file has been found
                                // integrals are in spatial orbitals in FCIDUMP file
                                // and stored in spin orbital basis
                                buhf = false;
                                H.initMatrixStorage(buhf);
                            }
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
                    j = std::stoi(parts[2]);
                    k = std::stoi(parts[3]);
                    l = std::stoi(parts[4]);
                    val = std::stod(parts[0]);

                    if (molpro_fcidump){
                        // Molpro FCIDUMP
                        // Molpro uses the following ordering for spin orbital FCIDUMPS
                        // where a denotes an alpha spin and b a beta spin orbital
                        // 1: (aa|aa) <-> <aa|aa>
                        // 2: (bb|bb) <-> <bb|bb>
                        // 3: (aa|bb) <-> <ab|ab>
                        // 4: <a|h|a>
                        // 5: <b|h|b>
                        // with delimiters of 0.00000000 0 0 0 0

                        // check for delimiter
                        delimiter = false;
                        if ((i==0)&&(k==0)&&(j==0)&&(l==0)&&(std::fabs(val)<1e-10)){
                            // this is a delimiter
                            // move to next type of spin combination
                            spin_types += 1;
                            delimiter = true;
                        }

                        if (!delimiter){

                            if (spin_types == 0){
                                // <aa|aa>
                                i = (i*2) - 1;
                                j = (j*2) - 1;
                                k = (k*2) - 1;
                                l = (l*2) - 1;
                            }
                            else if (spin_types == 1){
                                // <bb|bb>
                                i *= 2;
                                j *= 2;
                                k *= 2;
                                l *= 2;
                            }
                            else if (spin_types == 2){
                                // <ab|ab>
                                i = (i*2) - 1;
                                j *= 2;
                                k = (k*2) - 1;
                                l *= 2;
                            }
                            else if (spin_types == 3){
                                // <a|h|a>
                                i = (i*2) - 1;
                                k = (k*2) - 1;
                            }
                            else if (spin_types == 4){
                                // <b|h|b>
                                i *= 2;
                                k *= 2;
                            }
                            else if (spin_types == 5){
                                // core energy
                            }
                            else{
                                std::cerr << "WARNING !" << std::endl;
                                std::cerr << "Something went wrong while reading in Molpro FCIDUMP file " << file_name << std::endl;
                                exit(1);
                            }

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
                            else if ((i==0)&&(k==0)&&(j==0)&&(l==0)){
                                // core energy
                                // E_core
                                H.setMatrixElement(val);
                            }
                        }
                    }
                    else{
                        // normal FCIDUMP
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
                        else if ((i==0)&&(k==0)&&(j==0)&&(l==0)){
                            // core energy
                            // E_core
                            H.setMatrixElement(val);
                        }
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
                    else if ((i==0)&&(k==0)&&(j==0)&&(l==0)){
                        // core energy
                        // E_core
                        H.setMatrixElement(val);
                    }
                }
            }
            // in order to start reading in the integrals
            std::string end_string = parts[0];
            std::transform(end_string.begin(),end_string.end(),end_string.begin(), ::tolower);
            if ((end_string.substr(1,3) == "end") || (molpro_fcidump && (end_string == "/"))){
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
    bool bfirsta,bfirstb,bnewij,bnewa;
    int orbi,orbj,orba,orbb;
    int elecpairs,spinpair;
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

}
