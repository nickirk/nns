#include "Hamiltonian.hpp"
#include <random>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>


int Hamiltonian::getId(int i) const {
    // convert a spin orbital index to a spatial one if necessary
    int ind{0};

    if (spin_orbs){
        // stored in spin orbtials
        ind = i;
        
        return ind;
    }
    else{
        // stored in spatial orbitals
        ind = (i-1)/2 + 1;

        return ind;
    }

}

//---------------------------------------------------------------------------------------------------//


void Hamiltonian::initMatrixStorage(bool bspin_orbs){
    // are integrals stored in spin or spatial orbitals
    
    spin_orbs = false;
    spin_orbs = bspin_orbs;
    if (spin_orbs){
        std::cout << "Storing 1- and 2-electron integrals in spin orbitals" << std::endl;
    }
    else{
        std::cout << "Storing 1- and 2-electron integrals in spatial orbitals" << std::endl;
    }

}

//---------------------------------------------------------------------------------------------------//


void Hamiltonian::setMatrixElement(int p, int q, double newEntry){
    // <p|h|q> a+_p a_q
    int pq{0};

    if (p > q){
        pq = (p*(p-1))/2 + q;
    }
    else{
        pq = (q*(q-1))/2 + p;
    }

    oneBodyEntries[pq-1]=newEntry;
}

//---------------------------------------------------------------------------------------------------//

void Hamiltonian::setMatrixElement(int p, int q, int r, int s, double newEntry){
    // <pq|rs> a+_p a+_q a_s a_r

    int pr{0},qs{0};
    int pqrs{0};

    if (p > r){
        pr = (p*(p-1))/2 + r;
    }
    else{
        pr = (r*(r-1))/2 + p;
    }

    if (q > s){
        qs = (q*(q-1))/2 + s;
    }
    else{
        qs = (s*(s-1))/2 + q;
    }

    if (pr > qs){
        pqrs = (pr*(pr-1))/2 + qs;
    }
    else{
        pqrs = (qs*(qs-1))/2 + pr;
    }


    twoBodyEntries[pqrs-1]=newEntry;
}

//---------------------------------------------------------------------------------------------------//

double Hamiltonian::getMatrixElement(int p, int q) const {
    // <p|h|q> a+_p a_q
    int pq{0};

    if (p > q){
        pq = (p*(p-1))/2 + q;
    }
    else{
        pq = (q*(q-1))/2 + p;
    }

    double h_pq = 0.0;
    h_pq = oneBodyEntries[pq-1];

    return h_pq;
}

//---------------------------------------------------------------------------------------------------//

double Hamiltonian::getMatrixElement(int p, int q, int r, int s) const {

    // <pq|rs> a+_p a+_q a_s a_r

    int pr{0},qs{0};
    int pqrs{0};

    if (p > r){
        pr = (p*(p-1))/2 + r;
    }
    else{
        pr = (r*(r-1))/2 + p;
    }

    if (q > s){
        qs = (q*(q-1))/2 + s;
    }
    else{
        qs = (s*(s-1))/2 + q;
    }

    if (pr > qs){
        pqrs = (pr*(pr-1))/2 + qs;
    }
    else{
        pqrs = (qs*(qs-1))/2 + pr;
    }


    double u_pqrs = 0.0; 
    u_pqrs = twoBodyEntries[pqrs-1];

    return u_pqrs;
}

//---------------------------------------------------------------------------------------------------//


double Hamiltonian::operator()(detType const &alpha, detType const &beta) const{
    // get the Hamiltonian matrix element H_{alpha,beta} between two 
    // configuration alpha and beta
//  std::cout << "size d=" << d << " alpha size=" << alpha.size() << " beta size" << beta.size() << std::endl;
    if(static_cast<int>(alpha.size())!=d || d!=static_cast<int>(beta.size())){
      if(alpha.size()==beta.size()){
      	throw SizeMismatchError(d,alpha.size());
      }
      else{
      	throw SizeMismatchError(alpha.size(),beta.size());
      }
    }

    // first get the spin orbitals involved in the excitation leading
    // from alpha to beta (holes in alpha and particles in beta) 
    // and the parity between the two configurations which is the sign 
    // determined by the number of 
    std::vector<int> excitations;
    std::vector<int> holes;
    std::vector<int> same;
    double diff{0.0};
    bool parity{false};
    double paritysign{1.0};
    int nperm{0};


    // get the differences in occupied spin obitals
    // determine the parity between the configurations
    // this is equal to (-1)^g where g is the number of occupied spin 
    // orbitals between the positions of the spin orbitals involved 
    // in the excitation
    // spin orbitals are stored in range 1,...,M
    for (int i=0;i<d;++i){
        diff=static_cast<int>(alpha[i])-static_cast<int>(beta[i]);
        if (diff>0){
            // holes in alpha
            holes.push_back(i+1);
        }
        if (diff<0){
            // particles in beta
            excitations.push_back(i+1);
        }
        else{
            // these spin orbitals are present in both
            if (diff==0 && alpha[i]){
                same.push_back(i+1);
            }
        }
    }
    // holes and excitations are already ordered, i.e. 
    // holes[1] > holes[0] and excitations[1] > excitations[0]
    if (holes.size()!=excitations.size() || holes.size()>2){
        return 0.0;
    }

    // need to ensure that holes and particles of one excitation share the 
    // same spin for a double excitation
    if (holes.size()==2){
        if ((holes[0]%2)!=(excitations[0]%2)){
            // swap particles so that spins match
            int tmp{0};
            tmp = excitations[0];
            excitations[0] = excitations[1];
            excitations[1] = tmp;
        }
    }

    // get the parity which depends on the distance between the holes and 
    // particles
    // for a double excitation if a crossing of the two excitations occurs 
    // an additional permutation is needed
    nperm = 0;
    detType ex_alpha=alpha;
    if (holes.size()!=0){
        for (size_t i=0; i<holes.size(); ++i){
            nperm += getFermiSign(ex_alpha,(holes[i]-1),(excitations[i]-1));
            // for double excitations: if there is a crossing of the excitation this 
            // complicates everything -> best to simulate excitation
            annihilate(ex_alpha,(holes[i]-1));
            create(ex_alpha,(excitations[i]-1));
        }
    }
    parity = false;
    paritysign = 1.0;
    if ((nperm%2)==1){
        parity = true;
        paritysign = -1.0;
    }


    if (holes.size()==0){
        // diagonal Hamiltonian matrix element, i.e. alpha = beta
        double diagonalTerm{0.0};
        int indi{0},indj{0};
        // since the configurations are the same the parity is +1
        for(size_t i=0;i<same.size();++i){
            // \sum_i <i|h|i>
            indi = this->getId(same[i]);
            diagonalTerm += this->getMatrixElement(indi,indi);
        }
        int idn{0},idx{0};
        for(size_t i=0;(i<same.size()-1);++i){
            for(size_t j=(i+1);j<same.size();++j){
                // \sum_j>i <ij|ij> - <ij|ji>
                // need to ensure that alpha(j) > alpha(i) which 
                // is not guaranteed
                indi = this->getId(same[i]);
                indj = this->getId(same[j]);
                idx = std::max(indi,indj);
                idn = std::min(indi,indj);
                // Coulomb term <ij|ij>
                diagonalTerm += this->getMatrixElement(idn,idx,idn,idx);
                // exchange term <ij|ji>
                // only non-zero contribution when orbitals share the 
                // same spin
                if ((same[i]%2)==(same[j]%2)){
                    diagonalTerm -= this->getMatrixElement(idn,idx,idx,idn);
                }
            }
        }
        
        return diagonalTerm;
    }
    else if (holes.size()==1){
        // single excitations
        double singleTerm{0.0};
        int indi{0},indj{0},inda{0};


        // <i|h|a>
        if ((holes[0]%2)==(excitations[0]%2)){
            indi = this->getId(holes[0]);
            inda = this->getId(excitations[0]);
            singleTerm += this->getMatrixElement(indi,inda);
        }

        // \sum_j <ij|aj> - <ij|ja>
        for (size_t i=0;i<same.size();++i){
            indi = this->getId(holes[0]);
            indj = this->getId(same[i]);
            inda = this->getId(excitations[0]);
            // coulomb term <ij|aj>
            if ((holes[0]%2)==(excitations[0]%2)){
                singleTerm += this->getMatrixElement(indi,indj,inda,indj);
            }
            // exchange term <ij|ja>
            if (((holes[0]%2)==(excitations[0]%2))&&((holes[0]%2)==(same[i]%2))){
                singleTerm -= this->getMatrixElement(indi,indj,indj,inda);
            }
        }
        // self terms cancel 

        // parity
        if (parity){
            singleTerm *= paritysign;
        }

        return singleTerm;

    }
    else if (holes.size()==2){
        // double excitation
        double doubleTerm{0.0};
        int indi{0},indj{0},inda{0},indb{0};

        // no contribution from one-particle integrals

        // for excitation ij->ab
        // <ij|ab> - <ij|ba>
        // only non-zero contributions are those for which the spins of 
        // the orbitals are the same 
        // coulomb term <ij|ab>
        if (((holes[0]%2)==(excitations[0]%2))&&((holes[1]%2)==(excitations[1]%2))){
            indi = this->getId(holes[0]);
            indj = this->getId(holes[1]);
            inda = this->getId(excitations[0]);
            indb = this->getId(excitations[1]);
            doubleTerm = this->getMatrixElement(indi,indj,inda,indb);
        }
        else{
            doubleTerm = 0.0;
        }
        // only non-zero contributions are those for which the spins of 
        // the orbitals are the same 
        // exchange term <ij|ba>
        if (((holes[0]%2)==(excitations[1]%2))&&((holes[1]%2)==(excitations[0]%2))){
            indi = this->getId(holes[0]);
            indj = this->getId(holes[1]);
            inda = this->getId(excitations[0]);
            indb = this->getId(excitations[1]);
            doubleTerm -= this->getMatrixElement(indi,indj,indb,inda);
        }
        
        // parity
        if (parity){
            doubleTerm *= paritysign;
        }

        return doubleTerm;
    }
    else{
        // more than double excitation -> 0
        double Term{0.0};

        Term = 0.0;

        return Term;
    }
}

//---------------------------------------------------------------------------------------------------//

//for testing purposes - prints hamiltonian in N particle sector

void Hamiltonian::printMatrix(int N){
  std::cout<<"Printing Hamiltonian\n";
  std::vector<int> alpha, beta;
}

