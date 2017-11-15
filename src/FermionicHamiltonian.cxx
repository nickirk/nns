/*
 * fermionicHamiltonian.cxx
 *
 *  Created on: Nov 15, 2017
 *      Author: guther
 */

#include "FermionicHamiltonian.hpp"

FermionicHamiltonian::~FermionicHamiltonian() {
}

int FermionicHamiltonian::getFermiSign(detType const &alpha, int annihilatorIndex, int creatorIndex) const{
  //construct sign for conversion canonical shape (basisState 1(up),1(down),...,L(up),L(down)) to relative shape (a_exc^\dagger a_holes source)
  int fermiSign{1};
  int start{0}, end{0}, offset{0};
  if(annihilatorIndex>creatorIndex){
    start=creatorIndex;
    end=annihilatorIndex;
    offset=0;
  }
  else{
    start=annihilatorIndex;
    end=creatorIndex;
    offset=1;
  }
  for(int k=start+offset;k<end;++k){
    if(alpha[k]){
      fermiSign*=-1;
    }
  }
  return fermiSign;
}

FermionicHamiltonian generateFermiHubbard(int dim, double U, double t){
  FermionicHamiltonian H(dim);
  for(int i=0;i<dim-1;i+=2){
    //Hubbard interaction
    H.setMatrixElement(i,i+1,i,i+1,-U);
  }
  for(int i=0;i<dim-2;++i){
    //hopping, excluding the pbc term
    H.setMatrixElement(i,i+2,-t);
    H.setMatrixElement(i+2,i,-t);
  }
  for(int sigma=0;sigma<2;++sigma){
    //boundary term for pbc
    H.setMatrixElement(dim-1-sigma,1-sigma,-t);
    H.setMatrixElement(1-sigma,dim-1-sigma,-t);
  }
  return H;
}

