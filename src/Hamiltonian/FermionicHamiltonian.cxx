/*
 * fermionicHamiltonian.cxx
 *
 *  Created on: Nov 15, 2017
 *      Author: guther
 */

#include <random>
#include "FermionicHamiltonian.hpp"
#include "../HilbertSpace/Determinant.hpp"

namespace networkVMC{

FermionicHamiltonian::~FermionicHamiltonian() {
}

int FermionicHamiltonian::getFermiSign(detType const &alpha, int annihilatorIndex, int creatorIndex) const{
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

//int FermionicHamiltonian::getFermiSign(detType const &alpha, int annihilatorIndex, int creatorIndex) const{
//  //construct sign for conversion canonical shape (basisState 1(up),1(down),...,L(up),L(down)) to relative shape (a_exc^\dagger a_holes source)
//  int fermiSign{1};
//  int start{0}, end{0}, offset{0};
//  if(annihilatorIndex>creatorIndex){
//    start=creatorIndex;
//    end=annihilatorIndex;
//    offset=0;
//  }
//  else{
//    start=annihilatorIndex;
//    end=creatorIndex;
//    offset=1;
//  }
//  for(int k=start+offset;k<end;++k){
//    if(alpha[k]){
//      fermiSign*=-1;
//    }
//  }
//  return fermiSign;
//}

}

