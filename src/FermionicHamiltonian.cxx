/*
 * fermionicHamiltonian.cxx
 *
 *  Created on: Nov 15, 2017
 *      Author: guther
 */

#include "FermionicHamiltonian.hpp"
#include <random>

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

FermionicHamiltonian generateFermiHubbard(int dim, double U, double t){
  FermionicHamiltonian H(dim);
  H.initMatrixStorage(true);

  for(int i=2;i<dim+1;i+=2){
      //Hubbard interaction
      H.setMatrixElement(i,i-1,i,i-1,U);
      H.setMatrixElement(i-1,i,i-1,i,U);
      //hopping term, excluding pbc
      H.setMatrixElement(i,i+2,t);
      H.setMatrixElement(i-1,i+1,t);
  }
  for(int sigma=0;sigma<2;++sigma){
      //boundary term for pbc
      H.setMatrixElement(dim-sigma,2-sigma,t);
  }
  
  return H;
}

std::vector<detType> FermionicHamiltonian::getCoupledStates(detType const &source) const{
  int const d=source.size();
  std::vector<detType> coupledList;
  std::vector<int> spawnLeft, spawnRight;
  for(int i=0;i<d;++i){
    //search all sites from which one can hop to the left/right
    if(source[i]){
      if(source[(i+2)%d]==0){
	spawnRight.push_back(i);
      }
      if(source[(i-2+d)%d]==0){
	spawnLeft.push_back(i);
      }
    }
  }
  detType targetTmp=source;
  for (unsigned int i=0; i<spawnLeft.size(); ++i){
    targetTmp=source;
    annihilate(targetTmp,spawnLeft[i]);
    create(targetTmp,(spawnLeft[i]-2+d)%d);
    coupledList.push_back(targetTmp);
  }
  for (unsigned int i=0; i<spawnRight.size(); ++i){
    targetTmp=source;
    annihilate(targetTmp,spawnRight[i]);
    create(targetTmp,(spawnRight[i]+2)%d);
    coupledList.push_back(targetTmp);
  }
  int const numSpawn=spawnLeft.size()+spawnRight.size();
  return coupledList;
}

detType FermionicHamiltonian::getRandomCoupledState(detType const &source,  double &pGet) const{
  std::random_device rng;
  double const normalization=static_cast<double>(rng.max());
  int const d=source.size();
  detType target=source;
  //for hubbard
  //hubbard one-particle basis: 1(up), 1(down), 2(up), ..., L(up), L(down)
  std::vector<int> spawnLeft, spawnRight;
  for(int i=0;i<d;++i){
    //search all sites from which one can hop to the left/right
    if(source[i]){
      if(source[(i+2)%d]==0){
	spawnRight.push_back(i);
      }
      if(source[(i-2+d)%d]==0){
	spawnLeft.push_back(i);
      }
    }
  }
  int const numSpawn=spawnLeft.size()+spawnRight.size();
  int p{static_cast<int>(rng()/normalization*(numSpawn+1))};
  //pick one of those sites at random (including direction)
  pGet=1.0/static_cast<double>(numSpawn+1);
  if(p>=numSpawn){
    return source;
  }
  int const offset=spawnLeft.size();
  int newPos{0};
  if(p>=offset){
    newPos=spawnRight[p-offset];
    annihilate(target,newPos);
    create(target,(newPos+2)%d);
  }
  else{
    newPos=spawnLeft[p];
    annihilate(target,newPos);
    create(target,(newPos-2+d)%d);
  }
  //for ab-initio
  /*
  std::vector<int> holes, exc;
  for(int i=0;i<d;++i){
    if(source[i])
      exc.push_back(i);
    else
      holes.push_back(i);
  }
  int numHole=holes.size();
  int numExc=exc.size();
  if(numHole<2 || numExc<2){
    pGet=1;
    return source;
  }
  int p=static_cast<int>(rng()/normalization*numHole);
  int q=static_cast<int>(rng()/normalization*(numHole-1));
  if(q>=p){
    //ensure q!=p
    ++q;
  }
  int r=static_cast<int>(rng()/normalization*numExc);
  int s=static_cast<int>(rng()/normalization*(numExc-1));
  if(r>=s){
    //ensure r!=s
    ++r;
  }
  target.addParticle(holes[p]);
  target.addParticle(holes[q]);
  target.removeParticle(exc[r]);
  target.removeParticle(exc[s]);
  //simplest algorithm, assuming global couplings
  pGet=4.0/(numExc*numHole*(numHole-1)*(numExc-1));
  */
  return target;
}

