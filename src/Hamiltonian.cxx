#include "Hamiltonian.hpp"
#include <random>
#include <cmath>
#include <vector>
#include <iostream>

void Hamiltonian::setMatrixElement(int r, int s, double newEntry){
  //r is creation inedx, s is annihilation index
  oneBodyEntries[s+d*r]=newEntry;
}

//---------------------------------------------------------------------------------------------------//

void Hamiltonian::setMatrixElement(int p, int q, int r, int s, double newEntry){
  //twoBodyEntries are meant to be antisymmetrized
  //p,q are creation indices, r,s are annihilation indices
  //for input, order of p,q and r,s does not matter, all entries are generated anyway (might be optimized)
  //default input: r>s, p>q
  if(p==q || r==s){
    newEntry=0.0;
  }
  twoBodyEntries[s+r*d+q*d*d+p*d*d*d]=newEntry;
  twoBodyEntries[r+s*d+q*d*d+p*d*d*d]=-newEntry;
  twoBodyEntries[s+r*d+p*d*d+q*d*d*d]=-newEntry;
  twoBodyEntries[r+s*d+p*d*d+q*d*d*d]=newEntry;
}

//---------------------------------------------------------------------------------------------------//

double Hamiltonian::getMatrixElement(detType const &alpha, detType const &beta) const{
  if(static_cast<int>(alpha.size())!=d || d!=static_cast<int>(beta.size())){
    return 0.0;
  }
  std::vector<int> excitations;
  std::vector<int> holes;
  std::vector<int> same;
  double diff{0.0};
  for(int i=0;i<d;++i){
    diff=static_cast<int>(alpha[i])-static_cast<int>(beta[i]);
    if(diff>0){
      excitations.push_back(i);
      }
    if(diff<0){
      holes.push_back(i);
    }
    else{
      if(diff==0 && alpha[i]){
	same.push_back(i);
      }
    }
  }
  if(holes.size()!=excitations.size() || holes.size()>2){
    return 0.0;
  }
  if(holes.size()==0){
    double diagonalTerm{0.0};
    //there is no sign in the diagonal term since all appearing operators are bosonic (a_i^\dagger a_i)
    for(unsigned int i=0;i<same.size();++i){
      //all contributions if alpha==beta
      diagonalTerm+=oneBodyEntries[same[i]+d*same[i]];
      for(unsigned int j=i;j<same.size();++j){
	//use this ordering to eliminate fermionic sign (this is basically n_i*n_j)
	diagonalTerm+=twoBodyEntries[same[i]+d*same[j]+d*d*same[j]+d*d*d*same[i]];
      }
    }
    return diagonalTerm;
  }
  int fermiSign{-1};
  //Sign for conversion to canonical form
  if(holes.size()==2){
    fermiSign*=getFermiSign(beta,holes[0],excitations[0]);
    detType proxy=beta;
    annihilate(proxy,holes[0]);
    create(proxy,(excitations[0]));
    fermiSign*=getFermiSign(proxy,holes[1],excitations[1]);
    //it is always holes[1]>holes[0] and excitations[1]>excitations[0]
    //take into account fermi sign due to states with indices between holes[0] and holes[1]
    return fermiSign*twoBodyEntries[holes[0]+d*holes[1]+d*d*excitations[0]+d*d*d*excitations[1]];
  }
  double twoBodyTerm{0.0};
  int localSign{0};
  fermiSign=getFermiSign(beta,holes[0],excitations[0]);
  for(int k=0;k<d;++k){
    //sum up all second-order contributions. The correct sign is included in the definition of twoBodyEntries
    if(alpha[k]){
      if((k>holes[0] && k>excitations[0]) || (k<holes[0] && k<excitations[0]))
	localSign=-1;
      else
	localSign=1;
      twoBodyTerm+=localSign*twoBodyEntries[excitations[0]+d*k+d*d*holes[0]+d*d*d*k];
    }
  }
  return fermiSign*(oneBodyEntries[excitations[0]+d*holes[0]]+twoBodyTerm);
}

//---------------------------------------------------------------------------------------------------//

detType getRandomCoupledState(detType const &source, double &pGet){
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

//---------------------------------------------------------------------------------------------------//

//for testing purposes - prints hamiltonian in N particle sector

void Hamiltonian::printMatrix(int N){
  std::cout<<"Printing Hamiltonian\n";
  std::vector<int> alpha, beta;
}

int getFermiSign(detType const &alpha, int annihilatorIndex, int creatorIndex){
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

Hamiltonian generateHubbard(int dim, double U, double t){
  Hamiltonian H(dim);
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
