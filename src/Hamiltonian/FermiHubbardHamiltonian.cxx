/*
 * FermiHubbardHamiltonian.cxx
 *
 *  Created on: Jun 20, 2018
 *      Author: guther
 */

#include "FermiHubbardHamiltonian.hpp"
#include "../HilbertSpace/Determinant.hpp"
#include "ExcitationGenerators/RSHubbardExcitgen.hpp"

namespace networkVMC {

FermiHubbardHamiltonian::~FermiHubbardHamiltonian() {
}

std::vector<detType> FermiHubbardHamiltonian::getCoupledStates(detType const &source) const{
  int const d=source.size();
  std::vector<detType> coupledList;
  std::vector<int> spawnLeft, spawnRight;
  // this creates a list of all possible hoppings
  getRSHubSpawnLists(source, spawnLeft, spawnRight);

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
  return coupledList;
}

//---------------------------------------------------------------------------------------------------//

FermiHubbardHamiltonian generateFermiHubbard(int dim, double U, double t){
  FermiHubbardHamiltonian H(dim);
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

} /* namespace networkVMC */
