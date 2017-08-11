#include <vector>
#include <random>
#include <ModelSys.hpp>
#include <Sampler.hpp>

void sampler::generateList(std::vector<detType > &list) const{
  // this just repeatedly gets new random states and adds them to the list
  list = std::vector<detType >(numStates);
  for(int i=0;i<numStates;++i){
    cDet=getRandomDeterminant(cDet);
    list[i]=cDet;
  }
}

detType getRandomDeterminant(detType const &startingPoint) const{
  std::vector<detType > tempDets(0);
  // first, convert the determinant to an index
  int j = intcast(startingPoint);
  // get the range that corresponds to this index' row
  // here, we need the sorted representation. This is best done in the Hamiltonian class
  int lower = H.lowerPos(j);
  int upper = H.upperPos(j);
  int row,col;
  std::random_device rng;
  double const normalizer = static_cast<double>(rng.max());
  double value,p;
  double K=0.0;
  // compute the normalization sum_j |K_ij|
  for(int i=lower;i<=upper;++i){
    H.sparseAccess(i,row,col,value);
    K+=std::abs(value);
  }
  for(int i=lower;i<=upper;++i){
    p=rng()/normalizer;
    // get all coupled determinants and accept them into the temporary list
    // with probability K_ij/K
    H.sparseAccess(i,row,col,value);
    if(p<value/K){
      // here, we need the reverse of intcast, i.e. the conversion of the index
      // to the determinant. It shall just return lookuptable(i) for a given index i
      // (lookuptable contains the determinants for conversion to int)
      tempDets.push_back(detTypeCast(col));
    }
  }
  // pick a random determinant from the temporary list
  int const chosen=static_cast<int>(rng()/normalizer*tempDets.size());
  return tempDets[chosen];
}
