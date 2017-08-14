#ifndef Sampler_DEFINED
#define Sampler_DEFINED

#include <detType.hpp>
#include <Hamiltonian.hpp>

class sampler{
public:
  sampler(Hamiltonian const &H_, detType const &HF, Basis const &fullBasis_, int numStates_):H(H_),fullBasis(fullBasis_),
							     numStates(numStates_),cDet(HF){}
  // two functionalities: get a random coupled determinant and get an array of 
  // random coupled determinants
  void generateList(std::vector<detType > &list)const;
  detType getRandomDeterminant(detType const &startingPoint) const;
  // for ab-initio: introduce an overload of generateList for ab-initio hamiltonians
private:
  // Hamiltonian
  Hamiltonian &H;
  Basis &fullBasis;
  int numStates;
  // this tracks the last determinant
  mutable detType cDet;
};

inline double dblAbs(double x){if(x>0) return x;return -x;}

#endif
