#include <vector>
#include <random>
#include <ModelSys.hpp>

// type of the determinants
typedef std::vector<bool> detType;

class sampler{
public:
  sampler(Hamiltonian const &H_, int numStates_, detType HF):H(H_),numStates(numStates_),
							     cDet(HF){}
  // two functionalities: get a random coupled determinant and get an array of 
  // random coupled determinants
  void generateList(std::vector<detType > &list)const;
  detType getRandomDeterminant(detType const &startingPoint) const;
  // for ab-initio: introduce an overload of generateList for ab-initio hamiltonians
private:
  // Hamiltonian
  Hamiltonian &H;
  int numStates;
  // this tracks the last determinant
  mutable detType cDet;
};
