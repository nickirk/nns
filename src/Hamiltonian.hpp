#ifndef HAMILTONIAN_CLASS
#define HAMILTONIAN_CLASS

#include "Determinant.hpp"

class Hamiltonian{
 public:
  Hamiltonian(){}
  //dimension is the dimension of the single-particle hilbert space
  explicit Hamiltonian(int dimension):d(dimension),oneBodyEntries(std::vector<double>(d*d,0.0)),twoBodyEntries(std::vector<double>(d*d*d*d,0.0)){}
  void setMatrixElement(int r, int s, double newEntry);
  void setMatrixElement(int p, int q, int r, int s, double newEntry);
  double getMatrixElement(detType const &alpha, detType const &beta)const;
  void printMatrix(int N);
 private:
  int d;
  //these are the coefficients of the second quantized hamiltonian
  std::vector<double> oneBodyEntries;
  std::vector<double> twoBodyEntries;
};

detType getRandomCoupledState(detType const &source, double &p);
//basically counts the number of particles between annihilatorIndex and creatorIndex in alpha, yielding the sign of a^\dagger_i a_j (alpha) in the canonical operator ordering
int getFermiSign(detType const &alpha, int annihilatorIndex, int creatorIndex);
Hamiltonian generateHubbard(int dim, double U, double t);

#endif
