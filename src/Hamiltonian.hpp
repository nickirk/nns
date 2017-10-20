#ifndef HAMILTONIAN_CLASS
#define HAMILTONIAN_CLASS

#include "Determinant.hpp"

class Hamiltonian{
 public:
  Hamiltonian(){}
  //dimension is the dimension of the single-particle Hilbert space
  explicit Hamiltonian(int dimension):d(2*dimension),oneBodyEntries(std::vector<double>(d*d,0.0)),twoBodyEntries(std::vector<double>(d*d*d*d,0.0)){}
  void setMatrixElement(int r, int s, double newEntry);
  void setMatrixElement(int p, int q, int r, int s, double newEntry);
  double operator()(detType const &alpha, detType const &beta) const{return getMatrixElement(alpha, beta);}
  void printMatrix(int N);
 private:
  int d;
  //these are the coefficients of the second quantized hamiltonian
  std::vector<double> oneBodyEntries;
  std::vector<double> twoBodyEntries;
  double getMatrixElement(detType const &alpha, detType const &beta)const;
};

detType getRandomCoupledState(detType const &source, double &p);
//basically counts the number of particles between annihilatorIndex and creatorIndex in alpha, yielding the sign of a^\dagger_i a_j (alpha) in the canonical operator ordering
int getFermiSign(detType const &alpha, int annihilatorIndex, int creatorIndex);
Hamiltonian generateHubbard(int dim, double U, double t);

#endif
