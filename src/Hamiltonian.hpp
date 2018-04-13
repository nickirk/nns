#ifndef HAMILTONIAN_CLASS
#define HAMILTONIAN_CLASS

#include "Determinant.hpp"
#include "./utilities/Errors.hpp"

class Hamiltonian{
 public:
  // constructors
  Hamiltonian():d(0){}
  //dimension is the dimension of the single-particle Hilbert space
  explicit Hamiltonian(int dimension):d(dimension),donebodyint(((d*(d+1))/2)),dtwobodyint(((donebodyint*(donebodyint+1))/2)),oneBodyEntries(std::vector<double>(donebodyint,0.0)),twoBodyEntries(std::vector<double>(dtwobodyint,0.0)){}
  //explicit Hamiltonian(int dimension):d(dimension),oneBodyEntries(std::vector<double>(d*d,0.0)),twoBodyEntries(std::vector<double>(d*d*d*d,0.0)){}
  // set the one- and two-body integrals of the Hamiltonian operator
  void setMatrixElement(int p, int q, double newEntry);
  void setMatrixElement(int p, int q, int r, int s, double newEntry);
  // get the one- and two-body integrals of the Hamiltonian operator
  double getMatrixElement(int p, int q) const;
  double getMatrixElement(int p, int q, int r, int s) const;
  // convert a spin orbital to a spatial orbitals index
  int getId(int i) const;
  void initMatrixStorage(bool bspin_orbs);
  // get the Hamiltonian matrix element between two configurations alpha 
  // and beta
  double operator()(detType const &alpha, detType const &beta) const;
  void printMatrix(int N);
  virtual int getFermiSign(detType const &alpha, int annihilatorIndex, int creatorIndex) const =0;
  virtual detType getRandomCoupledState(detType const &source, double &p) const = 0;
  virtual std::vector<detType> getCoupledStates(detType const &source) const = 0;
  virtual ~Hamiltonian(){};
 private:
  // dimenstion of single-particle Hilbert space
  int d;
  // dimension of 1-e integrals
  int donebodyint;
  // dimenstion of 2-e integrals
  int dtwobodyint;
  // are the integrals stored in spin or spatial orbitals
  bool spin_orbs;
  //these are the coefficients of the second quantized hamiltonian
  // 1-e integrals <i|h|j>
  std::vector<double> oneBodyEntries;
  // 2-e integrals <ij|kl>
  std::vector<double> twoBodyEntries;
  int twoBodyIndex(int s, int r, int q, int p) const{
		return s + r * d + q * d * d + p * d * d * d;
  }

  int oneBodyIndex(int s, int r) const{
		return s + d * r;
  }
};

detType getRandomCoupledState(detType const &source, double &p);
std::vector<detType> getCoupledStates(detType const &source);
//basically counts the number of particles between annihilatorIndex and creatorIndex in alpha, yielding the sign of a^\dagger_i a_j (alpha) in the canonical operator ordering
int getFermiSign(detType const &alpha, int annihilatorIndex, int creatorIndex);

#endif
