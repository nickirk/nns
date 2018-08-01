#ifndef TwoBodyHamiltonian_CLASS
#define TwoBodyHamiltonian_CLASS

#include "../utilities/Errors.hpp"
#include "../utilities/TypeDefine.hpp"
#include "Hamiltonian.hpp"

namespace networkVMC{

// generic two-body hamiltonian (no specification if bosonic or fermionic)
// in principle, all two-body Hamiltonians could use this, but for lattice models, dedicated
// implementations are more efficient, so this is mainly for Ab-initio
class TwoBodyHamiltonian: public Hamiltonian{
 public:
  // constructors
  TwoBodyHamiltonian():d(0),donebodyint(0),dtwobodyint(0),spinOrbs(false),coreEnergy(0.0),
  linExactFlag(true),partExactFlag(true){};
  //dimension is the dimension of the single-particle Hilbert space
  explicit TwoBodyHamiltonian(int dimension):
		  Hamiltonian(),d(dimension),donebodyint(((d*(d+1))/2)),dtwobodyint(((donebodyint*(donebodyint+1))/2)),
		  spinOrbs(false),oneBodyEntries(std::vector<double>(donebodyint,0.0)),
		  twoBodyEntries(std::vector<double>(dtwobodyint,0.0)),coreEnergy(0.0),partExactFlag(true),linExactFlag(true){}
  //explicit Hamiltonian(int dimension):d(dimension),oneBodyEntries(std::vector<double>(d*d,0.0)),twoBodyEntries(std::vector<double>(d*d*d*d,0.0)){}
  // set the one- and two-body integrals of the Hamiltonian operator
  void setMatrixElement(int p, int q, double newEntry);
  void setMatrixElement(int p, int q, int r, int s, double newEntry);
  // set the core energy of the Hamiltonian operator
  void setMatrixElement(double newEntry);
  // get the one- and two-body integrals of the Hamiltonian operator
  double getMatrixElement(int p, int q) const;
  double getMatrixElement(int p, int q, int r, int s) const;
  // get the core energy of the Hamiltonian operator
  double getMatrixElement() const;
  // convert a spin orbital to a spatial orbitals index
  int getId(int i) const;
  // get the number of orbitals
  int getNumOrbs() const {return d;}
  void initMatrixStorage(bool bspin_orbs);
  // get the TwoBodyHamiltonian matrix element between two configurations alpha
  // and beta
  // can throw an InvalidDeterminantError or SizeMismatchError if alpha/beta do not make sense
  virtual double operator()(detType const &alpha, detType const &beta) const;

  void printMatrix(int N);
  // this generates all states coupled to source
  // note that - in contrast to excitation generation - this is really
  // a property of the TwoBodyHamiltonian
  virtual std::vector<detType> getCoupledStates(detType const &source) const = 0;
  virtual ~TwoBodyHamiltonian(){};

  bool linExact() const {return linExactFlag;}
  bool partExact() const {return partExactFlag;}

  // this is for setting defaults: choose a default excitgen depending
  // on the type of H
  virtual HType type() const {return Constant;}
 private:
  // dimenstion of single-particle Hilbert space
  int d;
  // dimension of 1-e integrals
  int donebodyint;
  // dimenstion of 2-e integrals
  int dtwobodyint;
  // are the integrals stored in spin or spatial orbitals
  bool spinOrbs;
  //these are the coefficients of the second quantized TwoBodyHamiltonian
  // 1-e integrals <i|h|j>
  std::vector<double> oneBodyEntries;
  // 2-e integrals <ij|kl>
  std::vector<double> twoBodyEntries;
  // core energy
  double coreEnergy;
  int twoBodyIndex(int s, int r, int q, int p) const;
  int oneBodyIndex(int s, int r) const;

protected:
  bool linExactFlag;
  bool partExactFlag;
  // the commutation relation of the underlying creation/annihilation operators goes in here
  virtual int getFermiSign(detType const &alpha, int annihilatorIndex, int creatorIndex) const =0;
};

detType getRandomCoupledState(detType const &source, double &p);
std::vector<detType> getCoupledStates(detType const &source);
//basically counts the number of particles between annihilatorIndex and creatorIndex in alpha, yielding the sign of a^\dagger_i a_j (alpha) in the canonical operator ordering
int getFermiSign(detType const &alpha, int annihilatorIndex, int creatorIndex);

}

#endif
