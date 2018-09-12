#ifndef TwoBodyHamiltonian_CLASS
#define TwoBodyHamiltonian_CLASS

#include "../utilities/Errors.hpp"
#include "../utilities/TypeDefine.hpp"
#include "Hamiltonian.hpp"

namespace networkVMC{

/**
* \class TwoBodyHamiltonian
* \brief generic two-body hamiltonian (no specification if bosonic or fermionic)
*/

// in principle, all two-body Hamiltonians could use this, but for lattice models, dedicated
// implementations are more efficient, so this is mainly for Ab-initio
class TwoBodyHamiltonian: public Hamiltonian{
   public:
  /// Sets the dimension of the local Hilbert space to 0
  TwoBodyHamiltonian():d(0),donebodyint(0),dtwobodyint(0),spinOrbs(false),coreEnergy(0.0),
  linExactFlag(true),partExactFlag(true){};
  /// \param[in] dimension The dimension of the single-particle Hilbert space
  explicit TwoBodyHamiltonian(int dimension):
		  Hamiltonian(),d(dimension),donebodyint(((d*(d+1))/2)),dtwobodyint(((donebodyint*(donebodyint+1))/2)),
		  spinOrbs(false),oneBodyEntries(std::vector<double>(donebodyint,0.0)),
		  twoBodyEntries(std::vector<double>(dtwobodyint,0.0)),coreEnergy(0.0),linExactFlag(true),partExactFlag(true){}
  /**
   * \brief Set the 1-e integrals of the Hamiltonian operator
   * \param[in] p,q Indices of the involved orbitals
   * \param[in] newEntry New value of the matrix element
   */
  void setMatrixElement(int p, int q, double newEntry);

  /**
   * \brief Set the 2-e integrals of the Hamiltonian operator
   * \param[in] p,q,r,s Indices of the involved orbitals
   * \param[in] newEntry New value of the matrix element
   */
  void setMatrixElement(int p, int q, int r, int s, double newEntry);
  /**
   * \brief Set the core energy of the Hamiltonian operator
   * \param[in] newEntry New value of the matrix element
   */
  void setMatrixElement(double newEntry);
  // get the one- and two-body integrals of the Hamiltonian operator

  /**
   * \brief Get the matrix element of a single-excitation p->q
   * \param[in] p,q Indices of orbitals involved
   * \return Matrix element (p|q)
   */
  double getMatrixElement(int p, int q) const;
  /**
   * \brief Get the matrix element of a single-excitation pq->rs
   * \param[in] p,q,r,s Indices of orbitals involved
   * \return Matrix element (pr|qs)
   */
  double getMatrixElement(int p, int q, int r, int s) const;
  /**
   * \brief Get the core energy of the Hamiltonian operator
   * \return Core matrix element
   */
  double getMatrixElement() const;
  /**
   * \brief Convert a spin orbital to a spatial orbitals index
   * \param[in] i Spin orbital index
   * \return Index of the orbital as used for matrix elements (spatial or spin, depending on Hamiltonian)
   */
  int getId(int i) const;
  /// \return The number of orbitals
  int getNumOrbs() const {return d;}
  /// \brief Initialize the buffers for the 1-e and 2-e integrals
  void initMatrixStorage(bool bspin_orbs);
  // get the TwoBodyHamiltonian matrix element between two configurations alpha and beta
  // can throw an InvalidDeterminantError or SizeMismatchError if alpha/beta do not make sense
  virtual double operator()(detType const &alpha, detType const &beta) const;

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
  /// dimenstion of single-particle Hilbert space
  int d;
  /// dimension of 1-e integrals
  int donebodyint;
  /// dimenstion of 2-e integrals
  int dtwobodyint;
  /// are the integrals stored in spin or spatial orbitals
  bool spinOrbs;
  //these are the coefficients of the second quantized TwoBodyHamiltonian
  /// 1-e integrals <i|h|j>
  std::vector<double> oneBodyEntries;
  /// 2-e integrals <ij|kl>
  std::vector<double> twoBodyEntries;
  /// core energy
  double coreEnergy;

  /**
   * \brief Indexing function for 2-e integrals
   * \param[in] s,r,q,p Indices of the integral (qr|sp)
   * \return Index of <sr|qp> in the contiguous array of integrals
   */
  int twoBodyIndex(int s, int r, int q, int p) const;
  /**
   * \brief Indexing function for 1-e integrals
   * \param[in] s,r Indices of the integral (s|r)
   * \return Index of (s|r) in the contiguous array of integrals
   */
  int oneBodyIndex(int s, int r) const;

protected:
  bool linExactFlag;
  bool partExactFlag;
  /**
   *  \brief Commutation relation of the underlying creation/annihilation operators
   *  \param[in] alpha Basis vector to apply operators to
   *  \param[in] annihilatorIndex Index of annihilation operator to apply
   *  \param[in] creatorIndex Index of creation operator to apply
   *  \return Sign of the resulting basis vector relative to the original one
   *
   *  When applying an excitation operator to a basis vector alpha, an additional sign pops up
   *  depending on the commutation relation of the creation/annihilation operators. This function computes
   *  said sign.
   */
  virtual int getFermiSign(detType const &alpha, int annihilatorIndex, int creatorIndex) const =0;
};

}

#endif
