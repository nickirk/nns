/*
 * EnergyEs.hpp
 *
 *  Created on: Jun 26, 2018
 *      Author: guther
 */

#ifndef SRC_COSTFUNCTIONS_ENERGYES_HPP_
#define SRC_COSTFUNCTIONS_ENERGYES_HPP_

#include "CostFunction.hpp"
#include "EnergyCFBaseClass.hpp"
#include "../utilities/DeepCpyUniquePtr.hpp"
#include <Eigen/Dense>

// set the default arguments
#include "EnergyEsForward.hpp"

namespace networkVMC {

class Hamiltonian;

/**
 * \class EnergyEs
 * \brief Energy-based CostFunction using a stochastic evaluation of the expectation value
 *
 * \tparam F Type of the parameters to optimize
 * \tparam coeffType Type of the vector coefficients of the input vector
 * The implemented CostFunction uses the expectation value of the energy expectation value
 * with respect to a stochastic sampling of the basis vectors
 */
template <typename F, typename coeffType>
class EnergyEs: public CostFunction<F, coeffType> {
  public:
	/// Type of the derivative
	using T=Eigen::Matrix<F, Eigen::Dynamic, 1>;
	/**
	 * \param H_ Hamiltonian used for obtaining the energy
	 * \param numCons_ Number of matrix elements taken into account per basis vector
	 */
	EnergyEs(Hamiltonian const &H_, int numCons_=20);
	virtual ~EnergyEs();
	virtual void setUpCF(SamplerType const &sT);

// The operations are actually performed by another EnergyCF (worker)
	T nabla(State<coeffType> const &input) const {return worker->nabla(input);}
	coeffType calc(State<coeffType> const &input) const {return worker->calc(input);};

	/**
	 * \brief Auxiliary function to get the norm
	 * \return The norm of the last input state
	 */
	double getNormalizer() const {return worker->getNormalizer();}

	// Allow for polymorphic copy
	virtual EnergyEs* clone() const {return new EnergyEs(*this);}

	// the energy estimators do need connections
	virtual int connectionsRequired() const {return numCons;}
private:
	/// Hamiltonian used for obtaining the energy
	Hamiltonian const& H;
	// this is not a stand-alone CF, the work is done by another,
	// owned CF
	/// Owned CostFunction doing the work
	DeepCpyUniquePtr<EnergyCFBaseClass<F, coeffType>> worker;
	/// number of matrix elements to be taken into account per basis vector
	int numCons;
};

} /* namespace networkVMC */

#endif /* SRC_COSTFUNCTIONS_ENERGYES_HPP_ */
