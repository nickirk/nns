/*
 * Solver.hpp
 *
 *  Created on: Jan 10, 2018
 *      Author: guther
 */

#ifndef SRC_SOLVER_HPP_
#define SRC_SOLVER_HPP_

#include <Eigen/Dense>

// Class doing the minimization in a force-field given by force by different schemes
class Solver {
public:
	explicit Solver(double eta);
	~Solver();
// gamma is a prefactor used in the iterations
	void setGamma(double eta){gamma = eta;}
// We may either only supply the force
	void update(Eigen::VectorXd  &w, Eigen::VectorXd const &force) const;
// Or also the second derivatives
	void update(Eigen::VectorXd  &w, Eigen::VectorXd const &force,
						   Eigen::VectorXcd const &ci, Eigen::MatrixXcd const &dcdw, int const &iteration) const;
private:
	double gamma;
};

#endif /* SRC_SOLVER_HPP_ */
