/*
 * Solver.hpp
 *
 *  Created on: Jan 10, 2018
 *      Author: guther
 */

#ifndef SRC_SOLVER_HPP_
#define SRC_SOLVER_HPP_

#include <Eigen/Dense>

class Solver {
public:
	explicit Solver(double eta);
	~Solver();
	void setGamma(double eta){gamma = eta;}
	void update(Eigen::VectorXd  &w, Eigen::VectorXd const &force) const;
	void update(Eigen::VectorXd  &w, Eigen::VectorXd const &force,
						   Eigen::VectorXcd const &ci, Eigen::MatrixXcd const &dcdw, int const &iteration) const;
private:
	double gamma;
};

#endif /* SRC_SOLVER_HPP_ */
