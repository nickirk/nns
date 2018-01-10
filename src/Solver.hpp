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
	Eigen::VectorXd update(Eigen::VectorXd const &w, Eigen::VectorXd const &force) const;
	Eigen::VectorXd update(Eigen::VectorXd const &w, Eigen::VectorXd const &force,
						   Eigen::VectorXd const &ci, Eigen::MatrixXd const &dcdw) const;
private:
	double gamma;
};

#endif /* SRC_SOLVER_HPP_ */
