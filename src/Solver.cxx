/*
 * Solver.cxx
 *
 *  Created on: Jan 10, 2018
 *      Author: guther
 */

#include "Solver.hpp"

Solver::Solver(double eta):gamma(eta) {

}

Solver::~Solver() {
}

Eigen::VectorXd Solver::update(Eigen::VectorXd const &w, Eigen::VectorXd const &force) const{
	return w-gamma*force;
}

Eigen::VectorXd Solver::update(Eigen::VectorXd const &w, Eigen::VectorXd const &force,
					   Eigen::VectorXd const &ci, Eigen::MatrixXd const &dcdw) const{
	Eigen::MatrixXd s, okokp;
	Eigen::VectorXd ok;
	double normalizer = ci.norm();
	ok = w.conjugate()*dcdw/normalizer;
	okokp = dcdw.conjugate()*dcdw/normalizer;
	s = okokp - ok.adjoint()*ok;
	return w-gamma*s.inverse()*force;
}
