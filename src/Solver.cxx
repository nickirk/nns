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
					   Eigen::VectorXcd const &ci, Eigen::MatrixXcd const &dcdw) const{
	Eigen::MatrixXd s, okokp;
	Eigen::VectorXcd ok;
	double normalizer = ci.norm();
	ok = dcdw*w.conjugate()/normalizer;
	okokp = dcdw.conjugate()*dcdw/normalizer;
	s = (okokp - ok*ok.adjoint()).real();
	return w-gamma*s.inverse()*force;
}
