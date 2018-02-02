/*
 * Solver.cxx
 *
 *  Created on: Jan 10, 2018
 *      Author: guther
 */

#include "Solver.hpp"
#include <iostream>

Solver::Solver(double eta):gamma(eta) {

}

Solver::~Solver() {
}

void Solver::update(Eigen::VectorXd &w, Eigen::VectorXd const &force) const{
	w-=gamma*force;
}

void Solver::update(Eigen::VectorXd &w, Eigen::VectorXd const &force,
					   Eigen::VectorXcd const &ci, Eigen::MatrixXcd const &dcdw, int const &iteration) const{
	Eigen::MatrixXd s, okokp;
	Eigen::VectorXcd ok;
	double normalizer = ci.norm();
	ok = dcdw*ci.conjugate()/normalizer;
	okokp = (dcdw*dcdw.adjoint()/normalizer).real();
	s = okokp - (ok*ok.adjoint()).real();
        double lamda = std::max(100*std::pow(0.999,iteration), 1e-2);
        s+=s.diagonal().asDiagonal()*lamda;
	w-=gamma*s.inverse()*force;
}
