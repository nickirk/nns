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
        std::cout << "force=" << std::endl;
        std::cout <<   force  << std::endl;
	w-=gamma*force;
}

void Solver::update(Eigen::VectorXd &w, Eigen::VectorXd const &force,
					   Eigen::VectorXcd const &ci, Eigen::MatrixXcd const &dcdw, int const &iteration) const{
	Eigen::MatrixXd s, okokp;
	Eigen::VectorXcd ok;
	double normalizer = ci.norm();
        //std::cout << "normalizer=" << normalizer << std::endl;
	ok = dcdw*ci.conjugate()/normalizer;
	okokp = (dcdw*dcdw.adjoint()/normalizer).real();
        //std::cout << "okokp=" << std::endl;       
        //std::cout << okokp << std::endl;       
        //std::cout << "ok*ok=" << std::endl;       
        //std::cout << (ok*ok.adjoint()).real() << std::endl;       
	s = okokp - (ok*ok.adjoint()).real();
        //Eigen::MatrixXd diag(1.1* Eigen::VectorXd::Ones(w.size()).asDiagonal());
        double lamda = std::max(100*std::pow(0.999,iteration), 1e-2);
        //s+=s.diagonal().asDiagonal()*lamda;
        s+=Eigen::VectorXd::Ones(w.size()).asDiagonal()*lamda;
        //std::cout << "s=" << std::endl;       
        //std::cout << s << std::endl;       
        //std::cout << "s^-1=" << std::endl;       
        //std::cout << s.inverse() << std::endl;       
        std::cout << " change of para=" << std::endl;
        std::cout <<   gamma*s.inverse()*force  << std::endl;
	w-=gamma*s.inverse()*force;
}
