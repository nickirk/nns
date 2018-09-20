/*
 * RBM.cxx
 *
 *  Created on: Apr 29, 2018
 *      Author: ghanem and Liao
 */

#include "RBM.hpp"
#include <new>

namespace networkVMC{

RBM::RBM(int sizeInput_, int sizeHidden_):
sizeHidden(sizeHidden_),
sizeInput(sizeInput_),
numPars(sizeInput_+sizeHidden_+sizeInput_*sizeHidden_),
a_offset(0),
b_offset(sizeInput_),
w_offset(sizeInput_+sizeHidden_),
pars_vec(numPars),
a(pars_vec.data()+a_offset, sizeInput_),
b(pars_vec.data()+b_offset, sizeHidden_),
w(pars_vec.data()+w_offset, sizeHidden_, sizeInput_){
	//pars_vec.real() = T::Ones(numPars)*0.01;
	//pars_vec.imag() = T::Ones(numPars)*0.01;
	pars_vec.setRandom(numPars);
	pars_vec.normalize();
	// pars_vec.real() = pars_vec.real().unaryExpr(&NormalDistribution);
	// pars_vec.imag() = pars_vec.imag().unaryExpr(&NormalDistribution);
	//a.setZero();
	//b.normalize();
	//w.normalize();
}

//---------------------------------------------------------------------------//

RBM::RBM(RBM const& source):
	a(nullptr, source.sizeInput),b(nullptr, source.sizeHidden),
	w(nullptr, source.sizeHidden, source.sizeInput){  // the maps have to be initialized, but leave them empty for now
	// delegate this to the copy-assignment
	*this = source;
}

//---------------------------------------------------------------------------//

RBM& RBM::operator=(RBM const &source){
	if(&source != this){
		// copy the data
		sizeHidden = source.sizeHidden;
		numPars = source.numPars;
		sizeInput = source.sizeInput;
		pars_vec = source.pars_vec;
		a_offset = source.a_offset;
		b_offset = source.b_offset;
		w_offset = source.w_offset;

		// then, set the maps to the copied datasets (we dont need to do this when moving data,
		// which is why we default the move ctor)

		// this is done using the placement-new syntax
		new (&a) Eigen::Map<paraVector>(pars_vec.data()+a_offset, sizeInput);
		new (&b) Eigen::Map<paraVector>(pars_vec.data()+b_offset, sizeHidden);
		new (&w) Eigen::Map<Eigen::Matrix<paraType, Eigen::Dynamic, Eigen::Dynamic> >
			(pars_vec.data()+w_offset, sizeHidden, sizeInput);

	}
	return *this;
}

//---------------------------------------------------------------------------//

RBM::~RBM(){}

//---------------------------------------------------------------------------//

paraVector const& RBM::pars() const
{
	return pars_vec;
}

//---------------------------------------------------------------------------//

coeffType RBM::getCoeff(detType const &det) const
{
	paraVector s = paraVector::Ones(sizeInput);

	for (int j=0; j<sizeInput; j++)
	{
		s(j) = det[j]?1.0:-1.0;
	}

	paraVector coshVec = (b+w*s).array().cosh();
	paraType psi = std::exp(a.dot(s));
	for(int i=0; i<sizeHidden;i++)
	{
		psi = psi * coshVec(i);
	}
	//if (verbatimCast(det)==15) psi*=1;
	//else psi*=0.2;

	return psi;
}

//---------------------------------------------------------------------------//


paraVector RBM::getDeriv(detType const &det) const
{
	paraVector dCdWk= paraVector::Zero(numPars);
	// do the mapping inside for loop, private
	Eigen::Map<paraVector> da(dCdWk.data()+a_offset, sizeInput);
	Eigen::Map<paraVector> db(dCdWk.data()+b_offset, sizeHidden);
	Eigen::Map<Eigen::Matrix<paraType, Eigen::Dynamic, Eigen::Dynamic>> dw(dCdWk.data()+w_offset, sizeHidden, sizeInput);
	paraVector s = paraVector::Ones(sizeInput);

	for (int j=0; j<sizeInput; j++)
	{
		s(j) = det[j]?1.0:-1.0;
	}

	paraVector coshVec = (b+w*s).array().cosh();

	paraType psi = std::exp(a.dot(s));
	for(int i=0; i<sizeHidden;i++)
	{
		psi = psi * coshVec(i);
	}

	da = psi*s;

	paraVector tanhVec = (b+w*s).array().tanh();
	db = psi*tanhVec;

	dw = psi*tanhVec*s.transpose();
	assert(dw.rows()==sizeHidden);
	assert(dw.cols()==sizeInput);
	return dCdWk;
}

//---------------------------------------------------------------------------//

}
