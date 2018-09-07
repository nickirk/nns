/*
 * RBM.cxx
 *
 *  Created on: Apr 29, 2018
 *      Author: ghanem and Liao
 */

#include "RBM.hpp"
#include <new>

namespace networkVMC{

template <typename F, typename coeffType>
RBM<F, coeffType>::RBM(int sizeInput_, int sizeHidden_):
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

template <typename F, typename coeffType>
RBM<F, coeffType>::RBM(RBM<F, coeffType> const& source):
	a(nullptr, source.sizeInput),b(nullptr, source.sizeHidden),
	w(nullptr, source.sizeHidden, source.sizeInput){  // the maps have to be initialized, but leave them empty for now
	// delegate this to the copy-assignment
	*this = source;
}

//---------------------------------------------------------------------------//

template <typename F, typename coeffType>
RBM<F, coeffType>& RBM<F, coeffType>::operator=(RBM<F, coeffType> const &source){
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
		new (&a) Eigen::Map<T>(pars_vec.data()+a_offset, sizeInput);
		new (&b) Eigen::Map<T>(pars_vec.data()+b_offset, sizeHidden);
		new (&w) Eigen::Map<Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> >
			(pars_vec.data()+w_offset, sizeHidden, sizeInput);

	}
	return *this;
}

//---------------------------------------------------------------------------//

template<typename F, typename coeffType>
RBM<F, coeffType>::~RBM(){}

//---------------------------------------------------------------------------//


template <typename F, typename coeffType>
typename RBM<F, coeffType>::T const& RBM<F, coeffType>::pars() const
{
	return pars_vec;
}

//---------------------------------------------------------------------------//


template <typename F, typename coeffType>
coeffType RBM<F, coeffType>::getCoeff(detType const &det) const
{
	T s = T::Ones(sizeInput);

	for (int j=0; j<sizeInput; j++)
	{
		s(j) = det[j]?1.0:-1.0;
	}

	T coshVec = (b+w*s).array().cosh();
	F psi = std::exp(a.dot(s));
	for(int i=0; i<sizeHidden;i++)
	{
		psi = psi * coshVec(i);
	}
	//if (verbatimCast(det)==15) psi*=1;
	//else psi*=0.2;

	return psi;
}

//---------------------------------------------------------------------------//


template <typename F, typename coeffType>
typename RBM<F, coeffType>::T RBM<F, coeffType>::getDeriv(detType const &det) const
{
	T dCdWk= T::Zero(numPars);
	// do the mapping inside for loop, private
	Eigen::Map<T> da(dCdWk.data()+a_offset, sizeInput);
	Eigen::Map<T> db(dCdWk.data()+b_offset, sizeHidden);
	Eigen::Map<Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic>> dw(dCdWk.data()+w_offset, sizeHidden, sizeInput);
	T s = T::Ones(sizeInput);

	for (int j=0; j<sizeInput; j++)
	{
		s(j) = det[j]?1.0:-1.0;
	}

	T coshVec = (b+w*s).array().cosh();

	F psi = std::exp(a.dot(s));
	for(int i=0; i<sizeHidden;i++)
	{
		psi = psi * coshVec(i);
	}

	da = psi*s;

	T tanhVec = (b+w*s).array().tanh();
	db = psi*tanhVec;

	dw = psi*tanhVec*s.transpose();
	assert(dw.rows()==sizeHidden);
	assert(dw.cols()==sizeInput);
	return dCdWk;
}

//---------------------------------------------------------------------------//


//instantiate
template class RBM<double, double>;

template class RBM<std::complex<double>, std::complex<double>>;

}
