/*
 * RBM.cpp
 *
 *  Created on: Apr 29, 2018
 *      Author: ghanem
 */

#include "RBM.hpp"
#include <iostream>

namespace networkVMC
{
    RBM::RBM(int sizeInput, int sizeHidden):
    sizeHidden(sizeHidden),
    sizeInput(sizeInput),
    numPars(sizeInput+sizeHidden+sizeInput*sizeHidden),
    a_offset(0),
    b_offset(sizeInput),
    w_offset(sizeInput+sizeHidden),
    pars_vec(numPars),
    a(pars_vec.data()+a_offset, sizeInput),
    b(pars_vec.data()+b_offset, sizeHidden),
    w(pars_vec.data()+w_offset, sizeHidden, sizeInput)
    {
        pars_vec.setRandom();
        a.normalize();
        b.normalize();
        w.normalize();

    }

    RBM::~RBM()
    {
        
    }

    Eigen::VectorXcd const& RBM::pars() const
    { 
        return pars_vec;
    }

    coeffType RBM::getCoeff(detType const &det) const
    {
        Eigen::VectorXcd s = Eigen::VectorXcd::Ones(sizeInput);

        for (int j=0; j<sizeInput; j++)
        {
            s(j) = det[j]?1.0:-1.0;
        }
        
        Eigen::VectorXcd coshVec = (b+w*s).array().cosh();
        coeffType psi = std::exp(a.dot(s));
        for(int i=0; i<sizeHidden;i++)
        {
            psi = psi * coshVec(i);
        }

        return psi;
    }

    void RBM::getDeriv(detType const &det, 
                       Eigen::Map<Eigen::VectorXcd> & da, 
                       Eigen::Map<Eigen::VectorXcd> & db, 
                       Eigen::Map<Eigen::MatrixXcd> & dw) const
    {
        Eigen::VectorXcd s = Eigen::VectorXcd::Ones(sizeInput);

        for (int j=0; j<sizeInput; j++)
        {
            s(j) = det[j]?1.0:-1.0;
        }

        Eigen::VectorXcd coshVec = (b+w*s).array().cosh();

        coeffType psi = std::exp(a.dot(s));
        for(int i=0; i<sizeHidden;i++)
        {
            psi = psi * coshVec(i);
        }

        da = psi*s;
        
        Eigen::VectorXcd tanhVec = (b+w*s).array().tanh();
        db = psi*tanhVec;

        dw = psi*tanhVec*s.transpose();
        assert(dw.rows()==sizeHidden);
        assert(dw.cols()==sizeInput);
    }

    Eigen::MatrixXcd RBM::calcdCdwSR(State const &outputState)
    {
        Eigen::MatrixXcd result(numPars, outputState.size());
        
        #pragma omp parallel for
        for(int i=0;i<outputState.size();i++)
        {
            Eigen::Map<Eigen::VectorXcd> da(result.data()+i*numPars+a_offset, sizeInput);
            Eigen::Map<Eigen::VectorXcd> db(result.data()+i*numPars+b_offset, sizeHidden);
            Eigen::Map<Eigen::MatrixXcd> dw(result.data()+i*numPars+w_offset, sizeHidden, sizeInput);;
            getDeriv(outputState.det(i), da, db, dw);
        }
        return result;
    }

    VecCType RBM::calcNablaPars(State const &input, nablaType const &outerDerivative)
    {
        Eigen::MatrixXcd dCdW=calcdCdwSR(input);
        Eigen::VectorXcd result= Eigen::VectorXcd::Zero(numPars);

        //There is no reduction for vectors of Eigen, so we use critical section instead.
        //Note: OpenMP 4.0 supports custom reduction. But it is not supported in my current gcc.
        #pragma omp parallel
        {
            Eigen::VectorXcd private_result= Eigen::VectorXcd::Zero(numPars);
            #pragma omp for
            for(int i=0;i<input.size();i++)
            {
                private_result += (outerDerivative[i]*dCdW.col(i)).conjugate();
            }

            #pragma omp critical
            result+=private_result;
        }

        return result;
    }

}