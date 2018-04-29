/*
 * RBM.cpp
 *
 *  Created on: Apr 29, 2018
 *      Author: ghanem
 */

#include "RBM.hpp"

namespace networkVMC
{
    RBM::RBM(int sizeInput, int sizeHidden)
    {
        a = Eigen::VectorXd::Random(sizeInput);
        b = Eigen::VectorXd::Random(sizeHidden);
        w = Eigen::MatrixXd::Random(sizeHidden, sizeInput);
        this->sizeHidden = sizeHidden;
        this->sizeInput = sizeInput;
    }

    coeffType RBM::getCoeff(detType const &det) const
    {
        Eigen::VectorXd s = Eigen::VectorXd::Ones(sizeInput);

        for (int j=0; j<sizeInput; j++)
        {
            s(j) = det[j]?1.0:-1.0;
        }

        Eigen::VectorXd coshVec = (b+w*s).array().cosh();
        double psi = std::exp(a.dot(s));
        for(int i=0; i<sizeHidden;i++)
        {
            psi = psi* 2* coshVec(i);//2 * cosh(b(i)+w.row(i).dot(s));
        }

        return psi;
    }

    void RBM::getDeriv(detType const &det, Eigen::VectorXd & da, Eigen::VectorXd & db, Eigen::MatrixXd & dw) const
    {
        Eigen::VectorXd s = Eigen::VectorXd::Ones(sizeInput);

        for (int j=0; j<sizeInput; j++)
        {
            s(j) = det[j]?1.0:-1.0;
        }

        Eigen::VectorXd coshVec = (b+w*s).array().cosh();

        double psi = std::exp(a.dot(s));
        for(int i=0; i<sizeHidden;i++)
        {
            psi = psi* 2* coshVec(i);//2 * cosh(b(i)+w.row(i).dot(s));
        }

        da = psi*s;
        
        Eigen::VectorXd tanhVec = (b+w*s).array().tanh();
        db = psi*tanhVec;

        dw = psi*tanhVec*s.transpose();
    }

    RBM::~RBM()
    {
        
    }
}