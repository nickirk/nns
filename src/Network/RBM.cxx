/*
 * RBM.cxx
 *
 *  Created on: Apr 29, 2018
 *      Author: ghanem and Liao
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
        for(size_t i=0;i<outputState.size();i++)
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
            for(size_t i=0;i<input.size();i++)
            {
                private_result += (outerDerivative[i]*dCdW.col(i)).conjugate();
            }

            #pragma omp critical
            result+=private_result;
        }

        return result;
    }

    VecCType RBM::calcNablaParsConnected(
	   State const &inputState,
	   nablaType const &dEdC
     ){
        int numDets = inputState.size();
        // spaceSize = size of sampled dets and their coupled ones
        int spaceSize = inputState.totalSize();
        Eigen::VectorXcd dEdW= Eigen::VectorXcd::Zero(numPars);
        // Eigen::VectorXcd dEdWTmp= Eigen::VectorXcd::Zero(numPars);
        Eigen::MatrixXcd dCdW = Eigen::MatrixXcd::Zero(numPars, spaceSize);
        std::vector<std::complex<double>> dedc=dEdC;
          #pragma omp for
          // fill up the matrix of dCdW, like in EnergyEsMarkov.cxx
          // reserve space and in the end use matrix*vector instead of
          // a summation
          for (int i=0; i < numDets; ++i){
            //need private dCtdW
            Eigen::VectorXcd dCtdW= Eigen::VectorXcd::Zero(numPars);
            // do the mapping inside for loop, private
            Eigen::Map<Eigen::VectorXcd> da(dCtdW.data()+a_offset, sizeInput);
            Eigen::Map<Eigen::VectorXcd> db(dCtdW.data()+b_offset, sizeHidden);
            Eigen::Map<Eigen::MatrixXcd> dw(dCtdW.data()+w_offset, sizeHidden, sizeInput);;
            //update vector dCidWk
            getDeriv(inputState.det(i), da, db, dw);
            // multiplication should be done by matrix vector product
            // fill up the dCdW matrix
            dCdW.col(i) << (2*dCtdW*dedc[i]).real();
            dedc[i] = 1;
            std::vector<detType> coupledDets = inputState.coupledDets(i);
            std::vector<coeffType > coupledCoeffs = inputState.coupledCoeffs(i);
            size_t coupledSize = inputState.coupledDets(i).size();
            size_t pos = inputState.locate(i);
            for (size_t j(0); j < coupledSize; ++j){
              //update dCjdW
              getDeriv(coupledDets[j], da, db, dw);
              // fill up the dCdW matrix with coupled dets contribution
              dCdW.col(numDets+pos+j) << (dCtdW);
              //dEdWTmp +=  dCtdW * dEdC[pos];
            }
          }
        // map std::vector of dEdC to Eigen Vector
        Eigen::VectorXcd dEdCEigen=Eigen::Map<Eigen::VectorXcd>(dedc.data(),spaceSize);
        // make it parallel. TODO
        dEdW = (dCdW * dEdCEigen).conjugate();
        return dEdW;
      }

}
