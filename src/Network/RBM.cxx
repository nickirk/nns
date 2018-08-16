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

    Eigen::VectorXcd RBM::getDeriv(detType const &det) const
    {
        Eigen::VectorXcd dCdWk= Eigen::VectorXcd::Zero(numPars);
        // do the mapping inside for loop, private
        Eigen::Map<Eigen::VectorXcd> da(dCdWk.data()+a_offset, sizeInput);
        Eigen::Map<Eigen::VectorXcd> db(dCdWk.data()+b_offset, sizeHidden);
        Eigen::Map<Eigen::MatrixXcd> dw(dCdWk.data()+w_offset, sizeHidden, sizeInput);
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
        return dCdWk;
    }

    Eigen::VectorXcd RBM::getMarkovDeriv(detType const &det) const
    {
        Eigen::VectorXcd dCdWk= Eigen::VectorXcd::Zero(numPars);
        // do the mapping inside for loop, private
        Eigen::Map<Eigen::VectorXcd> da(dCdWk.data()+a_offset, sizeInput);
        Eigen::Map<Eigen::VectorXcd> db(dCdWk.data()+b_offset, sizeHidden);
        Eigen::Map<Eigen::MatrixXcd> dw(dCdWk.data()+w_offset, sizeHidden, sizeInput);

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

        da = s;
        
        Eigen::VectorXcd tanhVec = (b+w*s).array().tanh();
        db = tanhVec;

        dw = tanhVec*s.transpose();
        assert(dw.rows()==sizeHidden);
        assert(dw.cols()==sizeInput);
        return dCdWk;
    }

    Eigen::MatrixXcd RBM::calcdCdwSR(State const &outputState)
    {
        Eigen::MatrixXcd result(numPars, outputState.size());
        
        #pragma omp parallel for
        for(size_t i=0;i<outputState.size();i++)
        {
            Eigen::VectorXcd dCtdW= Eigen::VectorXcd::Zero(numPars);
            // do the mapping inside for loop, private
            //update vector dCidWk
            dCtdW = getDeriv(outputState.det(i));
            result.col(i) << (dCtdW);
        }
        return result;
    }

    // deprecated due to the change in the structure of outerDerivative
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

    // used in conjunction with Markov chain Metropolis sampling
    VecCType RBM::calcNablaParsMarkovConnected(
	   State const &inputState,
	   nablaType const &dEdC,
     double const &energy
     ){
        int numDets = inputState.size();
        int spaceSize = inputState.totalSize();
        Eigen::VectorXcd dEdW= Eigen::VectorXcd::Zero(numPars);
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
            //update vector dCidWk
            dCtdW = getMarkovDeriv(inputState.det(i));
            // multiplication should be done by matrix vector product
            // fill up the dCdW matrix
            dCdW.col(i) << (dCtdW.conjugate()*inputState.coeff(i));
            dEdW -= energy * dCtdW.conjugate()/numDets; 
            //dedc[i] = 1;
            std::vector<detType> coupledDets = inputState.coupledDets(i);
            std::vector<coeffType > coupledCoeffs = inputState.coupledCoeffs(i);
            size_t coupledSize = inputState.coupledDets(i).size();
            size_t pos = inputState.locate(i);
            for (size_t j(0); j < coupledSize; ++j){
              // fill up the dCdW matrix with coupled dets contribution
              dCdW.col(numDets+pos+j) << (dCtdW.conjugate()*coupledCoeffs[j]);
              //dEdWTmp +=  dCtdW * dEdC[pos];
            }
          }
        // map std::vector of dEdC to Eigen Vector
        Eigen::VectorXcd dEdCEigen=Eigen::Map<Eigen::VectorXcd>(dedc.data(),spaceSize);
        // make it parallel. TODO
        dEdW += (dCdW * dEdCEigen);//.conjugate();
        return dEdW;
    }
    // The following function is used when use ListGen or fullSampler
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
            //update vector dCidWk
            dCtdW = getDeriv(inputState.det(i));
            //update vector dCidWk
            // multiplication should be done by matrix vector product
            // fill up the dCdW matrix
            dCdW.col(i) << (dCtdW.conjugate());
            //dedc[i] = 1;
            std::vector<detType> coupledDets = inputState.coupledDets(i);
            std::vector<coeffType > coupledCoeffs = inputState.coupledCoeffs(i);
            size_t coupledSize = inputState.coupledDets(i).size();
            size_t pos = inputState.locate(i);
            for (size_t j(0); j < coupledSize; ++j){
              //update dCjdW
              //dCtdW = getDeriv(coupledDets[j]);
              // fill up the dCdW matrix with coupled dets contribution
              dCdW.col(numDets+pos+j) << (dCtdW.conjugate());
              //dEdWTmp +=  dCtdW * dEdC[pos];
            }
          }
        // map std::vector of dEdC to Eigen Vector
        Eigen::VectorXcd dEdCEigen=Eigen::Map<Eigen::VectorXcd>(dedc.data(),spaceSize);
        // make it parallel. TODO
        dEdW = (dCdW * dEdCEigen);//.conjugate();
        return dEdW;
      }

}
