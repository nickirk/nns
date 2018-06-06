/*
 * RBM.hpp
 *
 *  Created on: Apr 29, 2018
 *      Author: ghanem
 */
#ifndef RBM_DEFINED
#define RBM_DEFINED

#include <Eigen/Dense>
#include "Parametrization.hpp"
#include "../utilities/State.hpp"
#include "../CostFunctions/CostFunction.hpp"

namespace networkVMC
{
    //Resttricted Boltzmann Machine
    class RBM: public Parametrization<VecCType>
    {
        public:
            RBM(int sizeInput, int sizeHidden);

            coeffType getCoeff(detType const &det) const;
            void getDeriv(detType const &det, Eigen::Map<Eigen::VectorXcd> & da, Eigen::Map<Eigen::VectorXcd> & db, Eigen::Map<Eigen::MatrixXcd> & dw) const;
            virtual VecCType const& pars() const;
            virtual VecCType calcNablaPars(State const &input, nablaType const &outerDerivative);
            virtual Eigen::MatrixXcd calcdCdwSR(State const &outputState);
            ~RBM();

        private:
            int sizeInput;
            int sizeHidden;        
            int numPars;
            int a_offset;
            int b_offset;
            int w_offset;
            Eigen::VectorXcd pars_vec;
            Eigen::Map<Eigen::VectorXcd> a;
            Eigen::Map<Eigen::VectorXcd> b;
            Eigen::Map<Eigen::MatrixXcd> w;            
    };
}

#endif //RBM_DEFINED