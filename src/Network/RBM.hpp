/*
 * RBM.hpp
 *
 *  Created on: Apr 29, 2018
 *      Author: ghanem and Liao
 */
#ifndef RBM_DEFINED
#define RBM_DEFINED

#include <Eigen/Dense>
#include "../CostFunctions/CostFunction.hpp"
#include "Parametrization.hpp"
#include "../utilities/State.hpp"
#include "../utilities/TypeDefine.hpp"

namespace networkVMC
{
    //Resttricted Boltzmann Machine
    class RBM: public Parametrization<VecCType>
    {
        public:
            RBM(int sizeInput, int sizeHidden);

            coeffType getCoeff(detType const &det) const;
            virtual Eigen::VectorXcd getDeriv(detType const &det) const;
            virtual Eigen::VectorXcd getSRDeriv(detType const &det) const;
            virtual VecCType const& pars() const;
            virtual VecCType calcNablaPars(State const &input, nablaType const &outerDerivative);
            virtual Eigen::MatrixXcd calcdCdwSR(State const &outputState);
            virtual VecCType calcNablaParsConnected(State const &inputState, nablaType const& dEdC);
            virtual VecCType calcNablaParsSRConnected(State const &inputState, nablaType const& dEdC, double const& energy);
            ~RBM();

        private:
            int sizeHidden;
            int sizeInput;
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
