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
#include "../math/MathFunctions.hpp"

namespace networkVMC
{
    //Resttricted Boltzmann Machine
    class RBM: public ClonableParametrization<VecCType,RBM>
    {
        public:
            RBM(int sizeInput_, int sizeHidden_);

            coeffType getCoeff(detType const &det) const;
<<<<<<< HEAD
            virtual Eigen::VectorXcd getDeriv(detType const &det) const;
            //virtual Eigen::VectorXcd getMarkovDeriv(detType const &det) const;
            virtual VecCType const& pars() const;
            virtual VecCType calcNablaPars(State const &input, nablaType const &outerDerivative);
            virtual Eigen::MatrixXcd calcdCdwSR(State const &outputState);
            virtual VecCType calcNablaParsConnected(State const &inputState, nablaType const& dEdC);
            //virtual VecCType calcNablaParsMarkovConnected(State const &inputState, nablaType const& dEdC, coeffType const& energy);
            ~RBM();
=======
            Eigen::VectorXcd getDeriv(detType const &det) const;
            Eigen::VectorXcd getMarkovDeriv(detType const &det) const;
            VecCType const& pars() const;
            VecCType calcNablaPars(State const &input, nablaType const &outerDerivative);
            Eigen::MatrixXcd calcdCdwSR(State const &outputState);
            VecCType calcNablaParsConnected(State const &inputState, nablaType const& dEdC);
            VecCType calcNablaParsMarkovConnected(State const &inputState, nablaType const& dEdC, double const& energy);
>>>>>>> 7e5214b8a261068e7eee1c2169c614655df65973

        private:
            int sizeHidden;
            int numPars;
            int sizeInput;
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
