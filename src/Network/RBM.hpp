/*
 * RBM.hpp
 *
 *  Created on: Apr 29, 2018
 *      Author: ghanem
 */
#ifndef RBM_DEFINED
#define RBM_DEFINED

#include <Eigen/Dense>
#include "../utilities/TypeDefine.hpp"

namespace networkVMC
{
    //Resttricted Boltzmann Machine
    class RBM
    {
        public:
            RBM(int sizeInput, int sizeHidden);

            coeffType getCoeff(detType const &det) const;
            void getDeriv(detType const &det, Eigen::VectorXd & da, Eigen::VectorXd & db, Eigen::MatrixXd & dw) const;

            ~RBM();

        private:
            Eigen::VectorXd a;
            Eigen::VectorXd b;
            Eigen::MatrixXd w;
            int sizeInput;
            int sizeHidden;
    };
}

#endif //RBM_DEFINED