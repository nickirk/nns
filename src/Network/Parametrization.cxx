/*
 * Parametrization.cxx
 *
 *  Created on: Aug 23, 2018
 *      Author: guther
 */


#include "Parametrization.hpp"
#include "calcNabla.hpp"
#include <fstream>
#include <iterator>
#include <string>
#include <utility>

namespace networkVMC{

template<typename F, typename coeffType>
void Parametrization<F,coeffType>::writeParsToFile(std::string file) const{
  T tmpPars = pars();
  std::vector<F> parsStd(tmpPars.data(), tmpPars.data()+tmpPars.size());
  std::ofstream fout(file);
  fout.precision(10);
  std::copy(parsStd.begin(), parsStd.end(),
  std::ostream_iterator<F>(fout, "\n"));
  fout.close();
};

//---------------------------------------------------------------------------//

template<typename F, typename coeffType>
void Parametrization<F,coeffType>::readParsFromFile(std::string file){
   std::ifstream input(file);
   if (!input){
      throw FileNotFound(file);
      exit(2);
   }

   T &innerPars = pars();
   std::vector<F> buffer{
     std::istream_iterator<F>(input),
     std::istream_iterator<F>() };
   assert (buffer.size() == pars().size());
   innerPars = Eigen::Map<T> (buffer.data(),buffer.size());
   //use the pars() with write access
 };

//---------------------------------------------------------------------------//

template<typename F, typename coeffType>
Parametrization<F, coeffType>::T Parametrization<F,coeffType>::calcNablaParsConnected(
  State<coeffType> const &inputState, Eigen::Matrix<F, Eigen::Dynamic ,1> const &dEdC){
	return calcNabla(inputState, dEdC, F(), [](detType const &det){return this->getDeriv(det);},
			PreFetched, getNumPars());
}

//---------------------------------------------------------------------------//

template<typename F, typename coeffType>
Parametrization<F, coeffType>::T Parametrization<F,coeffType>::calcNablaParsMarkovConnected(
		State<coeffType> const &inputState, Eigen::Matrix<F, Eigen::Dynamic ,1> const& dEdC, F const& energy){
	return calcNabla(inputState, dEdC, energy, [](detType const &det){return this->getMarkovDeriv(det);},
			Markov, getNumPars());
}

template class Parametrization<double, double>;
template class Parametrization<cType, cType>;

}

