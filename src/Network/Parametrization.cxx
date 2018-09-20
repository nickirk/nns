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

void Parametrization::writeParsToFile(std::string file) const{
  paraVector tmpPars = pars();
  std::vector<paraType> parsStd(tmpPars.data(), tmpPars.data()+tmpPars.size());
  std::ofstream fout(file);
  fout.precision(10);
  std::copy(parsStd.begin(), parsStd.end(),
  std::ostream_iterator<paraType>(fout, "\n"));
  fout.close();
};

//---------------------------------------------------------------------------//

void Parametrization::readParsFromFile(std::string file){
   std::ifstream input(file);
   if (!input){
      throw errors::FileNotFound(file);
   }

   paraVector &innerPars = pars();
   std::vector<paraType> buffer{
     std::istream_iterator<paraType>(input),
     std::istream_iterator<paraType>() };
   assert (buffer.size() == pars().size());
   innerPars = Eigen::Map<paraVector> (buffer.data(),buffer.size());
   //use the pars() with write access
 };

//---------------------------------------------------------------------------//

paraVector Parametrization::calcNablaParsConnected(
  State const &inputState, paraVector const &dEdC){
	return calcNabla(inputState, dEdC, paraType(), [this](detType const &det){return this->getDeriv(det);},
			PreFetched, getNumPars());
}

//---------------------------------------------------------------------------//

paraVector Parametrization::calcNablaParsMarkovConnected(
		State const &inputState, paraVector const& dEdC, paraType const& energy){
	return calcNabla(inputState, dEdC, energy, [this](detType const &det){return this->getMarkovDeriv(det);},
			Markov, getNumPars());
}

}

