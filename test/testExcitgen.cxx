#include "../src/Hamiltonian/ExcitationGenerators/AllExcitationGenerators.hpp"
#include "defaultSystem.hpp"
#include "../src/Hamiltonian/AbInitioHamiltonian.hpp"
#include "../src/utilities/SpinConfig.hpp"
#include "../src/HilbertSpace/Basis.hpp"
#include <assert.h>
#include <iostream>
#include <chrono>
#include <cmath>

using namespace networkVMC;

double eps{1e-10};
void testExcitgen(ExcitationGenerator &eg, detType const &HF);

void testProbUpdater(detType const &HF){
	ExcitationGenerator::ProbUpdater pBiasGen(HF);
	assert(std::fabs(pBiasGen.pParallel() - 0.333333) < 0.001);
}

void testExcitgen(ExcitationGenerator &eg, detType const &HF){
	// time the excitgen
	std::chrono::steady_clock::time_point tS = std::chrono::steady_clock::now();
	double prob{0.0};
	for(int i=0; i<10; ++i){
		auto excit = eg.generateExcitation(HF,prob);
		int exLvl = getExcitLvl(excit,HF);
		assert((exLvl == 1 or exLvl == 0 or exLvl == 2));
		double probCheck = eg.getExcitationProb(HF,excit);
		std::cout << "Generated lvl " << exLvl << std::endl;
		std::cout << "probCheck = " << probCheck << " vs " << prob << std::endl;
		assert(std::fabs(probCheck - prob) < eps);
	}
	eg.updateBiases();
	// how long did it take?
	std::chrono::duration<double> dT = std::chrono::duration_cast<std::chrono::duration<double> >(
			std::chrono::steady_clock::now()-tS);
	std::cout << "Excitgen took " << dT.count() << " seconds\n";
}

void testUniformExcitgen(detType const &HF){
	std::cout << "Testing uniform excitgen\n";
	UniformExcitgen ufg(HF);
	testExcitgen(ufg,HF);
}

void testHubbardExcitgen(detType const &HF){
	std::cout << "Testing RSHubbard excitgen\n";
	RSHubbardExcitgen hbg{};
	testExcitgen(hbg,HF);
}

void testLatticeExcitgen(LatticeHamiltonian const &H, detType const &HF){
	std::cout << "Testing lattice excitgen\n";
	LatticeExcitgen lg(H);
	testExcitgen(lg,HF);
}

void testWeightedExcitgen(Hamiltonian const &H, detType const &HF){
	std::cout << "Testing weighted excitgen\n";
	WeightedExcitgen weg(H,HF);
	// Note that due to some rounding (ask Lauretta about this)
	// the weighted excitgen can also create doubles for the
	// real-space hubbard in the default config
	testExcitgen(weg,HF);
}

void testExcitmatType(){
	ExcitmatType excitmat;
	assert(excitmat(1,1) == -1);
	excitmat(1,0) += 1;
	assert(excitmat(1,0) == 0);
}

int main(){
	// get default system parameters
	int numSites{32};
	auto defaultH = generateDefaultHubbard(numSites);
	// run the tests
	testExcitmatType();

	detType HF(2*numSites);
	// AFM initial state
	for(int i = 0; i < numSites; ++i){
		if(i%2){
			HF[2*i] = true;
		}
		else{
			HF[2*i+1] = true;
		}
	}
	if(numSites == 4) testProbUpdater(HF);

	testUniformExcitgen(HF);
	testHubbardExcitgen(HF);
	testLatticeExcitgen(defaultH,HF);
	// now test the weighted excitgen, this requires a 2-body Hamiltonian
	std::string file_name = "FCIDUMP";
	auto modelHam = readAbInitioHamiltonian(file_name);
	SpinConfig sC(2,2,2*19);
	Basis abBasis(sC,modelHam);
	testWeightedExcitgen(modelHam,abBasis.getDetByIndex(0));
}
