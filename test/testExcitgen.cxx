#include "../src/Hamiltonian/ExcitationGenerators/AllExcitationGenerators.hpp"
#include "defaultSystem.hpp"
#include <assert.h>
#include <iostream>
#include <cmath>

using namespace networkVMC;

double eps{1e-10};
void testExcitgen(ExcitationGenerator &eg, detType const &HF);

void testProbUpdater(detType const &HF){
	ExcitationGenerator::ProbUpdater pBiasGen(HF);
	assert(std::fabs(pBiasGen.pParallel() - 0.333333) < 0.001);
}

void testExcitgen(ExcitationGenerator &eg, detType const &HF){
	double prob{0.0};
	for(int i=0; i<10; ++i){
		auto excit = eg.generateExcitation(HF,prob);
		int exLvl = getExcitLvl(excit,HF);
		assert((exLvl == 1 or exLvl == 0 or exLvl == 2));
		double probCheck = eg.getExcitationProb(HF,excit);
		std::cout << "Generated lvl " << exLvl << std::endl;
		assert(std::fabs(probCheck - prob) < eps);
	}
	eg.updateBiases();
}

void testUniformExcitgen(detType const &HF){
	UniformExcitgen ufg(HF);
	testExcitgen(ufg,HF);
}

void testHubbardExcitgen(detType const &HF){
	RSHubbardExcitgen hbg{};
	// Note that due to some rounding (ask Lauretta about this)
	// the weighted excitgen can also create doubles for the
	// real-space hubbard in the default config
	testExcitgen(hbg,HF);
}

void testWeightedExcitgen(Hamiltonian const &H, detType const &HF){
	WeightedExcitgen weg(H,HF);
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
	int numSites{4};
	auto defaultH = generateDefaultHubbard(numSites);
	auto defaultBasis = generateDefaultBasis(numSites);

	// run the tests
	testExcitmatType();

	auto HF = defaultBasis.getDetByIndex(0);
	testProbUpdater(HF);

	testUniformExcitgen(HF);
	testHubbardExcitgen(HF);
	testWeightedExcitgen(defaultH,HF);
}
