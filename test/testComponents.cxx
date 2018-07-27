#include "testComponents.hpp"
#include "defaultSystem.hpp"
#include <iostream>

using namespace networkVMC;

// use ADAM to solve
void solveADAM(Parametrization<> &para, Sampler &msampler, Hamiltonian const &H){
	ADAM<> solver(0.01);
	EnergyEs eCF(H,-1);
	Trainer<> etr(para,msampler,solver,eCF,H);
	for(int i = 1; i < 1000; ++i){
		etr.train();
		std::cout << "Current energy: " << etr.getE()<<std::endl;
	}
}

// use the full sampler with the direct parametrization
void testDeterministicFullSampling(SpinConfig const &sC, Hamiltonian const &H){
	Basis basis(sC,H);
	DirectParametrization<> para(basis);
	auto HF = basis.getDetByIndex(0);
	FullSampler<> mySampler(H,basis,para);
	solveADAM(para,mySampler,H);
}

void testAdj(LatticeHamiltonian const &test){
	std::cout << test.size() << std::endl;
	for(int i = 0; i < test.size(); ++i){
		std::cout << "Adjacent to " << i << ": ";
		for(auto j:test.adjacents(i)){
			std::cout << j << " ";
		}
		std::cout<<std::endl;
	}
}
