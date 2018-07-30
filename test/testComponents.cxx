#include "testComponents.hpp"
#include "defaultSystem.hpp"

using namespace networkVMC;

// use the full sampler with the direct parametrization
void testDeterministicFullSampling(SpinConfig const &sC, Hamiltonian const &H){
	Basis basis(sC,H);
	DirectParametrization<> para(basis);
	auto HF = basis.getDetByIndex(0);
	FullSampler<> mySampler(H,basis,para);
	solveADAM(para,mySampler,H);
}

// RBM + Metropolis sampler, using Stochastic Reconfiguration
void testRBMMetropolis(networkVMC::SpinConfig const &sC, networkVMC::Hamiltonian const &H){
	Basis basis(sC,H);
	int const numHidden = 20;
	// RBM parametrization
	RBM network(sC.numSpinOrbs(),numHidden);
	// Metropolis Sampler (needs complex coeffs)
	MetropolisSampler<VecCType> mySampler(H,basis.getDetByIndex(0),network);
	solveADAM(network, mySampler, H);
}

void testAdj(LatticeHamiltonian const &test){
	std::cout << "Number of sites: " << test.size() << std::endl;
	for(int i = 0; i < test.size(); ++i){
		std::cout << "Adjacent to " << i << ": ";
		for(auto j:test.adjacents(i)){
			std::cout << j << " ";
		}
		std::cout<<std::endl;
	}
}

void testBasis(Basis const &basis){
	std::cout << "size= " << basis.size() << std::endl;
	for (int i=0; i < basis.size(); ++i){
	  detType det;
	  det = basis.getDetByIndex(i);
	  int index(basis.getIndexByDet(det));
	  std::vector<int> pos(getOccupiedPositions(det));
	  std::cout << "index= " << index << std::endl;
	  std::cout << "i = " << i  << std::endl;
	  for (int j=0; j < pos.size(); ++j){
		std::cout << pos[j] << "," ;
	  }
	  std::cout << std::endl;
	}
}
