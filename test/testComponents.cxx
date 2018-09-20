#include "testComponents.hpp"
#include "defaultSystem.hpp"

using namespace networkVMC;

// use ADAM to solve
void solveEs(networkVMC::Parametrization &para, networkVMC::Sampler &msampler,
		networkVMC::Hamiltonian const &H, networkVMC::Solver &solver){
	networkVMC::EnergyEs eCF(H,-1);
	networkVMC::Trainer etr(para,msampler,solver,eCF,H);
	for(int i = 1; i < 4000; ++i){
		etr.train();
		std::cout << "On iteration " << i << std::endl;
		std::cout << "Current energy: " << etr.getE()<<std::endl;
	}
}

//---------------------------------------------------------------------------//

void solveADAM(networkVMC::Parametrization &para, networkVMC::Sampler &msampler,
		networkVMC::Hamiltonian const &H){
	networkVMC::ADAM solver(0.01);
    solveEs(para,msampler,H,solver);
}

//---------------------------------------------------------------------------//

void solveSRec(networkVMC::Parametrization &para, networkVMC::Sampler &msampler, networkVMC::Hamiltonian const &H){
	networkVMC::StochasticReconfiguration solver(para,0.01);
	solveEs(para,msampler,H,solver);
}

//---------------------------------------------------------------------------//

// use the full sampler with the direct parametrization
void testDeterministicFullSampling(SpinConfig const &sC, Hamiltonian const &H){
	Basis basis(sC,H);
	DirectParametrization para(basis);
	auto HF = basis.getDetByIndex(0);
	FullSampler mySampler(H,basis,para);
	solveADAM(para,mySampler,H);
}

//---------------------------------------------------------------------------//

// RBM + Metropolis sampler, using Stochastic Reconfiguration
void testRBMMetropolis(networkVMC::SpinConfig const &sC, networkVMC::Hamiltonian const &H){
	Basis basis(sC,H);
	int const numHidden = 20;
	// RBM parametrization
	RBM network(sC.numSpinOrbs(),numHidden);
	// Metropolis Sampler (needs complex coeffs)
	MetropolisSampler mySampler(H,basis.getDetByIndex(0), basis, network);
	solveADAM(network, mySampler, H);
}

//---------------------------------------------------------------------------//

void testTrialMetropolis(SpinConfig const &sC, Hamiltonian const &H){
	int const numHidden =  20;
	// Base parametrization
	RBM basePara(sC.numSpinOrbs(), numHidden);
	// Trial Wf
	Basis basis(sC,H);
	DirectParametrization twf(basis);

	TrialWfPara varWF(basePara, twf);

	MetropolisSampler mySampler(H,basis.getDetByIndex(0), basis, varWF);
	solveADAM(varWF, mySampler, H);
}

//---------------------------------------------------------------------------//

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

//---------------------------------------------------------------------------//

void testBasis(Basis const &basis){
	std::cout << "size= " << basis.size() << std::endl;
	for (size_t i=0; i < basis.size(); ++i){
	  detType det;
	  det = basis.getDetByIndex(i);
	  int index(basis.getIndexByDet(det));
	  std::vector<int> pos(getOccupiedPositions(det));
	  std::cout << "index= " << index << std::endl;
	  std::cout << "i = " << i  << std::endl;
	  for (size_t j=0; j < pos.size(); ++j){
		std::cout << pos[j] << "," ;
	  }
	  std::cout << std::endl;
	}
}
