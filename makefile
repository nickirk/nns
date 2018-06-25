CXX=g++
LINKFLAGS=$(LDFLAGS) 
LDFLAGS=-std=c++11 -ggdb -Wall -fopenmp
EIGEN_PATH=lib/eigen/
EIGEN_FLAGS=-I$(EIGEN_PATH)
SOURCEFILES=Hamiltonian/AbInitioHamiltonian.cxx Hamiltonian/ExcitationGenerators/defaultExcitgensMap.cxx Hamiltonian/FermiHubbardHamiltonian.cxx Hamiltonian/ExcitationGenerators/WeightedExcitgen.cxx Hamiltonian/ExcitationGenerators/ExcitationGenerator.cxx Hamiltonian/ExcitationGenerators/RSHubbardExcitgen.cxx Hamiltonian/ExcitationGenerators/UniformExcitgen.cxx Hamiltonian/ExcitationGenerators/ProbUpdater.cxx Hamiltonian/ExcitationGenerators/WeightedSelector.cxx Hamiltonian/ExcitationGenerators/WeightedSelector.cxx Samplers/Sampler.cxx Network/Nnw.cxx Network/RBM.cxx Network/DirectParametrization.cxx HilbertSpace/Basis.cxx Hamiltonian/Hamiltonian.cxx HilbertSpace/Determinant.cxx CostFunctions/EnergyCF.cxx CostFunctions/NormCF.cxx CostFunctions/EnergyEs.cxx CostFunctions/EnergyEsPreFetched.cxx CostFunctions/EnergyCF.cxx Hamiltonian/FermionicHamiltonian.cxx Hamiltonian/BosonicHamiltonian.cxx Trainer.cxx Samplers/MetropolisSampler.cxx Samplers/ListGen.cxx CostFunctions/EnergyEsMarkov.cxx Network/Layers/Layer.cxx Network/Layers/DenseLayer.cxx Network/Layers/ConvLayer.cxx Network/Layers/InputLayer.cxx Network/LayerStructure.cxx math/MathFunctions.cxx Hamiltonian/SparseHMatrix.cxx Solvers/ADAM.cxx Solvers/AcceleratedGradientDescent.cxx Samplers/FullSampler.cxx
TESTFILES=defaultSystem.cxx testExcitgen.cxx testSubspaceCF.cxx testStateSort.cxx testNnw.cxx testRBM.cxx testDirect.cxx testSampler.cxx testBasis.cxx testEigen.cxx testAlg.cxx testCF.cxx testPreTrain.cxx testAbInitioHam.cxx testAlgAb.cxx testRandom.cxx testMetropolis.cxx
SRC=src
TST=test
BUILD=build
TSTBUILD=build/test
PROFILE=profile/build
DDIR=.depend
TSTDDIR=.depend/test
DIRECTORIES=$(DDIR) $(BUILD) $(TSTDDIR) $(TSTBUILD)
SOURCES=$(patsubst %,$(SRC)/%,$(SOURCEFILES)) $(patsubst %,$(TST)/%,$(TESTFILES))
OBJECTS=$(patsubst %.cxx,$(BUILD)/%.o,$(SOURCEFILES)) 
OPTIMOBJ=$(patsubst %.cxx,$(PROFILE)/%.o,$(SOURCEFILES))
DEPENDENCIES=$(patsubst %.cxx,$(DDIR)/%.d,$(SOURCEFILES)) $(patsubst %.cxx,$(TSTDDIR)/%.d,$(TESTFILES))
BTEST=basisTest
NNWTEST=nnwTest
RBMTEST=rbmTest
ALGTEST=algTest
MARKOVTEST=metropolisTest
STEST=samplerTest
HTEST=hamTest
ETEST=eigenTest
CFTEST=cfTest
PTTEST=ptTest
ABINHAMTEST=abinhamTest
ALGABTEST=algAbTest
DIRECTTST=directTest
SORTTST=sortTest
SCFTST=subspaceCFTest
BTESTOBJECT=$(TSTBUILD)/testBasis.o
NNWTESTOBJECT=$(TSTBUILD)/testNnw.o
RBMTESTOBJECT=$(TSTBUILD)/testRBM.o
MARKOVTESTOBJECT=$(TSTBUILD)/testMetropolis.o
ALGTESTOBJECT=$(TSTBUILD)/testAlg.o
STESTOBJECT=$(TSTBUILD)/testSampler.o
ETESTOBJECT=$(TSTBUILD)/testEigen.o
HTESTOBJECT=$(TSTBUILD)/testHam.o
CFTESTOBJECT=$(TSTBUILD)/testCF.o
PTTESTOBJECT=$(TSTBUILD)/testPreTrain.o
ABINHAMTESTOBJECT=$(TSTBUILD)/testAbInitioHam.o
ALGABTESTOBJECT=$(TSTBUILD)/testAlgAb.o
RNGTESTOBJECT=$(TSTBUILD)/testRandom.o
DIRECTTSTOBJECT=$(TSTBUILD)/testDirect.o
SORTTSTOBJECT=$(TSTBUILD)/testStateSort.o
SCFTSTOBJECT=$(TSTBUILD)/testSubspaceCF.o
DEFAULTSOBJECT=$(TSTBUILD)/defaultSystem.o
EXGENTSTOBJECT=$(TSTBUILD)/testExcitgen.o
TESTBASISEXEC=build/test/testBasis

$(DIRECTORIES):
	mkdir $@

$(DDIR)/%.d: $(SRC)/%.cxx $(DDIR)
	mkdir -p $(dir $@)
	$(CXX) $(EIGEN_FLAGS) $(LDFLAGS) -MM $< -o $@

$(BUILD)/%.o: $(SRC)/%.cxx
	mkdir -p $(dir $@)
	$(CXX) $(EIGEN_FLAGS) $(LDFLAGS) -c src/$*.cxx -o $@

$(TSTDDIR)/%.d: $(TST)/%.cxx $(TSTDDIR)
	mkdir -p $(dir $@)
	$(CXX) $(EIGEN_FLAGS) $(LDFLAGS) -MM $< -o $@

$(TSTBUILD)/%.o: $(TST)/%.cxx
	mkdir -p $(dir $@)
	$(CXX) $(EIGEN_FLAGS) $(LDFLAGS) -c $< -o $@


-include $(DEPENDENCIES)

testBasis: $(OBJECTS) $(BTESTOBJECT) 
	$(CXX) $^  $(LINKFLAGS) -o $(TSTBUILD)/$@ 

testAlg: $(OBJECTS) $(ALGTESTOBJECT)	
	$(CXX) $^  $(LINKFLAGS) -o $(TSTBUILD)/$@ 

testMetropolis: $(OBJECTS) $(MARKOVTESTOBJECT) $(DEFAULTSOBJECT)
	$(CXX)  $^ $(LINKFLAGS) -o $(TSTBUILD)/$@ 

testNnw: $(OBJECTS) $(NNWTESTOBJECT)
	$(CXX) $^ $(LINKFLAGS) -o $(TSTBUILD)/$@ 

testRBM: $(OBJECTS) $(RBMTESTOBJECT)
	$(CXX) $^ $(LINKFLAGS) -o $(TSTBUILD)/$@ 

testDirect: $(OBJECTS) $(DIRECTTSTOBJECT) $(DEFAULTSOBJECT)
	$(CXX) $^ $(LINKFLAGS) -o $(TSTBUILD)/$@

testSampler: $(OBJECTS) $(STESTOBJECT)
	$(CXX) $^ $(LINKFLAGS) -o $(TSTBUILD)/$@ 

testEigen: $(OBJECTS) $(ETESTOBJECT)
	$(CXX) $^ $(LINKFLAGS) -o $(TSTBUILD)/$@ 

testHam: $(OBJECTS) $(HTESTOBJECT)
	$(CXX)  $^ $(LINKFLAGS) -o $(TSTBUILD)/$@ 
testCF: $(OBJECTS) $(CFTESTOBJECT)
	$(CXX) $^ $(LINKFLAGS) -o $(TSTBUILD)/$@
testPT: $(OBJECTS) $(PTTESTOBJECT)
	$(CXX) $^  $(LINKFLAGS) -o $(TSTBUILD)/$@

testAbInitioHam: $(OBJECTS) $(ABINHAMTESTOBJECT)	
	$(CXX) $^ $(LINKFLAGS) -o $(TSTBUILD)/$@ 
testRandom: $(RNGTESTOBJECT)
	$(CXX)  $(LDFLAGS) $^ -o $(TSTBUILD)/$@
testSort: $(SORTTSTOBJECT) $(OBJECTS)
	$(CXX) $^ $(LINKFLAGS) -o $(TSTBUILD)/$@
testExGen: $(EXGENTSTOBJECT) $(OBJECTS) $(DEFAULTSOBJECT)
	$(CXX) $^ $(LINKFLAGS) -o $(TSTBUILD)/$@ 
testAlgAb: $(OBJECTS) $(ALGABTESTOBJECT)	
	$(CXX) $^  $(LINKFLAGS) -o $(TSTBUILD)/$@ 
clean:
	rm -r $(BUILD)/*                          
