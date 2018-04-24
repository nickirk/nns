CXX=g++
LDFLAGS=-std=c++11 -g -Wall
EIGEN_PATH=lib/eigen/
EIGEN_FLAGS=-I$(EIGEN_PATH)
SOURCEFILES=Sampler.cxx Nnw.cxx Basis.cxx Hamiltonian.cxx Determinant.cxx EnergyCF.cxx NormCF.cxx EnergyEstimator.cxx FermionicHamiltonian.cxx AbInitioHamiltonian.cxx BosonicHamiltonian.cxx Solver.cxx Trainer.cxx MetropolisSampler.cxx ListGen.cxx EnergyEsMarkov.cxx Layer.cxx DenseLayer.cxx ConvLayer.cxx InputLayer.cxx math/MathFunctions.cxx
TESTFILES=testNnw.cxx testSampler.cxx testBasis.cxx testEigen.cxx testAlg.cxx testCF.cxx testPreTrain.cxx testAbInitioHam.cxx testAlgAb.cxx testRandom.cxx testMarkov.cxx
SRC=src
TST=test
BUILD=build
TSTBUILD=build/test
DDIR=.depend
TSTDDIR=.depend/test
DIRECTORIES=$(DDIR) $(BUILD) $(TSTDDIR) $(TSTBUILD)
SOURCES=$(patsubst %,$(SRC)/%,$(SOURCEFILES)) $(patsubst %,$(TST)/%,$(TESTFILES))
OBJECTS=$(patsubst %.cxx,$(BUILD)/%.o,$(SOURCEFILES)) 
DEPENDENCIES=$(patsubst %.cxx,$(DDIR)/%.d,$(SOURCEFILES)) $(patsubst %.cxx,$(TSTDDIR)/%.d,$(TESTFILES))
BTEST=basisTest
NNWTEST=nnwTest
ALGTEST=algTest
MARKOVTEST=markovTest
STEST=samplerTest
HTEST=hamTest
ETEST=eigenTest
CFTEST=cfTest
PTTEST=ptTest
ABINHAMTEST=abinhamTest
ALGABTEST=algAbTest
BTESTOBJECT=$(TSTBUILD)/testBasis.o
NNWTESTOBJECT=$(TSTBUILD)/testNnw.o
MARKOVTESTOBJECT=$(TSTBUILD)/testMarkov.o
ALGTESTOBJECT=$(TSTBUILD)/testAlg.o
STESTOBJECT=$(TSTBUILD)/testSampler.o
ETESTOBJECT=$(TSTBUILD)/testEigen.o
HTESTOBJECT=$(TSTBUILD)/testHam.o
CFTESTOBJECT=$(TSTBUILD)/testCF.o
PTTESTOBJECT=$(TSTBUILD)/testPreTrain.o
ABINHAMTESTOBJECT=$(TSTBUILD)/testAbInitioHam.o
ALGABTESTOBJECT=$(TSTBUILD)/testAlgAb.o
RNGTESTOBJECT=$(TSTBUILD)/testRandom.o
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
	$(CXX)  $(LDLFLAGS) $^ -o $(TSTBUILD)/$@ 
	
testAlg: $(OBJECTS) $(ALGTESTOBJECT)	
	$(CXX)  $(LDLFLAGS) $^ -o $(TSTBUILD)/$@ 

testMarkov: $(OBJECTS) $(MARKOVTESTOBJECT)	
	$(CXX)  $(LDLFLAGS) $^ -o $(TSTBUILD)/$@ 

testNnw: $(OBJECTS) $(NNWTESTOBJECT)
	$(CXX)  $(LDLFLAGS) $^ -o $(TSTBUILD)/$@ 

testSampler: $(OBJECTS) $(STESTOBJECT)
	$(CXX)  $(LDLFLAGS) $^ -o $(TSTBUILD)/$@ 

testEigen: $(OBJECTS) $(ETESTOBJECT)
	$(CXX)  $(LDLFLAGS) $^ -o $(TSTBUILD)/$@ 

testHam: $(OBJECTS) $(HTESTOBJECT)
	$(CXX)  $(LDLFLAGS) $^ -o $(TSTBUILD)/$@ 
	
testCF: $(OBJECTS) $(CFTESTOBJECT)
	$(CXX)  $(LDLFLAGS) $^ -o $(TSTBUILD)/$@
	
testPT: $(OBJECTS) $(PTTESTOBJECT)
	$(CXX)  $(LDLFLAGS) $^ -o $(TSTBUILD)/$@

testAbInitioHam: $(OBJECTS) $(ABINHAMTESTOBJECT)	
	$(CXX)  $(LDLFLAGS) $^ -o $(TSTBUILD)/$@ 
testRandom: $(RNGTESTOBJECT)
	$(CXX)  $(LDLFLAGS) $^ -o $(TSTBUILD)/$@
                                            
testAlgAb: $(OBJECTS) $(ALGABTESTOBJECT)	
	$(CXX)  $(LDLFLAGS) $^ -o $(TSTBUILD)/$@ 
clean:
	rm -r $(BUILD)/*                          
                                           
