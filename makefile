CXX=g++
LDFLAGS=-std=c++11 -g -Wall
EIGEN_PATH=lib/eigen/
EIGEN_FLAGS=-I$(EIGEN_PATH)
SOURCEFILES=Sampler.cxx Nnw.cxx Basis.cxx Hamiltonian.cxx Determinant.cxx
TESTFILES=testNnw.cxx testSampler.cxx testBasis.cxx
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
STEST=samplerTest
HTEST=hamTest
ETEST=eigenTest
BTESTOBJECT=$(TSTBUILD)/testBasis.o
NNWTESTOBJECT=$(TSTBUILD)/testNnw.o
STESTOBJECT=$(TSTBUILD)/testSampler.o
ETESTOBJECT=$(TSTBUILD)/testEigen.o
HTESTOBJECT=$(TSTBUILD)/testHam.o

$(DIRECTORIES):
	mkdir $@

$(DDIR)/%.d: $(SRC)/%.cxx $(DDIR)
	$(CXX) $(EIGEN_FLAGS) $(LDFLAGS) -MM $< -o $@

$(BUILD)/%.o: $(SRC)/%.cxx
	$(CXX) $(EIGEN_FLAGS) $(LDFLAGS) -c src/$*.cxx -o $@

$(TSTDDIR)/%.d: $(TST)/%.cxx $(TSTDDIR)
	$(CXX) $(EIGEN_FLAGS) $(LDFLAGS) -MM $< -o $@

$(TSTBUILD)/%.o: $(TST)/%.cxx $(TSTBUILD)
	$(CXX) $(EIGEN_FLAGS) $(LDFLAGS) -c $< -o $@


-include $(DEPENDENCIES)

$(TSTDDIR)/testBasis: $(OBJECTS) $(BTESTOBJECT) 
	mkdir -p $(dir $@)
	$(CXX)  $(LDLFLAGS) $^ -o $@ 

$(TSTDDIR)/testNnw: $(OBJECTS) $(NNWTESTOBJECT)
	mkdir -p $(dir $@)
	$(CXX)  $(LDLFLAGS) $^ -o $@  

$(TSTDDIR)/testSampler: $(OBJECTS) $(STESTOBJECT)
	mkdir -p $(dir $@)
	$(CXX)  $(LDLFLAGS) $^ -o $@

$(TSTDDIR)/testEigen: $(OBJECTS) $(ETESTOBJECT)
	mkdir -p $(dir $@)
	$(CXX)  $(LDLFLAGS) $^ -o $@

$(TSTDDIR)/testHam: $(OBJECTS) $(HTESTOBJECT)
	mkdir -p $(dir $@)
	$(CXX)  $(LDLFLAGS)  $^ -o $@
                                            
                                            
                                           
