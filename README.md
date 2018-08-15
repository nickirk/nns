## Project goal

1. Exploring the possibility to use neural network to represent wavefunction and compute ground state energy in ab-initio systems;
2. Targeting strongly correlated systems (static correlations)
3. Exploring sampling scheme similar to DMC.

## Methods


Given the number of electrons and number of orbitals, constructing all the
 possible determinants and passing the total number of determinants as the
size of the Hamiltonian.

0. Constructing a three-layer neural network, with input layer, one hidden layer and the output layer.
The input layer are for the basis determinants and the output layer will give output of the coefficient of this determinant. The output layer has two neurons, giving the real and imaginary part of the coefficent, respectively.

1. Starting with a random determinant. Choosing
 randomly (with probability propotional to the offdiagonal terms'
magnitude and properly unbias it) an connection, obtaining the magnitudes
for the seed and its connection from the neural network. Decide to move to
 the connection from the seed with probability $`|\psi_i|^2/|\psi_{i-1}|^2`$. Repeat it for $`N`$ steps.
2. Store all the coefficents for each determiant obtained from last step and the inputs and activations of each neuron for each dterminant. (This will be the memory bottle neck).

3. Calculating
```math
E(\{W_i\})^{(t)}=\frac{<\Psi^{(t)}|H|\Psi^{(t)}>}{<\Psi^{(t)}|\Psi^{(t)}>}=\frac{\sum_{i}^N\sum_{j'}C_i^*C_j<i|H|j>}{\sum_i|C_i|^2}
=\left<C_jH_{ij}/C_i^*\right>_M
```
where the prime on the $`j`$ index means sum over the determinants which are connected to each $`i`$ determinant in the list. In the futuer, we should move to a stochastic version to sum over the connected determinants. In the last step, we adopt the Metropolis sampling and calculate the expectation value via the average of a Markov chain process.

4. The derivatives

```math
\delta_i^{(t)}=\frac{\partial E}{\partial W_i^{(t)}}
```

using back-propagation algorithm. $`t`$ is the training number.
The derivatives should also be computed from the Metropolis sampling. Note that the above derivative of energy wrt parameter should be a real number.
We also notice that we can write

```math
E(\{W_i\})=E(\{C_i(\{W_k\})\})
```

so we can write the derivative as a sum

```math
\frac{\partial E}{\partial W_i}=\sum_j\left(\frac{\partial E}{\partial C^*_j}\frac{\partial C_j^*}{\partial W_i}+\frac{\partial E}{\partial C_j}\frac{\partial C_j}{\partial W_i}\right)
=\sum_j\left(2Re(\frac{\partial E}{\partial C_j^*})\frac{\partial Re(C_j^*)}{\partial W_i}+2 Imag(\frac{\partial E}{\partial C_j^*}\frac{\partial Imag(C_j^*)}{\partial W_i})\right)
```

where

```math 
\frac{\partial E}{\partial C_j^*}=\frac{\partial <\Psi|H|\Psi>/<\Psi|\Psi>}{\partial C_j^*}
=\frac{\sum_iC_iH_{ji}-EC_j}{\sum_i|C_i|^2}
=\left<\frac{H_{ji}-E\delta_{ij}}{C^*_i}\right>_M
```

5. Updating the connection weights between neurons $`W_i^{(t+1)}=W_i^{(t)}-\alpha \delta_i^{(t)}`$. (More solvers can be used instead of simple gradient descend. Will add more in the future).

6. Going back to 1 and repeat the steps for $`M (t=0,1,2,...,M)`$ times (or until $`\delta_i^{(t)}<\epsilon`$ where $`\epsilon`$ is the tolerance). Or simple the energy converges.


## Naming of Variables

This is borrowed from CC4S code to make sure that when different contributors contribute, we can stay in more or less the same style.
- Class names are in CamelCase, starting with an upper case character. They should have noun character. (AmplitudeSolver, Algorithm, not Solve),
- Derived classes should end with the name of the super class, unless it is kind of obvious (BinaryFtodLoader)
- Methods are in camelCase, starting with a lower case character. They should have verb character. (solve, normalize, not loader)
- Instances of classes and types are named in camelCase, preferably according to the class name they are an instance of. (Parser parser, int maxIterations)
- Mathematical objects, such as numbers and tensors, may additionally be named according to the  ematical symbol nomenclature (int n, Tensor<> Vabij).
- Acronyms should be treated as a normal word. (FtodLoader, not FTODLoader)

## Doing test

The test files are under the test directory. In principle, each implemented functionality should have a test file and pass the test.

The library and the tests can be built with cmake, e.g.
```
mkdir build
cd build
cmake ../
```
Will set up the make environment, then, the single tests can be built with (e.g.)
```
make testAlg
```
The executables can be run from the test/ subdirectory of the build
directory. Note that those tests that need to read in a user-supplied
Hamiltonian require an FCIDUMP file (example files are supplied with
this program in the run/ directory).

The default build type is RELEASE and will use extensive compiler
optimization, to build in a debug mode, specify
```
-DCMAKE_BUILD_TYPE=Debug
```
when running cmake.
