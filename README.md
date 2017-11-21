## Project goal

1. Exploring the possibility to use neural network to represent wavefunction and compute ground state energy in ab-initio systems;
2. Targeting strongly correlated systems (static correlations)
3. Exploring sampling scheme similar to DMC.
4. Implement stochastic reconfiguration 
## Methods

Based on stochastic gradient descend, stochastically choose determinants 
to sample the wave function. However the most important ones are always 
included in each sampling (e.g. the HF determinant, if there are other 
determinants which are also important, we should always include them as 
well, static correlations), we only choose some other determinants 
stochastically. 

Given the number of electrons and number of orbitals, constructing all the
 possible determinants and passing the total number of determinants as the
size of the Hamiltonian.

1. Starting with a totaly random list of determinants. For each, choosing
 randomly (with probability propotional to the offdiagonal terms' 
magnitude and properly unbias it) an connection, obtaining the magnitudes
for the seed and its connection from the neural network. Decide to move to
 the connection from the seed with probability $`|\psi_i|^2/|\psi_seed|^2`$. 

3. Constructing a three-layer neural network, with input layer, one hidden layer and the output layer.
The input layer are for the basis determinants and the output layer will give output of the coefficient of this determinant.

4. Feeding each determinant $`i`$ in the list to the neural netwrok and obtaining the coefficient $`C_i`$.

5. Calculating 
```math
E({W_i})^{(t)}=<\Psi^{(t)}|H|\Psi^{(t)}>=\sum_{i}^N\sum_{j}C_i^*C_j<i|H|j>
```
and the derivatives 

```math
\delta_i^{(t)}=\frac{\partial E}{\partial W_i^{(t)}}
```
using back-propagation algorithm. $`t`$ is the training number.

6. Updating the connection weights between neurons $`W_i^{(t+1)}=W_i^{(t)}-\alpha \delta_i^{(t)}`$.

7. Going back to 1 and repeat the steps for $`M (t=0,1,2,...,M)`$ times (or until $`\delta_i^{(t)}<\epsilon`$ where $`\epsilon`$ is the tolerance).

## Structure of the code

```C++
class Hamiltonian {setting };
class Basis {constructing the basis};
class Sampler {sampling wavefunctions with determinants}
class NNW {neural network, back-propagation algorithm}
template 
getDets(){}
```

## Naming of Variables

- Class names are in CamelCase, starting with an upper case character. They should have noun character. (AmplitudeSolver, Algorithm, not Solve),
- Derived classes should end with the name of the super class, unless it is kind of obvious (BinaryFtodLoader)
- Methods are in camelCase, starting with a lower case character. They should have verb character. (solve, normalize, not loader)
- Instances of classes and types are named in camelCase, preferably according to the class name they are an instance of. (Parser parser, int maxIterations)
- Mathematical objects, such as numbers and tensors, may additionally be named according to the mathematical symbol nomenclature (int n, Tensor<> Vabij).
- Acronyms should be treated as a normal word. (FtodLoader, not FTODLoader)

