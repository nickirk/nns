## Project goal

1. Exploring the possibility to use neural network to represent wavefunction and compute ground state energy in ab-initio systems;
2. Targeting strongly correlated systems (static correlations)

## Methods

Based on stochastic gradient descend, stochastically choose determinants to sample the wave function. However the most important ones are always included in each sampling (e.g. the HF determinant, if there are other determinants which are also important, we should always include them as well, static correlations), we only choose some other determinants stochastically. 
0. Given the number of electrons and number of orbitals, constructing all the possible determinants and passing the total number of determinants as the size of the Hamiltonian.
1. Generating a Hamiltonian 
```math
\hat{K}=\hat{H}-E_{HF}\delta_{ij}
```
of which the diagonal terms starting with $`0`$ and increase as $`0+diagDelta*i*i`$ and offDiag terms are
randomly chosen to be small nonzero values. 
2. Creating a list of determinants. Starting from HF determinant as the first element of the list, choosing randomly the next determinant $`i`$ which is connected to the previous determinant based on the magnitude of the off-diagonal term as the second determinant in the list. Then starting from determinant $`i`$, repeat the above procedure until we have $`N`$ determinants in the list. 
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
7. Going back to 1 and repeat the steps for $`M (t=0,1,2,\dots,M)`$ times (or until $`\delta_i^k<\epsilon`$ where $`\epsilon`$ is the tolerance).

## Structure of the code

```C++
class Hamiltonian {setting };
class Basis {constructing the basis};
class Sampler {sampling wavefunctions with determinants}
class NNW {neural network, back-propagation algorithm}
template 
getDets(){}
```

