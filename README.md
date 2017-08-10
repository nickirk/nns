##Project goal

1. Exploring the possibility to use neural network to represent wavefunction and compute ground state energy in ab-initio systems;
2. Targeting strongly correlated systems (static correlations)

##Methods

Based on stochastic gradient descend, stochastically choose determinants to sample the wave function. However the most important ones are always included in each sampling, we only choose some other determinants stochastically. 
0. Constructing a three-layer neural network, with input layer, one hiden layer and the output layer.
The input layer are for the basis determinants and the output layer will give output of the coefficient of this determinant.
1. Starting from HF determinant, input it to the NNW and get the coefficient $C_0$.
2. Choosing randomly the next determinant $i$ which is connected to the previous determinant based on the magnitude of the off-diagonal term. Take this determinant $i$ as the input to the NNW and obtain the coefficient $C_i$.
3. Repeat step 2 until we obtain $N$ coefficients of $N$ determinants.
4. Calculating 
	$$E({W_i})^k=<\Psi^k|H|\Psi^k>=\sum_{i}^N\sum_{j}C_i^*C_j<i|H|j>$$
and the derivatives $\delta_i^k=\frac{\partial E}{\partial W_i}$ using back-propagation algorithm. 
5. Going back to 1. and repeat the steps for M times (or until $\delta_i^k< \epsilon$ where $\epsilon$ is the tolerance).

##Model test

