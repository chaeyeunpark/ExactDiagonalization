# Simple exact diagoanlization library
This library present an implementation of exact diagonalization (ED) method for translational invariant spin-$1/2$ Hamiltonians defined in some simple lattices. 

## Exact diagonalization for one dimensional system

First, let us consider an one dimensional translational invariance Hamiltonian. 
We will consider $Z_2$ parity symmetry later.

We first divide the full Hamiltonian $H$ as translational invariance components. Thus, we write

$$ H^{\alpha} = \sum_{j=0}^{N-1} h_j^{\alpha}$$

where we use $\alpha$ to group the Hamiltonian components.
Typically, we use each $\alpha$ to represent Pauli strings in the Hamiltonian.

For instance, the transverse field Ising model 

$$H = \sum_{j=0}^{N-1} \sigma^j_z \sigma^{j+1}_z + \sigma^j_x$$

can be divided by $h^{1}_j = \sigma^j_z \sigma^{j+1}_z$ and $h^2_j = \sigma^j_x$.

In the qubit notation, the basis in the computational basis can be represented by a bitstring 
$\sigma \in \mathbb{Z}_{2}^{N} $. 
For each $\sigma$, we define a translational invariance state with a given momentum $k \in [0,\cdots,N-1]$ as

$$|\sigma(k)\rangle = \frac{1}{\sqrt{N_\sigma}} \sum_{j=0}^{N-1} e^{2 \pi i jk/N} T^i|\sigma\rangle$$

where $T$ is the translation operator and $\mathcal{N_\sigma}$ defines the normalization.

As $|T^a \sigma(k) \rangle = e^{-2\pi i ak/N} | \sigma(k) \rangle$
, we only need one $\sigma$ among its translationally shifted versions. We thus define _representative_ as a configuration that has minimum value (when considered as a binary represenation of an integer) among $\\{\sigma, T\sigma, T^2 \sigma, \cdots \\}$.

We note that the normalization factor $N_\sigma$ is important for the calculation of the Hamiltonian.
In our implementation, we caculate the $N_\sigma$ using $R_\sigma \in \mathbb{N}^+$ which is the minimum positive integer such that $T^{R\sigma}|\sigma\rangle = | \sigma \rangle$.
Such $R_\sigma$ should satisfy $N | R_\sigma$ ($N$ is a multiple of $R_\sigma$; which follows from the division theorem).
Then we have

$$|\sigma(k)\rangle = \frac{1}{\sqrt{N_\sigma}} \sum_{i=0}^{N-1} e^{2 \pi i k/N} T^i|\sigma\rangle 
					= \frac{N/R_\sigma}{\sqrt{N_\sigma}} \sum_{i=0}^{R_\sigma-1} e^{2 \pi i k/N} T^i|\sigma\rangle $$

and $1 = \langle \sigma(k) | \sigma(k) \rangle$
yields

$$ N_\sigma = N^2/R_\sigma.$$




