# Simple exact diagoanlization library
This library present exact diagonalization implementaions for translational invariant spin-$1/2$ Hamiltonians defined in some simple lattices. 

##Exact diagonalization for one dimensional system

For simplicity, we only consider one dimensional translational invariance. We will consider $Z_2$ parity symmetry later.

We divide the full Hamiltonian $H$ as translational invariance components. Thus, we write
$$ H^{\alpha} = \sum_{j=0}^{N-1} h_j^{\alpha}$$
where we use $\alpha$ to indicate how they map the configuration.
In typical cases, each $\alpha$ represents Pauli strings in the Hamiltonian.

For instance, the transverse field Ising model 

$$H = \sum_{j=0}^{N-1} \sigma^j_z \sigma^{j+1}_z + \sigma^j_x$$

can be divided by $h^{1}_j = \sigma^j_z \sigma^{j+1}_z$ and $h^2_j = \sigma^j_x$.

In a qubit notation, the basis in the computational basis can be represented by a bitstring 
$\sigma \in \mathbb{Z}_{2}^{N} $. 
For each $\sigma$, we can define a translational invariance state with momentum $k \in [0,\cdots,N-1]$ as

$$|\sigma(k)\rangle = \frac{1}{\sqrt{N_\sigma}} \sum_{i=0}^{N-1} e^{2 \pi i k/N} T^i|\sigma\rangle$$

where $T$ is the translation operator and $\mathcal{N_\sigma}$ defines the normalization.

Before calculating hamiltonian using this, we have to calculate the normalization first. 
In our implementation, we find $R_\sigma \in \mathbb{N}^+$ which is the minimum positive integer that $T^{R\sigma}|\sigma\rangle = | \sigma \rangle$.
Such $R_\sigma$ should satisfy $N | R_\sigma$ ($N$ is a multiple of $R_\sigma$; which follows from the division theorem).
Then we have

$$|\sigma(k)\rangle = \frac{1}{\sqrt{N_\sigma}} \sum_{i=0}^{N-1} e^{2 \pi i k/N} T^i|\sigma\rangle 
					= \frac{N/R_\sigma}{\sqrt{N_\sigma}} \sum_{i=0}^{R_\sigma-1} e^{2 \pi i k/N} T^i|\sigma\rangle $$

and $1 = \langle \sigma(k) | \sigma(k) \rangle$
yields

$$ N_\sigma = N^2/R_\sigma.$$




