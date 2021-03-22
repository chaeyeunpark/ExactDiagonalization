# Simple exact diagoanlization library
This library presents exact diagonalization (ED) implementaions for translational invariant spin-$1/2$ Hamiltonians defined in several simple lattices. 

## Exact diagonalization for one dimensional system


To deal with a symmetry in ED, we divide the full Hamiltonian $H$ in symmetric subspaces. 
We will explain how to construct this subspace briefly in this section.

For simplicity, we only consider translational invariant system in one dimension. The $Z_2$ parity symmetry will be discussed later.
First, we write the Hamiltonian as
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

For each $\sigma \in \mathbb{Z}_2^N$, we define representative of $\sigma$ as $\mathrm{rep}(\sigma) = \min\\{\sigma, T \sigma, T^2 \sigma, \cdots, T^{N-1}\sigma\\}$.
Then for all representatives $\mathrm{R} = \\{\mathrm{rep}(\sigma)\\}$ and given $k$,
$\\{|\sigma(k)\rangle | \sigma \in \mathrm{R}\\}$
will be basis vectors for the translational invariant subspace with momentum $k$.




