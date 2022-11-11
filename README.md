# Spin Operator in Quantum Chemistry
The spin-matrix element in AO basis is easily calculated:
$$
\langle \mu \alpha | \hat{S}_x | \nu \beta \rangle = S_{\mu\nu} \sigma_{\alpha\beta}^x
$$
in which $\sigma_{\alpha\beta}^x$ is the two-component spinor Pauli matrix. The second
quantization formula is given by,
$$
\hat{S}_x = \sigma^{x}_{\alpha\beta} S_{\mu \nu}
\hat{a}^{\dagger}_{\mu \alpha} \hat{a}_{\nu \beta} 
$$

The overal expectation value $\langle \hat{S}_x \rangle$ is given by, 
$$
\langle \hat{S}_x \rangle = S_{\mu\nu} \sigma_{\alpha\beta}^x 
P_{\mu \nu}^{\alpha \beta}
$$
in which $P_{\mu \nu}^{\alpha \beta}$ is the GHF density matrix,
$$
P_{\mu \nu}^{\alpha \beta} = \sum_{i} C_{\mu }^{i \alpha} C_{\nu}^{i \beta}
$$
$C_{\mu}^{i \alpha}$ is the $\alpha$-block of the $i$-th GHF MO.

# Constrained Constrained Density Functional Theory Formula
The spin expectation value can be regarded as some kind of Mulliken population,
$$
\langle S_x \rangle_A =  \sum_{\alpha \beta} \sigma_{\alpha\beta}^x \sum_{\mu \in A} \sum_{\nu} P_{\mu \nu}^{\alpha \beta} S_{\nu \mu} 
$$
if we write a Lagrangian as, (omitting some sum symbols)
$$
\mathcal{L} = E[\mathbf{P}] + \sum_{x} \lambda_x (
    \sigma_{\alpha\beta}^x  \sum_{\mu \in A} P_{\mu \nu}^{\alpha \beta} S_{\nu \mu} - S_x^A)
$$
By manipulating the target expectation value, we can rotate the spin on some atom
to the desired value, and observe the effect on the total energy.