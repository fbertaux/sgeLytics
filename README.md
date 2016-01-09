# sgeLytics

Computational biology related. A small python module providing useful analytical functions for a (popular) stochastic gene expression model.

## Dependencies

Requires numpy and scipy.

## Model description

GeneON <-> GeneOFF (kon,koff)
GeneON -> mRNA + GeneON (ksm)
mRNA -> 0 (rm)
mRNA -> Prot + mRNA (ksp)
Prot -> 0 (rp)

## Equivalent parameterization

Giving the 6 parameters (kon, koff, ksm, rm, ksp, rp) is equivalent to give the 6 quantities (Ton, Toff, EM, HLm, EP, HLp) which are respectively the mean duration the gene stays ON, the mean duration the gene stays OFF, the mean mRNA level (of the steady-state distribution), the mRNA half-life, the mean protein level (of the steady-state distribution) and the protein half-life.