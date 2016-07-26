# sgeLytics

Computational biology related.  A small python module providing useful
analytical functions for a (popular) stochastic gene expression model.

## Dependencies

Requires `numpy` and `scipy`.

## Model description

We consider a variant of the two-state transcription bursting model (see
Shahrezaei et al., 2008 and many others...) in which translation and
protein degradation reaction are assumed to be deterministic (see Paszek
et al., 2007).

This assumption is suited when the average protein level is large, which
is the case for most proteins in mammalian cells.

| `GeneON <-> GeneOFF` (kon, koff) | stochastic    |
| `GeneON -> mRNA + GeneON` (ksm)  | stochastic    |
| `mRNA -> 0` (rm)                 | stochastic    |
| `mRNA -> Prot + mRNA` (ksp)      | deterministic |
| `Prot -> 0` (rp)                 | deterministic |

Analytical results on the corresponding steady-state distribution were
derived for example in (Paszek et al., 2007).  During my thesis, I also
derived an expression for the associated auto-correlation function of
protein level (given in Bertaux et al., 2014).

## Equivalent parameterizations

Giving the 6 "biochemical" parameters `(kon, koff, ksm, rm, ksp, rp)`
is equivalent to give the 6 quantities `(Ton, Toff, EM, HLm, EP, HLp)`
which are respectively the mean duration the gene stays ON, the mean
duration the gene stays OFF, the mean mRNA level (of the steady-state
distribution), the mRNA half-life, the mean protein level (of the
steady-state distribution) and the protein half-life.

`(kon=1/Toff, koff=1/Ton)` are also equivalent to `(EG=Ton/(Ton+Toff),
rg=kon+koff)` where `EG` is the average fraction of the time the gene is
ON and `rg` is rate constant associated with both ON and OFF switching.

Because the translation and protein degradation are deterministic, `ksp`
just sets the scale for protein levels, but their fluctuation properties
are independent of `ksp`.

## Finding parameters from protein level fluctuation properties

Because we were interested in finding parameterizations from protein
level fluctuations properties, we finally found that if it exists, a
single parameterization is compatible with a given protein level noise
(`CV`) and mixing time (half auto-correlation time) when the three
timescales `rg`, `rm` and `rp` are fixed.

See the function `defineModelFromCVTau_rg_rm_rp`.

## Examples of use

See `example.py`.
