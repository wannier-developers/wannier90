# Methodology

`wannier90` computes maximally-localised Wannier functions (MLWF)
following the method of Marzari and Vanderbilt (MV) [@marzari-prb97].
For entangled energy bands, the method of Souza, Marzari and Vanderbilt
(SMV) [@souza-prb01] is used. We introduce briefly the methods and key
definitions here, but full details can be found in the original papers
and in Ref. [@mostofi-cpc08].

First-principles codes typically solve the electronic structure of
periodic materials in terms of Bloch states, $\psi_{n{\bf k}}$. These
extended states are characterised by a band index $n$ and crystal
momentum ${\bf k}$. An alternative representation can be given in terms
of spatially localised functions known as Wannier functions (WF). The WF
centred on a lattice site ${\bf R}$, $w_{n{\bf R}}({\bf r})$, is written
in terms of the set of Bloch states as

$$
w_{n{\bf R}}({\bf r})=\frac{V}{(2\pi)^3}\int_{\mathrm{BZ}}
\left[\sum_{m} U^{({\bf k})}_{mn} \psi_{m{\bf k}}({\bf
    r})\right]e^{-\mathrm{i}{\bf k}.{\bf R}} \:\mathrm{d}{\bf k} \ ,
$$

where $V$ is the unit cell volume, the integral is over the Brillouin
zone (BZ), and $\mathbf{U}^{(\mathbf{k})}$ is a unitary matrix that
mixes the Bloch states at each ${\bf k}$. $\mathbf{U}^{(\mathbf{k})}$ is
not uniquely defined and different choices will lead to WF with varying
spatial localisations. We define the spread $\Omega$ of the WF as

$$
\Omega=\sum_n \left[\langle w_{n{\bf 0}}({\bf r})| r^2 | w_{n{\bf
      0}}({\bf r}) \rangle - | \langle w_{n{\bf 0}}({\bf r})| {\bf r}
      | w_{n{\bf 0}}({\bf r}) \rangle |^2 \right].
$$

The total spread
can be decomposed into a gauge invariant term $\Omega_{\rm I}$ plus a
term ${\tilde \Omega}$ that is dependant on the gauge choice
$\mathbf{U}^{(\mathbf{k})}$. ${\tilde \Omega}$ can be further divided
into terms diagonal and off-diagonal in the WF basis, $\Omega_{\rm D}$
and $\Omega_{\rm OD}$,

$$
\Omega=\Omega_{\rm I}+{\tilde \Omega}=\Omega_{\rm I}+\Omega_{\rm
  D}+\Omega_{\rm OD}
$$

where

$$
\Omega_{{\rm I}}=\sum_n \left[\langle w_{n{\bf 0}}({\bf r})| r^2 | w_{n{\bf
      0}}({\bf r}) \rangle - \sum_{{\bf R}m} \left| \langle w_{m{\bf
      R}}({\bf r})| {\bf r} | w_{n{\bf 0}}({\bf r}) \rangle \right| ^2
      \right]
$$

$$
\Omega_{\rm D}=\sum_n \sum_{{\bf R}\neq{\bf 0}} |\langle w_{n{\bf
    R}}({\bf r})| {\bf r} | w_{n{\bf 0}}({\bf r}) \rangle|^2
$$

$$
\Omega_{\rm OD}=\sum_{m\neq n} \sum_{{\bf R}} |\langle w_{m{\bf R}}({\bf
  r})| {\bf r} | w_{n{\bf 0}}({\bf r}) \rangle |^2
$$

The MV method
minimises the gauge dependent spread $\tilde{\Omega}$ with respect the
set of $\mathbf{U}^{(\mathbf{k})}$ to obtain MLWF.

`wannier90` requires two ingredients from an initial electronic
structure calculation.

1.  The overlaps between the cell periodic part of the Bloch states
    $|u_{n{\bf k}}\rangle$
    
    $$
    M_{mn}^{(\bf{k,b})}=\langle u_{m{\bf k}}|u_{n{\bf k}+{\bf b}}\rangle,
    $$

    where the vectors ${\bf b}$, which connect a given k-point with its
    neighbours, are determined by `wannier90` according to the
    prescription outlined in Ref. [@marzari-prb97].

2.  As a starting guess the projection of the Bloch states
    $|\psi_{n\bf{k}}\rangle$ onto trial localised orbitals
    $|g_{n}\rangle$
    
    $$
    A_{mn}^{(\bf{k})}=\langle \psi_{m{\bf k}}|g_{n}\rangle,
    $$

Note that $\mathbf{M}^{(\mathbf{k},\mathbf{b})}$,
$\mathbf{A}^{(\mathbf{k})}$ and $\mathbf{U}^{(\mathbf{k})}$ are all
small, $N \times N$ matrices (see the following note) that are independent of the basis set
used to obtain the original Bloch states.

!!! note

    Technically, this is true for the case of an isolated group of $N$
    bands from which we obtain $N$ MLWF. When using the disentanglement
    procedure of Ref. [@souza-prb01], $\mathbf{A}^{(\mathbf{k})}$, for
    example, is a rectangular matrix. See
    Section [Entangled Energy Bands](#entangled-energy-bands).

To date, `wannier90` has been used in combination with electronic codes
based on plane-waves and pseudopotentials (norm-conserving and
ultrasoft [@vanderbilt-prb90]) as well as mixed basis set techniques
such as FLAPW [@posternak-prb02].

## Entangled Energy Bands

The above description is sufficient to obtain MLWF for an isolated set
of bands, such as the valence states in an insulator. In order to obtain
MLWF for entangled energy bands we use the "disentanglement" procedure
introduced in Ref. [@souza-prb01].

We define an energy window (the "outer window"). At a given k-point
$\bf{k}$, $N^{({\bf k})}_{{\rm win}}$ states lie within this energy
window. We obtain a set of $N$ Bloch states by performing a unitary
transformation amongst the Bloch states which fall within the energy
window at each k-point:

$$
| u_{n{\bf k}}^{{\rm opt}}\rangle = \sum_{m\in N^{({\bf k})}_{{\rm win}}}
U^{{\rm dis}({\bf k})}_{mn} | u_{m{\bf k}}\rangle
$$

where
$\bf{U}^{{\rm dis}({\bf k})}$ is a rectangular
$N^{({\bf k})}_{{\rm win}} \times N$ matrix (see the following note). The set of
$\bf{U}^{{\rm dis}({\bf k})}$ are obtained by minimising the gauge
invariant spread $\Omega_{{\rm I}}$ within the outer energy window. The
MV procedure can then be used to minimise $\tilde{\Omega}$ and hence
obtain MLWF for this optimal subspace.

!!! note
    As ${\bf U}^{{\rm dis}({\bf k})}$ is a rectangular matrix this is
    a unitary operation in the sense that $({\bf U}^{{\rm
     dis}({\bf k})})^{\dagger}{\bf U}^{{\rm dis}({\bf k})}={\bf 1}_N$.

It should be noted that the energy bands of this optimal subspace may
not correspond to any of the original energy bands (due to mixing
between states). In order to preserve exactly the properties of a system
in a given energy range (e.g., around the Fermi level) we introduce a
second energy window. States lying within this inner, or "frozen",
energy window are included unchanged in the optimal subspace.
