# Electronic transport calculations with the `BoltzWann` module

By setting `boltzwann=TRUE`, `postw90` will call the
`BoltzWann` routines to calculate some transport coefficients using the
Boltzmann transport equation in the relaxation time approximation.

In particular, the transport coefficients that are calculated are: the
electrical conductivity $\mathrm{\bm{\sigma}}$, the Seebeck coefficient
$\mathrm{\bm{S}}$ and the coefficient $\mathrm{\bm{K}}$ (defined below;
it is the main ingredient of the thermal conductivity).

The list of parameters of the `BoltzWann` module are summarized in
<!-- TODO: link the table in sec:wannier -->
Table [\[parameter_keywords_bw\]](#parameter_keywords_bw){reference-type="ref"
reference="parameter_keywords_bw"}. 
An example of a Boltzmann transport
calculation can be found in the `wannier90` Tutorial.

!!!Note
    By default, the code assumes to be working with a 3D bulk
    material, with periodicity along all three spatial directions. If you
    are interested in studying 2D systems, set the correct value for the
    `boltz_2d_dir` variable 
    <!-- TODO: linke with sec in wannier -->
    (see Sec. [\[sec:boltz2ddir\]](#sec:boltz2ddir){reference-type="ref"
    reference="sec:boltz2ddir"} for the documentation). 
    This is important
    for the evaluation of the Seebeck coefficient.

Please cite the following paper [@pizzi-cpc14] when publishing results
obtained using the `BoltzWann` module:

> G. Pizzi, D. Volja, B. Kozinsky, M. Fornari, and N. Marzari,
> *BoltzWann: A code for the evaluation of thermoelectric and electronic
> transport properties with a maximally-localized Wannier functions
> basis*,
> Comp. Phys. Comm. 185, 422 (2014), DOI:10.1016/j.cpc.2013.09.015.

### Theory {#sec:boltzwann-theory}

The theory of the electronic transport using the Boltzmann transport
equations can be found for instance in
Refs. [@ziman-book72; @grosso-book00; @mahan-itc06]. Here we briefly
summarize only the main results.

The current density $\mathrm{\bm{J}}$ and the heat current (or energy
flux density) $\mathrm{\bm{J}}_Q$ can be written, respectively, as

$$
\begin{equation}
\begin{aligned}
  \mathrm{\bm{J}}   &= \mathrm{\bm{\sigma}}(\mathrm{\bm{E}} - \mathrm{\bm{S}} \mathrm{\bm{\nabla }}T) \\
  \mathrm{\bm{J}}_Q &= T \mathrm{\bm{\sigma }}\mathrm{\bm{S}} \mathrm{\bm{E}} - \mathrm{\bm{K}} \mathrm{\bm{\nabla }}T,
\end{aligned}
\end{equation}
$$ 

where the electrical conductivity
$\mathrm{\bm{\sigma}}$, the Seebeck coefficient $\mathrm{\bm{S}}$ and
$\mathrm{\bm{K}}$ are $3\times 3$ tensors, in general.

!!! note
    the thermal conductivity $\mathrm{\bm{\kappa}}$ (actually, the
    electronic part of the thermal conductivity), which is defined as the
    heat current per unit of temperature gradient in open-circuit
    experiments (i.e., with $\mathrm{\bm{J}}=0$) is not precisely
    $\mathrm{\bm{K}}$, but
    $\mathrm{\bm{\kappa }}= \mathrm{\bm{K}}-\mathrm{\bm{S}} \mathrm{\bm{\sigma }}\mathrm{\bm{S}} T$
    (see for instance Eq. (7.89) of Ref. [@ziman-book72] or Eq. (XI-57b) of
    Ref. [@grosso-book00]). The thermal conductivity $\mathrm{\bm{\kappa}}$
    can be then calculated from the $\mathrm{\bm{\sigma}}$,
    $\mathrm{\bm{S}}$ and $\mathrm{\bm{K}}$ tensors output by the code.

These quantities depend on the value of the chemical potential $\mu$ and
on the temperature $T$, and can be calculated as follows:

$$
\begin{equation}
\label{eq:boltz-sigmas}
\begin{aligned}
  \mathrm{[\bm{\sigma}]}_{ij}(\mu,T)&=e^2 \int_{-\infty}^{+\infty} d\varepsilon \left(-\frac {\partial f(\varepsilon,\mu,T)}{\partial \varepsilon}\right)\Sigma_{ij}(\varepsilon), \\
  [\mathrm{\bm{\sigma }}\mathrm{\bm{S}}]_{ij}(\mu,T)&=\frac e T \int_{-\infty}^{+\infty} d\varepsilon \left(-\frac {\partial f(\varepsilon,\mu,T)}{\partial \varepsilon}\right)(\varepsilon-\mu)\Sigma_{ij}(\varepsilon),\\
  [\mathrm{\bm{K}}]_{ij}(\mu,T)&=\frac 1 T \int_{-\infty}^{+\infty} d\varepsilon \left(-\frac {\partial f(\varepsilon,\mu,T)}{\partial \varepsilon}\right)(\varepsilon-\mu)^2 \Sigma_{ij}(\varepsilon),
\end{aligned}
\end{equation}
$$ 

where $[\mathrm{\bm{\sigma }}\mathrm{\bm{S}}]$ denotes
the product of the two tensors $\mathrm{\bm{\sigma}}$ and
$\mathrm{\bm{S}}$, $f(\varepsilon,\mu,T)$ is the usual Fermi--Dirac
distribution function

$$
\begin{equation}
f(\varepsilon,\mu,T) = \frac{1}{e^{(\varepsilon-\mu)/K_B T}+1}
\end{equation}
$$ 

and
$\Sigma_{ij}(\varepsilon)$ is the Transport Distribution Function (TDF)
tensor, defined as

$$
\begin{equation}
\Sigma_{ij}(\varepsilon) = \frac 1 V \sum_{n,\mathrm{\bm{k}}} v_i(n,\mathrm{\bm{k}}) v_j(n,\mathrm{\bm{k}}) \tau(n,\mathrm{\bm{k}}) \delta(\varepsilon - E_{n,k}).
\end{equation}
$$

In the above formula, the sum is over all bands $n$ and all states
$\mathrm{\bm{k}}$ (including spin, even if the spin index is not
explicitly written here). $E_{n,\mathrm{\bm{k}}}$ is the energy of the
$n-$th band at $\mathrm{\bm{k}}$, $v_i(n,\mathrm{\bm{k}})$ is the $i-$th
component of the band velocity at $(n,\mathrm{\bm{k}})$, $\delta$ is the
Dirac's delta function, $V = N_k \Omega_c$ is the total volume of the
system ($N_k$ and $\Omega_c$ being the number of $k$-points used to
sample the Brillouin zone and the unit cell volume, respectively), and
finally $\tau$ is the relaxation time. In the *relaxation-time
approximation* adopted here, $\tau$ is assumed as a constant, i.e., it
is independent of $n$ and $\mathrm{\bm{k}}$ and its value (in fs) is
read from the input variable `boltz_relax_time`.

### Files

#### `seedname_boltzdos.dat`

OUTPUT. Written by `postw90` if `boltz_calc_also_dos` is `true`. Note
that even if there are other general routines in `postw90` which
specifically calculate the DOS, it may be convenient to use the routines
in `BoltzWann` setting `boltz_calc_also_dos = true` if one must also
calculate the transport coefficients. In this way, the (time-demanding)
band interpolation on the $k$ mesh is performed only once, resulting in
a much shorter execution time.

The first lines are comments (starting with \# characters) which
describe the content of the file. Then, there is a line for each energy
$\varepsilon$ on the grid, containing a number of columns. The first
column is the energy $\varepsilon$. The following is the DOS at the
given energy $\varepsilon$. The DOS can either be calculated using the
adaptive smearing scheme (see the following note) if `boltz_dos_adpt_smr` is `true`, or using
a "standard" fixed smearing, whose type and value are defined by
`boltz_dos_smr_type` and `boltz_dos_smr_fixed_en_width`, respectively.
If spin decomposition is required (input flag `spin_decomp`), further
columns are printed, with the spin-up projection of the DOS, followed by
spin-down projection.

!!! note

    Note that in `BoltzWann` the adaptive (energy) smearing scheme
    also implements a simple adaptive $k-$mesh scheme: if at any given
    $k$ point one of the band gradients is zero, then that $k$ point is
    replaced by 8 neighboring $k$ points. Thus, the final results for
    the DOS may be slightly different with respect to that given by the
    `dos` module.

#### `seedname_tdf.dat`

OUTPUT. This file contains the Transport Distribution Function (TDF)
tensor $\mathrm{\bm{\Sigma}}$ on a grid of energies.

The first lines are comments (starting with \# characters) which
describe the content of the file. Then, there is a line for each energy
$\varepsilon$ on the grid, containing a number of columns. The first is
the energy $\varepsilon$, the followings are the components if
$\mathrm{\bm{\Sigma}}(\varepsilon)$ in the following order:
$\Sigma_{xx}$, $\Sigma_{xy}$, $\Sigma_{yy}$, $\Sigma_{xz}$,
$\Sigma_{yz}$, $\Sigma_{zz}$. If spin decomposition is required (input
flag `spin_decomp`), 12 further columns are provided, with the 6
components of $\mathrm{\bm{\Sigma}}$ for the spin up, followed by those
for the spin down.

The energy $\varepsilon$ is in eV, while $\mathrm{\bm{\Sigma}}$ is in
$\displaystyle\frac{1}{\hbar^2}\cdot{\mathrm{eV}\cdot\mathrm{fs}}\cdot
{\mathring{\mathrm{A}}^{-1}}$.

#### `seedname_elcond.dat`

OUTPUT. This file contains the electrical conductivity tensor
$\mathrm{\bm{\sigma}}$ on the grid of $T$ and $\mu$ points.

The first lines are comments (starting with \# characters) which
describe the content of the file. Then, there is a line for each
$(\mu,T)$ pair, containing 8 columns, which are respectively: $\mu$,
$T$, $\sigma_{xx}$, $\sigma_{xy}$, $\sigma_{yy}$, $\sigma_{xz}$,
$\sigma_{yz}$, $\sigma_{zz}$. (The tensor is symmetric).

The chemical potential is in eV, the temperature is in K, and the
components of the electrical conductivity tensor ar in SI units, i.e. in
1/$\Omega$/m.

#### `seedname_sigmas.dat`

OUTPUT. This file contains the tensor
$\mathrm{\bm{\sigma}}\mathrm{\bm{S}}$, i.e. the product of the
electrical conductivity tensor and of the Seebeck coefficient as defined
by Eq. $\eqref{eq:boltz-sigmas}$, on the grid of $T$ and $\mu$ points.

The first lines are comments (starting with \# characters) which
describe the content of the file. Then, there is a line for each
$(\mu,T)$ pair, containing 8 columns, which are respectively: $\mu$,
$T$, $(\sigma S)_{xx}$, $(\sigma S)_{xy}$, $(\sigma S)_{yy}$,
$(\sigma S)_{xz}$, $(\sigma S)_{yz}$, $(\sigma S)_{zz}$. (The tensor is
symmetric).

The chemical potential is in eV, the temperature is in K, and the
components of the tensor ar in SI units, i.e. in A/m/K.

#### `seedname_seebeck.dat`

OUTPUT. This file contains the Seebeck tensor $\mathrm{\bm{S}}$ on the
grid of $T$ and $\mu$ points.

Note that in the code the Seebeck coefficient is defined as zero when
the determinant of the electrical conductivity $\mathrm{\bm{\sigma}}$ is
zero. If there is at least one $(\mu, T)$ pair for which
$\det \mathrm{\bm{\sigma}}=0$, a warning is issued on the output file.

The first lines are comments (starting with \# characters) which
describe the content of the file. Then, there is a line for each
$(\mu,T)$ pair, containing 11 columns, which are respectively: $\mu$,
$T$, $S_{xx}$, $S_{xy}$, $S_{xz}$, $S_{yx}$, $S_{yy}$, $S_{yz}$,
$S_{zx}$, $S_{zy}$, $S_{zz}$.

NOTE: therefore, the format of the columns of this file is different
from the other three files (elcond, sigmas and kappa)!

The chemical potential is in eV, the temperature is in K, and the
components of the Seebeck tensor ar in SI units, i.e. in V/K.

#### `seedname_kappa.dat`

OUTPUT. This file contains the tensor $\mathrm{\bm{K}}$ defined in
Sec. [Theory](#sec:boltzwann-theory) on the grid of $T$ and $\mu$ points.

The first lines are comments (starting with \# characters) which
describe the content of the file. Then, there is a line for each
$(\mu,T)$ pair, containing 8 columns, which are respectively: $\mu$,
$T$, $K_{xx}$, $K_{xy}$, $K_{yy}$, $K_{xz}$, $K_{yz}$, $K_{zz}$. (The
tensor is symmetric).

The chemical potential is in eV, the temperature is in K, and the
components of the $\mathrm{\bm{K}}$ tensor are the SI units for the
thermal conductivity, i.e. in W/m/K.
