# `wannier90` as a library

This is a description of the interface between any external program and
the wannier code. There are two subroutines: `wannier_setup` and
`wannier_run`. Calling `wannier_setup` will return information required
to construct the $M_{mn}^{(\mathbf{k,b})}$ overlaps
(Ref. [@marzari-prb97], Eq. (25)) and
$A_{mn}^{(\mathbf{k})}=\left\langle
  \psi_{m\mathbf{k}}|g_{n}\right\rangle$ projections
(Ref. [@marzari-prb97], Eq. (62); Ref. [@souza-prb01], Eq. (22)). Once
the overlaps and projection have been computed, calling `wannier_run`
activates the minimisation and plotting routines in `wannier90`.

**IMPORTANT NOTE:** the library mode ONLY works in serial. Please call
it from a serial code, or if compiled in parallel, make sure to run it
from a single MPI process.

You can find a minimal example of how the library mode can be used among
the tests, in the file `test-suite/library-mode-test/test_library.F90`
in the Wannier90 git repository.

### Subroutines

#### `wannier_setup`

**`wannier_setup(seed_name,mp_grid,num_kpts,real_lattice,recip_lattice,`\
`              kpt_latt,num_bands_tot,num_atoms,atom_symbols,atoms_cart,`\
`              gamma_only,spinors,nntot,nnlist,nncell,num_bands,num_wann,proj_site,`\
`              proj_l,proj_m,proj_radial,proj_z,proj_x,proj_zona,`\
`              exclude_bands,proj_s,proj_s_qaxis)`**

Conditions:

-   $\verb#num_kpts# = \verb#mp_grid(1)# \times \verb#mp_grid(2)#
    \times \verb#mp_grid(3)#$.

-   $\verb#num_nnmax# = 12$

This subroutine returns the information required to determine the
required overlap elements $M_{mn}^{(\mathbf{k,b})}$ and projections
$A_{mn}^{(\mathbf{k})}$, i.e., `M_matrix` and `A_matrix`, described in
Section [6.1.2](#wannier_run){reference-type="ref"
reference="wannier_run"}.

For the avoidance of doubt, `real_lattice(1,2)` is the $y-$component of
the first lattice vector $\mathbf{A}_{1}$, etc.

The list of nearest neighbours of a particular k-point `nkp` is given by
`nnlist(nkp,1:nntot)`.

Additionally, the parameter `shell_list` may be specified in the
`wannier90` input file.

#### `wannier_run`

**`wannier_run(seed_name,mp_grid,num_kpts,real_lattice,recip_lattice,`\
`            kpt_latt,num_bands,num_wann,nntot,num_atoms,atom_symbols,`\
`            atoms_cart,gamma_only,M_matrix_orig,A_matrix,eigenvalues,`\
`            U_matrix,U_matrix_opt,lwindow,wann_centres,wann_spreads,`\
`            spread`)**

-   `character(len=*), intent(in) :: seed_name`\
    The seedname of the current calculation.

-   `integer, dimension(3), intent(in) :: mp_grid`\
    The dimensions of the Monkhorst-Pack k-point grid.

-   `integer, intent(in) :: num_kpts`\
    The number of k-points on the Monkhorst-Pack grid.

-   `real(kind=dp), dimension(3,3),` ` intent(in) :: real_lattice`\
    The lattice vectors in Cartesian co-ordinates in units of Angstrom.

-   `real(kind=dp), dimension(3,3), intent(in) :: recip_lattice`\
    The reciprical lattice vectors in Cartesian co-ordinates in units of
    inverse Angstrom.

-   `real(kind=dp), dimension(3,num_kpts),` ` intent(in) :: kpt_latt`\
    The positions of the k-points in fractional co-ordinates relative to
    the reciprocal lattice vectors.

-   `integer, intent(in) :: num_bands`\
    The total number of bands to be processed.

-   `integer, intent(in) :: num_wann`\
    The number of MLWF to be extracted.

-   `integer, intent(in) :: nntot`\
    The number of nearest neighbours for each k-point.

-   `integer, intent(in) :: num_atoms`\
    The total number of atoms in the system.

-   `character(len=20), dimension(num_atoms),`
    ` intent(in) :: atom_symbols`\
    The elemental symbols of the atoms.

-   `real(kind=dp), dimension(3,num_atoms),` `intent(in) :: atoms_cart`\
    The positions of the atoms in Cartesian co-ordinates in Angstrom.

-   `logical, intent(in) :: gamma_only`\
    Set to `.true.` if the underlying electronic structure calculation
    has been performed with only $\Gamma$-point sampling and, hence, if
    the Bloch eigenstates that are used to construct
    $A_{mn}^{(\mathbf{k})}$ and $M_{mn}^{\mathbf{(k,b)}}$ are real.

-   `complex(kind=dp),`
    ` dimension(num_bands,num_bands,nntot,num_kpts),`\
    `                  intent(in) :: M_matrix`\
    The matrices of overlaps between neighbouring periodic parts of the
    Bloch eigenstates at each k-point, $M_{mn}^{(\mathbf{(k,b)})}$
    (Ref. [@marzari-prb97], Eq. (25)).

-   `complex(kind=dp), dimension(num_bands,num_wann,num_kpts),`\
    `                  intent(in) :: A_matrix`\
    The matrices describing the projection of `num_wann` trial orbitals
    on `num_bands` Bloch states at each k-point, $A_{mn}^{(\mathbf{k})}$
    (Ref. [@marzari-prb97], Eq. (62); Ref. [@souza-prb01], Eq. (22)).

-   `real(kind=dp), dimension(num_bands,num_kpts),`
    `intent(in) :: eigenvalues`\
    The eigenvalues $\varepsilon_{n\mathbf{k}}$ corresponding to the
    eigenstates, in eV.

-   `complex(kind=dp), dimension(num_wann,num_wann,num_kpts),`\
    `                  intent(out) :: U_matrix`\
    The unitary matrices at each k-point (Ref. [@marzari-prb97],
    Eq. (59))

-   `complex(kind=dp), dimension(num_bands,num_wann,num_kpts),`\
    `               optional, intent(out) :: U_matrix_opt`\
    The unitary matrices that describe the optimal sub-space at each
    k-point (see Ref. [@souza-prb01], Section IIIa). The array is packed
    (see below)

-   `logical, dimension(num_bands,num_kpts), optional, intent(out) :: lwindow`\
    The element `lwindow(nband,nkpt)` is `.true.` if the band `nband`
    lies within the outer energy window at kpoint `nkpt`.

-   `real(kind=dp), dimension(3,num_wann), optional, intent(out) :: wann_centres`\
    The centres of the MLWF in Cartesian co-ordinates in Angstrom.

-   `real(kind=dp), dimension(num_wann), optional, intent(out) :: wann_spreads`\
    The spread of each MLWF in Å$^{2}$.

-   `real(kind=dp), dimension(3), optional, intent(out) ::` `spread`\
    The values of $\Omega$, $\Omega_{\mathrm{I}}$ and $\tilde{\Omega}$
    (Ref. [@marzari-prb97], Eq. (13)).

Conditions:

-   $\verb#num_wann# \le \verb#num_bands#$

-   $\verb#num_kpts# = \verb#mp_grid(1)# \times \verb#mp_grid(2)#
    \times \verb#mp_grid(3)#$.

If $\verb#num_bands# = \verb#num_wann#$ then `U_matrix_opt` is the
identity matrix and `lwindow=.true.`

For the avoidance of doubt, `real_lattice(1,2)` is the $y-$component of
the first lattice vector $\mathbf{A}_{1}$, etc.

$$\begin{aligned}
\verb#M_matrix(m,n,nn,nkp)# & = & \left\langle u_{m\mathbf{k}} |
u_{n\mathbf{k+b}}\right\rangle\\
\verb#A_matrix(m,n,nkp)# & = &
\left\langle \psi_{m\mathbf{k}}|g_{n}\right\rangle\\
\verb#eigenvalues(n,nkp)# &=& \varepsilon_{n\mathbf{k}}
\end{aligned}$$ where $$\begin{aligned}
\mathbf{k} &=&\verb#kpt_latt(1:3,nkp)#\\
\mathbf{k+b}&=& \verb#kpt_latt(1:3,nnlist(nkp,nn))# +
\verb#nncell(1:3,nkp,nn)# 
\end{aligned}$$ and $\left\{|g_{n}\rangle\right\}$ are a set of initial
trial orbitals. These are typically atom or bond-centred Gaussians that
are modulated by appropriate spherical harmonics.

Additional parameters should be specified in the `wannier90` input file.
