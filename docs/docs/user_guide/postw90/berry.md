# Overview of the `berry` module

The `berry` module of `postw90` is called by setting ` berry = true` and
choosing one or more of the available options for `berry_task`. The
routines in the `berry` module which compute the $k$-space Berry
curvature, orbital magnetization and spin Hall conductivity are also
called when `kpath = true` and `kpath_task = {curv,morb,shc}`, or when
`kslice = true` and `kslice_task = {curv,morb,shc}`.

### Background: Berry connection and curvature

The Berry connection is defined in terms of the cell-periodic Bloch
states $\vert u_{n{\bf k}}\rangle=e^{-i{\bf k}\cdot{\bf r}}\vert
\psi_{n{\bf k}}\rangle$ as
$${\bf A}_n({\bf k})=\langle u_{n{\bf k}}\vert i\bm{\nabla}_{\bf k}\vert
u_{n{\bf k}}\rangle,$$ and the Berry curvature is the curl of the
connection,
$$\bm{\Omega}_n({\bf k})=\bm{\nabla}_{\bf k}\times {\bf A}_n({\bf k})=
-{\rm Im}
\langle \bm{\nabla}_{\bf k} u_{n{\bf k}}\vert \times
\vert\bm{\nabla}_{\bf k} u_{n{\bf k}}\rangle.$$ These two quantities
play a central role in the description of several electronic properties
of crystals [@xiao-rmp10]. In the following we will work with a matrix
generalization of the Berry connection,
$${\bf A}_{nm}({\bf k})=\langle u_{n{\bf k}}\vert i\bm{\nabla}_{\bf k}\vert
u_{m{\bf k}}\rangle={\bf A}_{mn}^*({\bf k}),
\label{eq:berry-connection-matrix}$$ and write the curvature as an
antisymmetric tensor, $$\label{eq:curv}
\Omega_{n,\alpha\beta}({\bf k}) =\epsilon_{\alpha\beta\gamma}
\Omega_{n,\gamma}({\bf k})=-2{\rm Im}\langle 
\nabla_{k_\alpha} u_{n\bf k}\vert \nabla_{k_\beta} u_{n\bf k}\rangle.$$

### `berry_task=kubo`: optical conductivity and joint density of states 

The Kubo-Greenwood formula for the optical conductivity of a crystal in
the independent-particle approximation reads
$$\sigma_{\alpha\beta}(\hbar\omega)=\frac{ie^2\hbar}{N_k\Omega_c}
\sum_{\bf k}\sum_{n,m}
\frac{f_{m{\bf k}}-f_{n{\bf k}}}
     {\varepsilon_{m{\bf k}}-\varepsilon_{n{\bf k}}}
\frac{\langle\psi_{n{\bf k}}\vert v_\alpha\vert\psi_{m{\bf k}}\rangle
      \langle\psi_{m{\bf k}}\vert v_\beta\vert\psi_{n{\bf k}}\rangle}
{\varepsilon_{m{\bf k}}-\varepsilon_{n{\bf k}}-(\hbar\omega+i\eta)}.$$
Indices $\alpha,\beta$ denote Cartesian directions, $\Omega_c$ is the
cell volume, $N_k$ is the number of $k$-points used for sampling the
Brillouin zone, and $f_{n{\bf k}}=f(\varepsilon_{n{\bf k}})$ is the
Fermi-Dirac distribution function. $\hbar\omega$ is the optical
frequency, and $\eta>0$ is an adjustable smearing parameter with units
of energy.

The off-diagonal velocity matrix elements can be expressed in terms of
the connection matrix [@blount-ssp62], $$\label{eq:velocity_mat}
\langle\psi_{n{\bf k}}\vert {\bf v} \vert\psi_{m{\bf k}}\rangle=
-\frac{i}{\hbar}(\varepsilon_{m{\bf k}}-\varepsilon_{n{\bf k}})
{\bf A}_{nm}({\bf k})\,\,\,\,\,\,\,\,(m\not= n).$$ The conductivity
becomes $$\begin{aligned}
\label{eq:sig-bz}
\sigma_{\alpha\beta}(\hbar\omega)&=
\frac{1}{N_k}\sum_{\bf k}\sigma_{{\bf k},\alpha\beta}(\hbar\omega)\\
\label{eq:sig-k}
\sigma_{{\bf k},\alpha\beta}(\hbar\omega)&=\frac{ie^2}{\hbar\Omega_c}\sum_{n,m}
(f_{m{\bf k}}-f_{n{\bf k}})
\frac{\varepsilon_{m{\bf k}}-\varepsilon_{n{\bf k}}}
{\varepsilon_{m{\bf k}}-\varepsilon_{n{\bf k}}-(\hbar\omega+i\eta)}
A_{nm,\alpha}({\bf k})A_{mn,\beta}({\bf k}).
\end{aligned}$$

Let us decompose it into Hermitian (dissipative) and anti-Hermitean
(reactive) parts. Note that $$\label{eq:lorentzian}
\overline{\delta}(\varepsilon)=\frac{1}{\pi}{\rm Im}
\left[\frac{1}{\varepsilon-i\eta}\right],$$ where $\overline{\delta}$
denotes a "broadended" delta-function. Using this identity we find for
the Hermitean part $$\label{eq:sig-H}
\sigma_{{\bf k},\alpha\beta}^{\rm H}(\hbar\omega)=-\frac{\pi e^2}{\hbar\Omega_c}
\sum_{n,m}(f_{m{\bf k}}-f_{n{\bf k}})
(\varepsilon_{m{\bf k}}-\varepsilon_{n{\bf k}})
A_{nm,\alpha}({\bf k})A_{mn,\beta}({\bf k})
\overline{\delta}(\varepsilon_{m{\bf k}}-\varepsilon_{n{\bf k}}-\hbar\omega).$$
Improved numerical accuracy can be achieved by replacing the Lorentzian
([\[eq:lorentzian\]](#eq:lorentzian){reference-type="ref"
reference="eq:lorentzian"}) with a Gaussian, or other shapes. The
analytical form of $\overline{\delta}(\varepsilon)$ is controlled by the
keyword `[kubo_]smr_type`.

The anti-Hermitean part of
Eq. ([\[eq:sig-k\]](#eq:sig-k){reference-type="ref"
reference="eq:sig-k"}) is given by $$\label{eq:sig-AH}
\sigma_{{\bf k},\alpha\beta}^{\rm AH}(\hbar\omega)=\frac{ie^2}{\hbar\Omega_c}
\sum_{n,m}(f_{m{\bf k}}-f_{n{\bf k}})
{\rm Re}\left[ \frac{\varepsilon_{m{\bf k}}-\varepsilon_{n{\bf k}}}
                    {\varepsilon_{m{\bf k}}-\varepsilon_{n{\bf k}}
                     -(\hbar\omega+i\eta)}
\right]
A_{nm,\alpha}({\bf k})A_{mn,\beta}({\bf k}).$$ Finally the joint density
of states is $$\label{eq:jdos}
\rho_{cv}(\hbar\omega)=\frac{1}{N_k}\sum_{\bf k}\sum_{n,m}
f_{n{\bf k}}(1-f_{m{\bf k}})
\overline{\delta}(\varepsilon_{m{\bf k}}-\varepsilon_{n{\bf k}}-\hbar\omega).$$

Equations ([\[eq:lorentzian\]](#eq:lorentzian){reference-type="ref"
reference="eq:lorentzian"}--[\[eq:jdos\]](#eq:jdos){reference-type="ref"
reference="eq:jdos"}) contain the parameter $\eta$, whose value can be
chosen using the keyword\
` [kubo_]smr_fixed_en_width`. Better results can often be achieved by
adjusting the value of $\eta$ for each pair of states, i.e.,
$\eta\rightarrow \eta_{nm\bf k}$. This is done as follows (see
description of the keyword `adpt_smr_fac`)
$$\eta_{nm{\bf k}}=\alpha\vert \bm{\nabla}_{\bf k}
(\varepsilon_{m{\bf k}}-\varepsilon_{n{\bf k}})\vert \Delta k.$$

The energy eigenvalues $\varepsilon_{n\bf k}$, band velocities
$\bm{\nabla}_{\bf k}\varepsilon_{n{\bf k}}$, and off-diagonal Berry
connection ${\bf A}_{nm}({\bf k})$ entering the previous four equations
are evaluated over a $k$-point grid by Wannier interpolation, as
described in Refs. [@wang-prb06; @yates-prb07]. After averaging over the
Brillouin zone, the Hermitean and anti-Hermitean parts of the
conductivity are assembled into the symmetric and antisymmetric tensors
$$\begin{aligned}
\sigma^{\rm S}_{\alpha\beta}&=
{\rm Re}\sigma^{\rm H}_{\alpha\beta}+i{\rm Im}\sigma^{\rm AH}_{\alpha\beta}\\
\sigma^{\rm A}_{\alpha\beta}&=
{\rm Re}\sigma^{\rm AH}_{\alpha\beta}+i{\rm Im}\sigma^{\rm H}_{\alpha\beta},
\end{aligned}$$ whose independent components are written as a function
of frequency onto nine separate files.

### `berry_task=ahc`: anomalous Hall conductivity

The antisymmetric tensor $\sigma^{\rm A}_{\alpha\beta}$ is odd under
time reversal, and therefore vanishes in non-magnetic systems, while in
ferromagnets with spin-orbit coupling it is generally nonzero. The
imaginary part ${\rm Im}\sigma^{\rm H}_{\alpha\beta}$ describes magnetic
circular dichroism, and vanishes as $\omega\rightarrow
0$. The real part ${\rm Re}\sigma^{\rm AH}_{\alpha\beta}$ describes the
anomalous Hall conductivity (AHC), and remains finite in the static
limit.

The intrinsic dc AHC is obtained by setting $\eta=0$ and $\omega=0$ in
Eq. ([\[eq:sig-AH\]](#eq:sig-AH){reference-type="ref"
reference="eq:sig-AH"}). The contribution from point ${\bf k}$ in the
Brillouin zone is
$$\sigma^{\rm AH}_{{\bf k},\alpha\beta}(0)=\frac{2e^2}{\hbar\Omega_c}
\sum_{n,m}f_{n\bf k}(1-f_{m\bf k})
{\rm Im}\langle \nabla_{k_\alpha} u_{n\bf k}\vert u_{m\bf k}\rangle
\langle u_{m\bf k}\vert\nabla_{k_\beta} u_{n\bf k}\rangle,$$ where we
replaced $f_{n\bf k}-f_{m\bf k}$ with
$f_{n\bf k}(1-f_{m\bf k})-f_{m\bf k}(1-f_{n\bf k})$.

This expression is not the most convenient for *ab initio* calculations,
as the sums run over the complete set of occupied and empty states. In
practice the sum over empty states can be truncated, but a relatively
large number should be retained to obtain accurate results. Using the
resolution of the identity $1=\sum_m \vert u_{m\bf
  k}\rangle \langle u_{m\bf k}\vert$ and noting that the term
$\sum_{n,m}f_{n\bf k}f_{m\bf k}(\ldots)$ vanishes identically, we arrive
at the celebrated formula for the intrinsic AHC in terms of the Berry
curvature, $$\begin{aligned}
\label{eq:ahc}
\sigma^{\rm AH}_{\alpha\beta}(0)&=\frac{e^2}{\hbar}
\frac{1}{N_k\Omega_c}\sum_{\bf k}(-1)\Omega_{\alpha\beta}({\bf k}),\\
%\sum_n (-1)f_{n\bf k}\Omega_{n,\alpha\beta}({\bf k}).
\label{eq:curv-occ}
\Omega_{\alpha\beta}({\bf k})&=\sum_n f_{n\bf k}\Omega_{n,\alpha\beta}({\bf k}).
\end{aligned}$$ Note that only *occupied* states enter this expression.
Once we have a set of Wannier functions spanning the valence bands
(together with a few low-lying conduction bands, typically)
Eq. ([\[eq:ahc\]](#eq:ahc){reference-type="ref" reference="eq:ahc"}) can
be evaluated by Wannier interpolation as described in
Refs. [@wang-prb06; @lopez-prb12], with no truncation involved.

### `berry_task=morb`: orbital magnetization

The ground-state orbital magnetization of a crystal is given
by [@xiao-rmp10; @ceresoli-prb06] $$\begin{aligned}
\label{eq:morb}
{\bf M}^{\rm orb}&=\frac{e}{\hbar}
%\int_{\rm BZ}\frac{d{\bf k}}{(2\pi)^3}
\frac{1}{N_k\Omega_c}\sum_{\bf k}{\bf M}^{\rm orb}({\bf k}),\\
\label{eq:morb-k}
{\bf M}^{\rm orb}({\bf k})&=
\sum_n\,\frac{1}{2}f_{n{\bf k}}\,
{\rm Im}\,\langle \bm{\nabla}_{\bf k}u_{n{\bf k}}\vert
\times
\left(H_{\bf k}+\varepsilon_{n{\bf k}}-2\varepsilon_F\right)
\vert \bm{\nabla}_{\bf k}u_{n{\bf k}}\rangle,
\end{aligned}$$ where $\varepsilon_F$ is the Fermi energy. The
Wannier-interpolation calculation is described in Ref. [@lopez-prb12].
Note that the definition of ${\bf M}^{\rm orb}({\bf k})$ used here
differs by a factor of $-1/2$ from the one in Eq. (97) and Fig. 2 of
that work.

### `berry_task=shc`: spin Hall conductivity {#sec:shc}

The Kubo-Greenwood formula for the intrinsic spin Hall conductivity
(SHC) of a crystal in the independent-particle approximation reads
[@qiao-prb2018; @ryoo-prb2019; @guo-prl2008] $$\label{eq:kubo_shc}
\sigma_{\alpha\beta}^{\text{spin}\gamma}(\omega) =  \frac{\hbar}{\Omega_c N_k}
\sum_{\bm{k}}\sum_{n} f_{n\bm{k}} \\
\sum_{m \neq n}
\frac{2\operatorname{Im}[\langle n\bm{k}| \hat{j}_{\alpha}^{\gamma}|m\bm{k}\rangle
	\langle m\bm{k}| -e\hat{v}_{\beta}|n\bm{k}\rangle]}
{(\epsilon_{n\bm{k}}-\epsilon_{m\bm{k}})^2-(\hbar\omega +i\eta)^2}.$$
The spin current operator $\hat{j}_{\alpha}^{\gamma}=
\frac{1}{2}\{\hat{s}_{\gamma},\hat{v}_{\alpha}\}$ where the spin
operator $\hat{s}_{\gamma}=\frac{\hbar}{2}\hat{\sigma}_{\gamma}$.
Indices $\alpha,\beta$ denote Cartesian directions, $\gamma$ denotes the
direction of spin, commonly $\alpha = x, \beta = y, \gamma = z$.
$\Omega_c$ is the cell volume, $N_k$ is the number of $k$-points used
for sampling the Brillouin zone, and
$f_{n{\bf k}}=f(\varepsilon_{n{\bf k}})$ is the Fermi-Dirac distribution
function. $\hbar\omega$ is the optical frequency, and $\eta>0$ is an
adjustable smearing parameter with unit of energy.

The velocity matrix element in the numerator is the same as
Eq. ([\[eq:velocity_mat\]](#eq:velocity_mat){reference-type="ref"
reference="eq:velocity_mat"}), so the only unknown quantity is the spin
current matrix
$\langle n\bm{k}| \hat{j}_{\alpha}^{\gamma}|m\bm{k}\rangle$. We can use
Wannier interpolation technique to efficiently calculate this matrix,
and there are two derivation according to the degree of approximation. A
noteworthy difference is the way in which two *ab-initio* matrix
elements are evaluated,
$$\langle u_{n{\bf k}}\vert\sigma_\gamma H_{\bf k}\vert u_{m{\bf k}+{\bf b}}\rangle, \langle u_{n{\bf k}}\vert\sigma_\gamma \vert u_{m{\bf k}+{\bf b}}\rangle, \gamma = x, y, z$$
These are evaluated by `pw2wannier90` using Ryoo's method. In contrast,
Qiao's method does not require `pw2wannier90`, but it assumes an
approximation
$1\approx\sum_{ l\in ab-initio{\rm \,bands}}|u_{l\bm{k}}\rangle \langle u_{l\bm{k}}|$.
You can choose which method to evaluate this value with `shc_method` in
the input file. For a full derivation please refer to
Ref. [@qiao-prb2018] or Ref. [@ryoo-prb2019].

The Eq. ([\[eq:kubo_shc\]](#eq:kubo_shc){reference-type="ref"
reference="eq:kubo_shc"}) can be further separated into band-projected
Berry curvature-like term $$\label{eq:kubo_shc_berry}
\Omega_{n,\alpha\beta}^{\text{spin}\gamma}(\bm{k}) = {\hbar}^2 \sum_{
	m\ne n}\frac{-2\operatorname{Im}[\langle n\bm{k}| 
	\frac{1}{2}\{\hat{\sigma}_{\gamma},\hat{v}_{\alpha}\}|m\bm{k}\rangle
	\langle m\bm{k}| \hat{v}_{\beta}|n\bm{k}\rangle]}
{(\epsilon_{n\bm{k}}-\epsilon_{m\bm{k}})^2-(\hbar\omega+i\eta)^2},$$
$k$-resolved term which sums over occupied bands
$$\label{eq:kubo_shc_berry_sum}
\Omega_{\alpha\beta}^{\text{spin}\gamma}(\bm{k}) = \sum_{n}
f_{n\bm{k}} \Omega_{n,\alpha\beta}^{\text{spin}\gamma}(\bm{k}),$$ and
the SHC is $$\sigma_{\alpha\beta}^{\text{spin}\gamma}(\omega) = 
-\frac{e^2}{\hbar}\frac{1}{\Omega_c N_k}\sum_{\bm{k}}
\Omega_{\alpha\beta}^{\text{spin}\gamma}(\bm{k}).$$ The unit of the
$\Omega_{n,\alpha\beta}^{\text{spin}\gamma}(\bm{k})$ is
$\text{length}^{2}$ (Angstrom$^2$ or Bohr$^2$, depending on your choice
of `berry_curv_unit` in the input file), and the unit of
$\sigma_{\alpha\beta}^{\text{spin}\gamma}$ is $(\hbar/e)$S/cm (the unit
is written in the header of the output file). The case of $\omega=0$
corresponds to direct current (dc) SHC while that of $\omega\ne0$
corresponds to alternating current (ac) SHC or frequency-dependent SHC.
Note in some papers
Eq. ([\[eq:kubo_shc_berry\]](#eq:kubo_shc_berry){reference-type="ref"
reference="eq:kubo_shc_berry"}) is called as spin Berry curvature.
However, it was pointed out by Ref. [@Gradhand_2012] that this name is
misleading, so we use a somewhat awkward name "Berry curvature-like
term" to refer to
Eq. ([\[eq:kubo_shc_berry\]](#eq:kubo_shc_berry){reference-type="ref"
reference="eq:kubo_shc_berry"}). The $k$-resolved term
Eq. ([\[eq:kubo_shc_berry_sum\]](#eq:kubo_shc_berry_sum){reference-type="ref"
reference="eq:kubo_shc_berry_sum"}) can be used to draw `kslice` plot,
and the band-projected Berry curvature-like term
Eq. ([\[eq:kubo_shc_berry\]](#eq:kubo_shc_berry){reference-type="ref"
reference="eq:kubo_shc_berry"}) can be used to color the `kpath` plot.

Same as the case of optical conductivity, the parameter $\eta$ contained
in the
Eq. ([\[eq:kubo_shc_berry\]](#eq:kubo_shc_berry){reference-type="ref"
reference="eq:kubo_shc_berry"}) can be chosen using the keyword
`[kubo_]smr_fixed_en_width`. Also, adaptive smearing can be employed by
the keyword `[kubo_]adpt_smr` (see Examples 29 and 30 in the Tutorial).

Please cite the following paper [@qiao-prb2018] or  [@ryoo-prb2019] when
publishing SHC results obtained using this method:

> Junfeng Qiao, Jiaqi Zhou, Zhe Yuan, and Weisheng Zhao,\
> *Calculation of intrinsic spin Hall conductivity by Wannier
> interpolation*,\
> Phys. Rev. B. 98, 214402 (2018), DOI:10.1103/PhysRevB.98.214402.

or

> Ji Hoon Ryoo, Cheol-hwan Park, and Ivo Souza,\
> *Computation of intrinsic spin Hall conductivities from first
> principles using maximally localized Wannier functions*,\
> Phys. Rev. B. 99, 235113 (2019), DOI:10.1103/PhysRevB.99.235113.

### `berry_task=sc`: shift current

The shift-current contribution to the second-order response is
characterized by a frequency-dependent third-rank tensor [@sipe-prb00]
$$\label{eq:shiftcurrent}
\begin{split}
\sigma^{abc}(0;\omega,-\omega)=&-\frac{i\pi e^3}{4\hbar^2 \Omega_c N_k}
\sum_{\bm{k}} \sum_{n,m}(f_{n\bm{k}}-f_{m\bm{k}})
\times
\left(r^b_{ mn}(\bm{k})r^{c;a}_{nm}(\bm{k}) + r^c_{mn}(\bm{k})r^{b;a}_{ nm}(\bm{k})\right)\\
&\times \left[\delta(\omega_{mn\bm{k}}-\omega)+\delta(\omega_{nm\bm{k}}-\omega)\right],
\end{split}$$ where $a,b,c$ are spatial indexes and
$\omega_{mn\bm{k}}=(\epsilon_{n\bm{k}}-\epsilon_{m\bm{k}})/\hbar$. The
expression in
Eq. [\[eq:shiftcurrent\]](#eq:shiftcurrent){reference-type="ref"
reference="eq:shiftcurrent"} involves the dipole matrix element
$$\label{eq:r}
r^a_{ nm}(\bm{k})=(1-\delta_{nm})A^a_{ nm}(\bm{k}),$$ and its
*generalized derivative* $$\label{eq:gen-der}
r^{a;b}_{nm}(\bm{k})=\partial_{k_{b}} r^a_{nm}(\bm{k})
-i\left(A^b_{nn}(\bm{k})-A^b_{ mm}(\bm{k})\right)r^a_{ nm}(\bm{k}).$$
The first-principles evaluation of the above expression is technically
challenging due to the presence of an extra $k$-space derivative. The
implementation in `wannier90` follows the scheme proposed in
Ref. [@ibanez-azpiroz_ab_2018], following the spirit of the
Wannier-interpolation method for calculating the AHC [@wang-prb06] by
reformulating $k\cdot p$ perturbation theory within the subspace of
wannierized bands. This strategy inherits the practical advantages of
the sum-over-states approach, but without introducing the truncation
errors usually associated with this procedure [@sipe-prb00].

As in the case of the optical conductivity, a broadened delta function
can be applied in
Eq. [\[eq:shiftcurrent\]](#eq:shiftcurrent){reference-type="ref"
reference="eq:shiftcurrent"} by means of the parameter $\eta$ (see
Eq. [\[eq:lorentzian\]](#eq:lorentzian){reference-type="ref"
reference="eq:lorentzian"}) using the keyword
`[kubo_]smr_fixed_en_width`, and adaptive smearing can be employed using
the keyword `[kubo_]adpt_smr`.

Please cite Ref. [@ibanez-azpiroz_ab_2018] when publishing shift-current
results using this method.

### `berry_task=kdotp`: $k\cdot p$ coefficients {#sec:kdotp}

Consider a Hamiltonian $$\label{eq:H}
H=H^{0}+H^{\prime}$$ where the eigenvalues $E_{n}$ and eigenfunctions
$\vert n\rangle$ of $H^{0}$ are known, and $H^{\prime}$ is a
perturbation. In a nutshell, quasi-degenerate perturbation theory
assumes that the set of eigenfunctions of $H^0$ can be divided into
subsets $A$ and $B$ that are weakly coupled by $H^{\prime}$, and that we
are only interested in subset $A$. This theory asserts that a
transformed Hamiltonian $\tilde{H}$ exists within subspace $A$ such that
$$\label{eq:pert-exp}
\tilde{H}=\tilde{H}^{0}+\tilde{H}^{1}+\tilde{H}^{2} + \cdots$$ where
$\tilde{H}^{j}$ contain matrix elements of $H^{\prime}$ to the $j$th
power. According to Appendix B of Ref [@winkler_spin-orbit_2003], the
first three terms are $$\begin{aligned}
\label{eq:pert-matelem0}
& \tilde{H}^{0}_{mm'} = H^{0}_{mm'},\\
\label{eq:pert-matelem1}
& \tilde{H}^{1}_{mm'} = H^{'}_{mm'},\\
\label{eq:pert-matelem2}
& \tilde{H}^{2}_{mm'} = \dfrac{1}{2}\sum_{l}H^{'}_{ml}H^{'}_{lm'}
\left( 
\dfrac{1}{E_{m}-E_{l}}+\dfrac{1}{E_{m'}-E_{l}}
\right),
\end{aligned}$$ where $m,m'\in A$ and $l\in B$. The approximation
$\tilde{H}\sim \tilde{H}^{0}+\tilde{H}^{1}$ amounts to truncating $H$ to
the $A$ subspace. By further including $\tilde{H}^{2}$, the coupling to
the $B$ subspace is incorporated approximately, "renormalizing" the
elements of the truncated matrix.

We adopt the notation described in Sec. III.B of Ref. [@wang-prb06]. We
shift the origin of $k$ space to the point where the band edge (or some
other band extremum of interest) is located, and Taylor expand around
that point the Wannier-gauge Hamiltonian, $$\label{eq:HW-exp}
H^{(W)}(\bm{k})=H^{(W)}(0)
+\sum_{a}H_{a}^{(W)}(0)k_{a}
+\dfrac{1}{2}\sum_{ab}H_{ab}^{(W)}(0)k_{a}k_{b}
+ \mathcal{O}(k^{3})$$ where $a,b=x,y,z$, and $$\begin{aligned}
&H_{a}^{(W)}(0)=\left. \dfrac{\partial H^{(W)}(\bm{k})}{\partial k_{a}}\right\rvert_{\bm{k}=0}\\
&H_{ab}^{(W)}(0)=\left. \dfrac{\partial^{2} H^{(W)}(\bm{k})}{\partial k_{a}\partial k_{b}}\right\rvert_{\bm{k}=0}
\end{aligned}$$

We now apply to $H^{(W)}(\bm{k})$ a similarity transformation $U(0)$
that diagonalizes $H^{(W)}(0)$, and call the transformed Hamiltonian
$H(\bm{k})$, $$\label{eq:Hbar}
H(\bm{k})=\overbrace{\overline{H}}^{H^{0}} + \overbrace{\sum_{a}\overline{H}_{a}k_{a}
+\dfrac{1}{2}\sum_{ab}\overline{H}_{ab}k_{a}k_{b}}^{H^{\prime}} + \mathcal{O}(k^{3}),$$
where we introduced the notation
$$\overline{\mathcal{O}}=U^{\dagger}(0)\mathcal{O}^{(W)}(0)U(0),$$ and
applied it to $\mathcal{O}=H,{H}_{a},{H}_{ab}$. We can now apply
quasi-degenerate perturbation theory by choosing the diagonal matrix
$\overline{H}$ as our $H^{0}$, and the remaining (nondiagonal) terms in
Eq. [\[eq:Hbar\]](#eq:Hbar){reference-type="ref" reference="eq:Hbar"} as
$H^{\prime}$. Collecting terms in
Eq. ([\[eq:pert-exp\]](#eq:pert-exp){reference-type="ref"
reference="eq:pert-exp"}) up to second order in $k$ we get
$$\label{eq:Htilde}
 \tilde{H}_{mm'}(\bm{k}) = 
\overline{H}_{mm'} + \sum_{a} \left(\overline{H}_{a}\right)_{mm'}k_{a}
 + \dfrac{1}{2}\sum_{a,b}\left[
\left(\overline{H}_{ab}\right)_{mm'} + \left({T}_{ab}\right)_{mm'}
\right]k_{a}k_{b}+ \mathcal{O}(k^{3}),$$ where $m,m'\in A$ and we have
defined the virtual-transition matrix $$\label{eq:Tab}
\left({T}_{ab}\right)_{mm'}=\sum_{l\in B}
\left(\overline{H}_{a}\right)_{ml}\left(\overline{H}_{b}\right)_{lm'} 
 \times
\left( 
\dfrac{1}{E_{m}-E_{l}}+\dfrac{1}{E_{m'}-E_{l}}
\right)
= 
\left({T}_{ab}\right)_{m'm}^{*}.$$ (The $T_{ab}$ term in Eq.
[\[eq:Htilde\]](#eq:Htilde){reference-type="ref" reference="eq:Htilde"}
gives an Hermitean contribution to $\tilde{H}(\bm{k})$ only after
summing over $a$ and $b$, whereas the other terms are Hermitean already
before summing.)

The implementation in `wannier90` follows the scheme proposed in
Ref. [@ibanez-azpiroz-ArXiv2019], and outputs $\overline{H}_{mm'}$ in
`seedname-kdotp_0.dat`, $\left(\overline{H}_{a}\right)_{mm'}$ in
`seedname-kdotp_1.dat`, and
$\left[\left(\overline{H}_{ab}\right)_{mm'} + \left({T}_{ab}\right)_{mm'}\right]/2$
in `seedname-kdotp_2.dat`.

Please cite Ref. [@ibanez-azpiroz-ArXiv2019] when publishing $k\cdot p$
results using this method.

### Needed matrix elements

All the quantities entering the formulas for the optical conductivity
and AHC can be calculated by Wannier interpolation once the Hamiltonian
and position matrix elements $\langle {\bf 0}n\vert H\vert
{\bf R}m\rangle$ and $\langle {\bf 0}n\vert {\bf r}\vert {\bf
  R}m\rangle$ are known [@wang-prb06; @yates-prb07]. Those matrix
elements are readily available at the end of a standard MLWF calculation
with `wannier90`. In particular, $\langle {\bf
  0}n\vert {\bf r}\vert {\bf R}m\rangle$ can be calculated by Fourier
transforming the overlap matrices in Eq. (1.7),
$$\langle u_{n{\bf k}}\vert u_{m{\bf k}+{\bf b}}\rangle.$$ Further
Wannier matrix elements are needed for the orbital
magnetization [@lopez-prb12]. In order to calculate them using Fourier
transforms, one more piece of information must be taken from the
$k$-space *ab-initio* calculation, namely, the matrices
$$\langle u_{n{\bf k}+{\bf b}_1}\vert
H_{\bf k}\vert u_{m{\bf k}+{\bf b}_2}\rangle$$ over the *ab-initio*
$k$-point mesh [@lopez-prb12]. These are evaluated by `pw2wannier90`,
the interface routine between ` pwscf` and `wannier90`, by adding to the
input file ` seedname.pw2wan` the line $${\tt
%\begin{quote}
write\_uHu = .true.
%\end{quote}
}$$ The calculation of spin Hall conductivity needs the spin matrix
elements
$$\langle u_{n{\bf k}}\vert \sigma_\gamma \vert u_{m{\bf k}}\rangle, 
\gamma = x, y, z$$ from the *ab-initio* $k$-point mesh. These are also
evaluated by `pw2wannier90` by adding to the input file
` seedname.pw2wan` the line $${\tt
	%\begin{quote}
	write\_spn = .true.
	%\end{quote}
}$$ If one uses Ryoo's method to calculate spin Hall conductivity, the
further matrix elements are needed: $$\langle u_{n{\bf k}}\vert
\sigma_\gamma H_{\bf k}\vert u_{m{\bf k}+{\bf b}}\rangle, \langle u_{n{\bf k}}\vert
\sigma_\gamma \vert u_{m{\bf k}+{\bf b}}\rangle,
\gamma = x, y, z$$ and these are evaluated by adding to the input file
` seedname.pw2wan` the lines $${\tt
	write\_sHu = .true.
}$$ $${\tt
	write\_sIu = .true.
}$$
