# Overview of the `gyrotropic` module

The `gyrotropic` module of `postw90` is called by setting
` gyrotropic = true` and choosing one or more of the available options
for `gyrotropic_task`. The module computes the quantities, studied in
[@tsirkin-arxiv17], where more details may be found.

### `gyrotropic_task=-d0`: the Berry curvature dipole 

The traceless dimensionless tensor 

$$
\begin{equation}
\label{eq:D_ab}
D_{ab}=\int[d{\bm k}]\sum_n
\frac{\partial E_n}{\partial{k_a}}
\Omega_n^b
\left(-\frac{\partial f_0}{\partial E}\right)_{E=E_n},
\end{equation}
$$

### `gyrotropic_task=-dw`: the finite-frequency generalization of the Berry curvature dipole 

$$
\begin{equation}
\label{eq:D-tilde}
\widetilde{D}_{ab}(\omega)=\int[d{\bm k}]\sum_n
\frac{\partial E_n}{\partial{k_a}}\widetilde\Omega^b_n(\omega)
\left(-\frac{\partial f_0}{\partial E}\right)_{E=E_n},
\end{equation}
$$

where $\widetilde{\bm\Omega}_{{\bm k}n}(\omega)$ is a finite-frequency
generalization of the Berry curvature: 

$$
\begin{equation}
\label{eq:curv-w}
\widetilde{\bm\Omega}_{{\bm k}n}(\omega)=-
\sum_m\,\frac{\omega^2_{{\bm k}mn}}{\omega^2_{{\bm k}mn}-\omega^2}
{\rm Im}\left({\bm A}_{{\bm k}nm}\times{\bm A}_{{\bm k}mn}\right)
\end{equation}
$$

Contrary to the Berry curvature, the divergence of
$\tilde{\bm\Omega}_{{\bm k}n}(\omega)$ is generally nonzero. As a
result, $\widetilde{D}(\omega)$ can have a nonzero trace at finite
frequencies, $\tilde{D}_\|\neq-2\tilde{D}_\perp$ in Te.

### `gyrotropic_task=-C`: the ohmic conductivity 

In the constant relaxation-time approximation the ohmic conductivity is
expressed as $\sigma_{ab}=(2\pi e\tau/\hbar)C_{ab}$, with

$$
\begin{equation}
\label{eq:C_ab}
C_{ab}=\frac{e}{h}\int[d{\bm k}]\sum_n\,
\frac{\partial E_n}{\partial{k_a}} \frac{\partial E_n}{\partial{k_b}}
\left(-\frac{\partial f_0}{\partial E}\right)_{E=E_n}
\end{equation}
$$

a positive
quantity with units of surface current density (A/cm).

### `gyrotropic_task=-K`: the kinetic magnetoelectric effect (kME) 

A microscopic theory of the intrinsic kME effect in bulk crystals was
recently developed [@yoda-sr15; @zhong-prl16].

The response is described by 

$$
\begin{equation}
\label{eq:K_ab}
K_{ab}=\int[d{\bm k}]\sum_n\frac{\partial E_n}{\partial{k_a}} m_n^b 
\left(-\frac{\partial f_0}{\partial E}\right)_{E=E_n},
\end{equation}
$$ 

which has the
same form as Eq. $\eqref{eq:D_ab}$, but with the Berry curvature replaced by the
intrinsic magnetic moment ${\bm m}_{{\bm k}n}$ of the Bloch electrons,
which has the spin and orbital components given by [@xiao-rmp10]

$$
\begin{equation}
\label{eq:m-spin}
\begin{aligned}
m^{\rm spin}_{{\bm k}n}&=-\frac{1}{2}g_s\mu_{\rm B} \langle\psi_{{\bm k}
      n}\vert\bf \sigma\vert\psi_{{\bm k}n}\rangle\\
% \label{eq:m-orb}
{\bm m}^{\rm orb}_{{\bm k}n}&=\frac{e}{2\hbar}{\rm Im}
\langle{\bm\partial}_{\bm k}u_{{\bm k}n}\vert\times
(H_{\bm k}-E_{{\bm k}n})\vert{\bm\partial}_{\bm k}u_{{\bm k}n}\rangle,
\end{aligned}
\end{equation}
$$ 

where $g_s\approx 2$ and we chose $e>0$.

### `gyrotropic_task=-dos`: the density of states 

The density of states is calculated with the same width and type of
smearing, as the other properties of the `gyrotropic` module

### `gyrotropic_task=-noa`: the interband contributionto the natural optical activity 

Natural optical rotatory power is given by [@ivchenko-spss75]

$$
\begin{equation}
\label{eq:rho-c}
\rho_0(\omega)=\frac{\omega^2}{2c^2}{\rm Re}\,\gamma_{xyz}(\omega).
\end{equation}
$$

for light propagating ling the main symmetry axis of a crystal $z$. Here
$\gamma_{xyz}(\omega)$ is an anti-symmetric (in $xy$) tensor with units
of length, which has both inter- and intraband contributions.

Following Ref. [@malashevich-prb10] for the interband contribution we
writewe write, with $\partial_c\equiv\partial/\partial k_c$,

$$
\begin{equation}
\label{eq:gamma-inter}
\begin{gathered}
{\rm Re}\,\gamma_{abc}^{\mathrm{inter}}(\omega)=\frac{e^2}{\varepsilon_0\hbar^2}
\int[d{\bm k}]
\sum_{n,l}^{o,e}\,
\Bigl[ \frac{1}{\omega_{ln}^2-\omega^2} 
{\rm Re}\left(A_{ln}^bB_{nl}^{ac}-A_{ln}^aB_{nl}^{bc}\right) \\
-\frac{3\omega_{ln}^2-\omega^2}{(\omega_{ln}^2-\omega^2)^2} 
\partial_c(E_l+E_n){\rm Im}\left(A_{nl}^aA_{ln}^b\right)   
\Bigr].
\end{gathered}
\end{equation}
$$ 

The summations over $n$ and $l$ span the occupied ($o$)
and empty ($e$) bands respectively, $\omega_{ln}=(E_l-E_n)/\hbar$, and
${\bm A}_{ln}({\bm k})$ is given by
[Berry Eq. (3)](../berry/#mjx-eqn:eq:berry-connection-matrix)
<!-- Eq. $\eqref{eq:berry-connection-matrix}$.  -->
Finally, the matrix
$B_{nl}^{ac}$ has both orbital and spin contributions given by

$$
\begin{equation}
\label{eq:B-ac-orb}
B_{nl}^{ac\,({\rm orb})}=
  \langle u_n\vert(\partial_aH)\vert\partial_c u_l\rangle
 -\langle\partial_c u_n\vert(\partial_aH)\vert u_l\rangle
\end{equation}
$$

and

$$
\begin{equation}
\label{eq:B-ac-spin}
B_{nl}^{ac\,({\rm spin})}=-\frac{i\hbar^2}{m_e}\epsilon_{abc}
\langle u_n\vert\sigma_b\vert u_l\rangle.
\end{equation}
$$ 

The spin matrix elements
contribute less than 0.5% of the total $\rho_0^{\rm inter}$ of Te.
Expanding $H=\sum_m \vert u_m\rangle E_m \langle u_m\vert$ we obtain for
the orbital matrix elements

$$
\begin{equation}
\label{eq:Bnl-sum}
B_{nl}^{ac\,({\rm orb})}=-i\partial_a(E_n+E_l)A_{nl}^c \sum_m \Bigl\{ (E_n-E_m) A_{nm}^aA_{ml}^c -(E_l-E_m) A_{nm}^cA_{ml}^a \Bigr\}.
\end{equation}
$$ 

This reduces the calculation of $B^{\text{(orb)}}$
to the evaluation of band gradients and off-diagonal elements of the
Berry connection matrix. Both operations can be carried out efficiently
in a Wannier-function basis following Ref. [@yates-prb07].

### `gyrotropic_task=-spin`: compute also the spin component of NOA and KME 

Unless this task is specified, only the orbital contributions are
calcuated in NOA and KME, thus contributions from
Eqs. $\eqref{eq:m-spin}$
and $\eqref{eq:B-ac-spin}$ are omitted.
