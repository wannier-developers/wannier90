# Some notes on the interpolation

In `wannier90` v.2.1, a new flag `use_ws_distance` has been introduced
(and it is set to `.true.` by default since version v3.0). Setting it to
`.false.` reproduces the "standard" behavior of `wannier90` in v.2.0.1
and earlier, while setting it to `.true.` changes the interpolation
method as described below. In general, this allows a smoother
interpolation, helps reducing (a bit) the number of $k-$points required
for interpolation, and reproduces the band structure of large supercells
sampled at $\Gamma$ only (setting it to `.false.` produces instead flat
bands, which might instead be the intended behaviour for small molecules
carefully placed at the centre of the cell).

The core idea rests on the fact that the Wannier functions
$w_{n\bm{\mathrm{R}}}(\bm{\mathrm{r}})$ that we build from
$N\times M\times L$ $k-$points are actually periodic over a supercell of
size $N\times M\times L$, but when you use them to interpolate you want
them to be *zero* outside this supercell. In 1D it is pretty obvious
want we mean here, but in 3D what you really want that they are zero
outside the Wigner--Seitz cell of the $N\times M\times L$ superlattice.

The best way to impose this condition is to check that every real-space
distance that enters in the $R\to k$ Fourier transform is the shortest
possible among all the $N\times M\times L-$periodic equivalent copies.

If the distances were between unit cells, this would be trivial, but the
distances are between Wannier functions which are not centred on
$\bm{\mathrm{R}}=0$. Hence, when you want to consider the matrix element
of a generic operator $\bm{\mathrm{O}}$ (i.e., the Hamiltonian)
$\langle w_{i\bm{\mathrm{0}}}(\bm{\mathrm{r}})|\bm{\mathrm{O}}|w_{j\bm{\mathrm{R}}}(\bm{\mathrm{r}})\rangle$
you must take in account that the centre $\bm{\mathrm{\tau}}_i$ of
$w_{i\bm{\mathrm{0}}}(\bm{\mathrm{r}})$ may be very far away from
$\bm{\mathrm{0}}$ and the centre $\bm{\mathrm{\tau}}_j$ of
$w_{j\bm{\mathrm{R}}}(\bm{\mathrm{r}})$ may be very far away from
$\bm{\mathrm{R}}$.

There are many way to find the shortest possible distance between
$w_{i\bm{\mathrm{0}}}(\bm{\mathrm{r}})$ and
$w_{j\bm{\mathrm{R}}}(\bm{\mathrm{r}}-\bm{\mathrm{R}})$, the one used
here is to consider the distance
$\bm{\mathrm{d}}_{ij\bm{\mathrm{R}}} = \bm{\mathrm{\tau}}_i - (\bm{\mathrm{\tau}}_j+\bm{\mathrm{R}})$
and all its superlattice periodic equivalents
$\bm{\mathrm{d}}_{ij\bm{\mathrm{R}}}+ \bm{\mathrm{\tilde R}}_{nml}$,
with
$\bm{\mathrm{\tilde R}}_{nml} = (Nn\bm{\mathrm{a}}_1 + Mm\bm{\mathrm{a}}_2 + Ll\bm{\mathrm{a}}_3)$
and $n,l,m = {-L,-L+1,...0,...,L-1,L}$, with $L$ controlled by the
parameter `ws_search_size`.

Then,

1.  if
    $\bm{\mathrm{d}}_{ij\bm{\mathrm{R}}}+ \bm{\mathrm{\tilde R}}_{nml}$
    is inside the $N\times M \times L$ super-WS cell, then it is the
    shortest, take it and quit

2.  if it is outside the WS, then it is not the shortest, throw it away

3.  if it is on the border/corner of the WS then it is the shortest, but
    there are other choices of $(n,m,l)$ which are equivalent, find all
    of them

In all distance comparisons, a small but finite tolerance is considered,
which can be controlled with the parameter `ws_distance_tol`.

Because of how the Fourier transform is defined in the `wannier90` code
(not the only possible choice) it is only
$\bm{\mathrm{R}}+\bm{\mathrm{\tilde R}}_{nml}$ that enters the
exponential, but you still have to consider the distance among the
actual centres of the Wannier functions. Using the centres of the
unit-cell to which the Wannier functions belong is not enough (but is
easier, and saves you one index).

Point 3 is not stricly necessary, but using it helps enforcing the
symmetry of the system in the resulting band structure. You will get
some small but evident symmetry breaking in the band plots if you just
pick one of the equivalent $\bm{\mathrm{\tilde R}}$ vectors.

Note that in some cases, all this procedure does absolutely nothing, for
instance if all the Wannier function centres are very close to 0 (e.g.,
a molecule carefully placed in the periodic cell).

In some other cases, the effect may exist but be imperceptible. E.g., if
you use a very fine grid of $k-$points, even if you don't centre each
functions perfectly, the periodic copies will still be so far away that
the change in centre applied with $\tt use\_ws\_distance$ does not
matter.

When instead you use few $k-$points, activating the
$\tt use\_ws\_distance$ may help a lot in avoiding spurious oscillations
of the band structure even when the Wannier functions are well
converged.
