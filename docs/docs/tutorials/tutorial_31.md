# 31: Platinum &#151; Selected columns of density matrix algorithm for spinor wavefunctions

Note: This tutorial requires a recent version of the `pw2wannier90.x`
post-processing code of `Quantum ESPRESSO` (v6.3 or above).

- Outline: *For bulk crystalline platinum with spin-orbit coupling,
    generate the $A_{mn}$ matrices via the selected columns of density
    matrix (SCDM) algorithm and the corresponding spinor-MLWFs. To
    better understand the input files and the results of these
    calculations, it is crucial that the Reader has familiarized with
    the concepts and methods explained in Ref. [@LinLin-ArXiv2017]. More
    info on the keywords related to the SCDM method may be found in the
    user_guide.*

    This tutorial focuses on the use of the SCDM method for
    spin-noncollinear systems. For the overview of the use of SCDM
    method to spinless systems, please refer to this [tutorial](tutorial_27.md).

- Directory: `tutorials/tutorial31/` *Files can be downloaded from
    [here](https://github.com/wannier-developers/wannier90/tree/develop/tutorials/tutorial31)*

    The input files for this tutorials are similar to the ones in Tutorial
    [29](tutorial_29.md), except that a coarser k-point grid is used and that
    the keywords related to `postw90.x` are removed.

- Input Files:

    - `Pt.scf` *The `pwscf` input file for the ground
        state calculation*

    - `Pt.nscf` *The `pwscf` input file to obtain Bloch
        states on a uniform grid*

    - `Pt.pw2wan` *The input file for `pw2wannier90` with keywords
        related to the SCDM method*

    - `Pt.win` *The `wannier90` input file*

We will compute 18 localized WFs. Since the band structure of platinum
is metallic, the low-lying bands are entangled with other high-energy
bands, and the columns of the density matrix are not exponentially
localized by construction. Thus, we use a modified density matrix
[@LinLin-ArXiv2017], with the function $f(\varepsilon_{n,\mathbf{k}})$
defined as a complementary error function. Refer to Tutorial [27](tutorial_27.md)
for the definition of the modified density matrix and the functional form of
$f(\varepsilon_{n,\mathbf{k}})$.

1. Run `pwscf` to obtain the ground state of platinum

    ```bash title="Terminal"
    pw.x < Pt.scf > scf.out
    ```

2. Run `pwscf` to obtain the Bloch states on a uniform
    $7\times 7\times 7$ $k$-point grid

    ```bash title="Terminal"
    pw.x < Pt.nscf > nscf.out
    ```

3. Inspect the `Pt.win` input file and make sure that the
    `auto_projections` flag is set to `.true.`. Also, make sure that no
    projection block is present.

4. Run `wannier90` to generate a list of the required overlaps (written
    into the `Pt.nnkp` file)

    ```bash title="Terminal"
    wannier90.x -pp Pt
    ```

5. Inspect the `Pt.nnkp` file and make sure you find the
    `auto_projections` block and that no projections have been written
    in the `projections` block.

6. Inspect the `Pt.pw2wan` input file. You will find four SCDM-related
    keywords: `scdm_proj`, `scdm_entanglement`, `scdm_mu` and
    `scdm_sigma`. In particular, the keyword `scdm_proj` will instruct
    `pw2wannier90.x` to use the SCDM method when generating the $A_{mn}$
    matrix. The remaining three keywords defines the formula and
    parameters to define the function $f(\varepsilon_{n\mathbf{k}})$
    (see Ref. [@LinLin-ArXiv2017] and Tutorial [27](tutorial_27.md)).

7. Run `pw2wannier90` to compute the overlap between Bloch states and
    the projections via the SCDM method (written in the `Pt.mmn` and
    `Pt.amn` respectively).

    ```bash title="Terminal"
    pw2wannier90.x < Pt.pw2wan > pw2wan.out
    ```

8. Inspect the `pw2wan.out` output file. Compared to the spinless case,
    you will find the following two additional lines.

    ```vi title="Output file"
             Number of pivot points with spin up  :     9
             Number of pivot points with spin down:     9
    ```

    These lines give information on the pivots obtained by the QR
    decomposition with column pivoting (QRCP) in the SCDM algorithm.
    Each pivot determines a point in the real-space grid and a spin
    state. The basis of the spin state is determined by the basis used
    in the electronic structure code. In `pwscf`, the basis
    states are spin up and down states along the Cartesian $z$-axis.

9. Run `wannier90` to compute the MLWFs

    ```bash title="Terminal"
    wannier90.x Pt
    ```
