# 35: Silicon — Projectability-disentangled Wannier functions with custom projectors

## Outline

Obtain MLWFs for silicon using projectability disentanglement, with additional
$3d$ projectors to describe high-energy conduction bands.
For more details on the methodology, see Ref.[@Qiao2023-pdwf].

## Input files
<!-- markdownlint-disable code-block-style -->

- Directory: [`tutorial/tutorial35/`](https://github.com/wannier-developers/wannier90/tree/develop/tutorials/tutorial35)

- `silicon.scf` The `pw.x` input file for ground state calculation

    ??? quote "silicon.scf"

        ```fortran linenums="1" title="Input file"
        --8<-- "tutorials/tutorial35/silicon.scf"
        ```

- `silicon.bands` The `pw.x` input file for band structure calculation

    ??? quote "silicon.bands"

        ```fortran linenums="1" title="Input file"
        --8<-- "tutorials/tutorial35/silicon.bands::45"
        ...
        ...
        ```

- `silicon.bandsx` The `bands.x` input file for extracting band structure eigenvalues

    ??? quote "silicon.bandsx"

        ```fortran linenums="1" title="Input file"
        --8<-- "tutorials/tutorial35/silicon.bandsx"
        ```

- `silicon.nscf` The `pw.x` input file to obtain Bloch states on a uniform grid

    ??? quote "silicon.nscf"

        ```fortran linenums="1" title="Input file"
        --8<-- "tutorials/tutorial35/silicon.nscf::45"
        ...
        ...
        ```

- `silicon.pw2wan` Input file for `pw2wannier90.x`

    ??? quote "silicon.pw2wan"

        ```fortran linenums="1" title="Input file"
        --8<-- "tutorials/tutorial35/silicon.pw2wan"
        ```

- `silicon.win` The `wannier90.x` input file

    ??? quote "silicon.win"

        ```fortran linenums="1" title="Input file"
        --8<-- "tutorials/tutorial35/silicon.win::50"
        ...
        ...
        ```

- `silicon_bandsdiff.gnu` The gnuplot script to compare DFT and Wannier bands

    ??? quote "silicon_bandsdiff.gnu"

        ```gnuplot linenums="1" title="Gnuplot script"
        --8<-- "tutorials/tutorial35/silicon_bandsdiff.gnu"
        ```

## Steps

1. Run `pw.x` to obtain the ground state of silicon

    ```bash title="Terminal"
    pw.x < silicon.scf > scf.out
    ```

2. Run `pw.x` to obtain the band structure of silicon

    ```bash title="Terminal"
    pw.x < silicon.bands > bands.out
    ```

3. Run `bands.x` to obtain a `silicon.bands.dat` file containing the band
    structure of silicon

    ```bash title="Terminal"
    bands.x < silicon.bandsx > bandsx.out
    ```

4. Run `pw.x` to obtain the Bloch states on a uniform k-point grid

    ```bash title="Terminal"
    pw.x < silicon.nscf > nscf.out
    ```

5. Run `pw.x` to generate a list of the required overlaps (written into the
    `silicon.nnkp` file).

    !!! note

        See `win` input file, no need to specify initial projections,
        they are chosen from the pseudo-atomic orbitals inside the
        `ext_proj/Si.dat` file.

    ```bash title="Terminal"
    wannier90.x -pp silicon
    ```

6. Run `pw2wannier90.x` to compute the overlap between Bloch states and
    the projections for the starting guess (written in the `silicon.mmn`
    and `silicon.amn` files).

    ```bash title="Terminal"
    pw2wannier90.x < silicon.pw2wan > pw2wan.out
    ```

7. Run `pw.x` to compute the MLWFs.

    ```bash title="Terminal"
    wannier90.x silicon
    ```

8. Run `gnuplot` to compare DFT and Wannier-interpolated bands, this
    will generate a PDF file `silicon_bandsdiff.pdf`, see
    Fig.[Bands comparison](#fig:silicon_bandsdiff).

    ```bash title="Terminal"
    ./silicon_bandsdiff.gnu
    ```

    <figure markdown="span" id="fig:silicon_bandsdiff">
    ![Bands diff](./silicon_bandsdiff_spd.webp){width="500"}
    <figcaption markdown="span">Comparison of DFT and Wannier bands for silicon.
    </figcaption>
    </figure>

9. (Optional) Clean up all output files

    ```bash title="Terminal"
    make clean
    ```

## Further ideas

1. Try changing the `atom_proj_exclude` in `silicon.pw2wan` file, i.e.,
    these commented lines

    ```fortran linenums="10" title="Input file" hl_lines="5"
    --8<-- "tutorials/tutorial35/silicon.pw2wan:10:14"
    ```

    !!! hint

        You need to set `num_wann = 8` in the `silicon.win` file as well, to
        be consistent with the reduced number of projectors.

2. Now that $3d$ projectors provide us a larger space for optimization,
    you can try increasing the `dis_froz_max` to freeze higher energy
    bands, if you are targeting at reproducing those eigenvalues.

    !!! note

        The `dis_proj_min/max` and `dis_froz_min/max` can be
        enabled simultaneously: the union of inner energy window and
        high-projectability states will be freezed, and the union of states
        outside outer energy window and having low projectability will be
        discarded. Thus, you can still use energy window to make sure
        near-Fermi energy states are well reproduced, and use
        "projectability window" to selectively freeze atomic-like states in
        the conduction region.

3. The default `dis_proj_max = 0.95` might not freeze all the states
    you want, try changing this value and see the band interpolation
    results. For other materials, it might worth trying decreasing this
    value to freeze more states.
