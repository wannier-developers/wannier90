# 24: Tellurium &#151; gyrotropic effects {#tellurium-gyrotropic-effects .unnumbered}

-   Outline: *Calculate the gyrotropic effects in trigonal right-handed
    Te* Similar to the calculations of [@tsirkin-arxiv17]

-   Directory: `tutorials/tutorial24/` *Files can be downloaded from [here](https://github.com/wannier-developers/wannier90/tree/develop/tutorials/tutorial24)*

-   Input files

    -   `Te.scf` *The `pwscf` input file for ground state
        calculation*

    -   `Te.nscf` *The `pwscf` input file to obtain Bloch
        states on a uniform grid*

    -   `Te.pw2wan` *The input file for `pw2wannier90`*

    -   `Te.win` *The `wannier90` input file*

To make things easy, the tutorial treats Te without spin-orbit

1.  Run `pwscf` to obtain the ground state of tellurium

    ```bash title="Terminal"
    pw.x < Te.scf > scf.out
    ```

2.  Run `pwscf` to obtain the Bloch states on a uniform
    `3x3x4` k-point grid

    ```bash title="Terminal"
    pw.x < Te.nscf > nscf.out
    ```

3.  Run `wannier90` to generate a list of the required overlaps (written
    into the `Te.nnkp` file).

    ```bash title="Terminal"
    wannier90.x -pp Te
    ```

4.  Run `pw2wannier90` to compute:

    -   The overlaps $\langle u_{n{\bf k}}\vert u_{m{\bf k}+{\bf
                  b}}\rangle$ (written in the `Te.mmn` file)

    -   The projections for the starting guess (written in the ` Te.amn`
        file)

    -   The matrix elements $\langle u_{n{\bf k}+{\bf b}_1}\vert
              H_{\bf k}\vert u_{m{\bf k}+{\bf b}_2}\rangle$ (written in
        the `Te.uHu` file)

    -   The spin matrix elements $\langle \psi_{n{\bf
                k}}\vert \sigma_i\vert \psi_{m{\bf k}}\rangle$ (would be
        written in the `Te.spn` file, but only if spin-orbit is
        included, which is not the case for the present tutorial)

    ```bash title="Terminal"
    pw2wannier90.x < Te.pw2wan > pw2wan.out
    ```

5.  Run `wannier90` to compute the MLWFs.

    ```bash title="Terminal"
    wannier90.x Te
    ```

6.  Add the following lines to the `wannier90.win` file:

    ```vi title="Input file"
    gyrotropic=true
    gyrotropic_task=-C-dos-D0-Dw-K
    fermi_energy_step=0.0025
    fermi_energy_min=5.8
    fermi_energy_max=6.2
    gyrotropic_freq_step=0.0025
    gyrotropic_freq_min=0.0
    gyrotropic_freq_max=0.1
    gyrotropic_smr_fixed_en_width=0.01
    gyrotropic_smr_max_arg=5
    gyrotropic_degen_thresh=0.001
    gyrotropic_box_b1=0.2 0.0 0.0
    gyrotropic_box_b2=0.0 0.2 0.0
    gyrotropic_box_b3=0.0 0.0 0.2
    gyrotropic_box_center=0.33333 0.33333 0.5
    gyrotropic_kmesh=50 50 50 
    ```

7.  Run `postw90` 
    to compute the gyrotropic properties: tensors $D$, $\widetilde{D}$,
    $K$, $C$ (See the User Guide):


    ```bash title="Terminal"
    postw90.x Te # (1)!
    mpirun -np 8 postw90.x Te  # (2)!
    ```

    1.   serial execution
    2.   example of parallel execution with 8 MPI processes


    The integration in the $k$-space is limited to a small area around
    the H point. Thus it is valid only for Fermi levels near the band
    gap. And one needs to multiply the results by 2, to account for the
    H' point. To integrate over the entire Brillouin zone, one needs to
    remove the `gyrotropic_box_`$\ldots$ parameters

8.  Now change the above lines to

    ```vi title="Input file"
    gyrotropic=true
    gyrotropic_task=-NOA
    fermi_energy=5.95
    gyrotropic_freq_step=0.0025
    gyrotropic_freq_min=0.0
    gyrotropic_freq_max=0.3
    gyrotropic_smr_fixed_en_width=0.01
    gyrotropic_smr_max_arg=5
    gyrotropic_band_list=4-9
    gyrotropic_kmesh=50 50 50
    ```

    and compute the interband natural optical activity

    ```bash title="Terminal"
    postw90.x Te # (1)!
    mpirun -np 8 postw90.x Te  # (2)!
    ```

    1.   serial execution
    2.   example of parallel execution with 8 MPI processes
