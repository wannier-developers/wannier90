# 23: Silicon &#151; $G_0W_0$ bands structure interpolation {#silicon-g_0w_0-bands-structure-interpolation .unnumbered}

Note: This tutorial requires a recent version of the `ypp`
post-processing code of `yambo`.

-   Outline: *Interpolate the bands structure of silicon obtained from
    many-body perturbation theory at the $G_0W_0$ level. Using the
    `yambo` code, the quasi-particle corrections (QP) are summed to
    Kohn-Sham eigenvalues, while the wavefunctions remain the same. *

-   Directory: `tutorials/tutorial23/` *Files can be downloaded from [here](https://github.com/wannier-developers/wannier90/tutorials/tutorial23)*

-   Input Files

    -    `silicon.scf` *The `pwscf` input file for the
        ground state calculation*

    -    `silicon.nscf ` *The `pwscf` input file to obtain
        Bloch states on a uniform grid*

    -    `silicon.gw.nscf ` *The `pwscf` input file to
        obtain Bloch states on a reduced grid with many empty bands*

    -    `silicon.pw2wan` *The input file for `pw2wannier90`*

    -    `silicon.win` *The `wannier90` input file*

    -    `silicon.gw.win` *The `wannier90` input file* (for the $G_0W_0$
        step)

    -    `yambo.in` *The `yambo` input file*

    -    `ypp.in` *The `ypp` input file*

1.  Copy the input files from the `INPUT directory` into a working
    directory (e.g. `WORK`)

2.  Run `pwscf` to obtain the ground state charge of
    silicon

    ```bash title="Terminal"
    pw.x < silicon.scf > scf.out
    ```

3.  Run `pwscf` to obtain the Bloch states reduced grid. We
    use a 8x8x8 with many bands (many empty bands are needed to perform
    a $G_0W_0$ with `yambo`)

    ```bash title="Terminal"
    pw.x < silicon.gw.nscf > nscf.gw.out
    ```

4.  Use the `k_mapper.py` utility to find the indexes of a 4x4x4 uniform
    grid into the 8x8x8 reduced grid

    ```bash title="Terminal"
    ./k_mapper.py 4 4 4 "../tutorials/tutorial23/WORK/nscf.gw.out"
    ```

    Use the output to complete the `yambo.in` input file (you also need
    to specify how many bands you want to compute the QP
    corrections, here you can use all the bands from 1 to 14). Then, you
    should have obtained something like:

    ```vi title="Output file"
     1| 1| 1|14|
    3| 3| 1|14|
    5| 5| 1|14|
    13| 13| 1|14|
    ...
    ```

5.  Enter the `si.save` directory and run `p2y`. A `SAVE` folder is
    created, you can move it up in the `/WORK/` directory.

6.  Run a $G_0W_0$ calculation from the `/WORK/` directory (remember, we
    are using a 8x8x8 grid but computing QP corrections only on a 4x4x4
    grid)

    ```bash title="Terminal"
    yambo 
    ```

7.  Run `pwscf` to obtain the Bloch states on a uniform
    k-point grid

    ```bash title="Terminal"
    pw.x < silicon.nscf > nscf.out
    ```

8.  Run `wannier90` to generate a list of the required overlaps (written
    into the `silicon.nnkp` file).

    ```bash title="Terminal"
    wannier90.x -pp silicon
    ```

9.  Run `pw2wannier90` to compute the overlap between Bloch states, the
    projections for the starting guess (written in the `silicon.mmn` and
    `silicon.amn` respectively).

    ```bash title="Terminal"
    pw2wannier90.x < silicon.pw2wan > pw2wan.out
    ```

10. Run `wannier90` to compute the MLWFs.

    ```bash title="Terminal"
    wannier90.x silicon
    ```

    At this point, you should have obtained the interpolated valence
    bands for silicon at the DFT level.

11. Run a `ypp` calculation (just type `ypp`)

    ```bash title="Terminal"
    ypp
    ```

    You should obtain a file `silicon.gw.unsorted.eig` which contains
    the QP corrections on a uniform 4x4x4 grid.

12. Run the gw2wannier90.py script to reorder, align and correct all
    matrices and files using the QP corrections

    ```bash title="Terminal"
    ../../../utility/gw2wannier90.py silicon mmn amn
    ```

13. Run `wannier90` to compute the MLWFs.

    ```bash title="Terminal"
    wannier90.x silicon.gw
    ```

    At this point, you should have obtained the interpolated valence
    bands for silicon at the $G_0W_0$ level.

After you completed the tutorial for the valence bands only, you can
repeat the final steps to interpolate also some conduction bands using
disentanglement (the code is already present as comments in the input
files).


