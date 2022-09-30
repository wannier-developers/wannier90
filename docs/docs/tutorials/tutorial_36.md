# 36: Silicon &#151; Valence and low-lying conduction states

## Marzari-Vanderbilt functional

- Outline: *Obtain MLWFs for the valence and low-lying conduction-band
    states of Si by minimizing the Marzari-Vanderbilt functional.
    Plot the interpolated bandstructure.*

- Directory: [`tutorial/tutorial36/`](https://github.com/wannier-developers/wannier90/tree/develop/tutorials/tutorial36)

- Input Files

    - `silicon.scf` *The `pwscf` input file for ground
        state calculation*

    - `silicon.nscf` *The `pwscf` input file to obtain
        Bloch states on a uniform grid*

    - `silicon.pw2wan` *Input file for `pw2wannier90`*

    - `silicon.win` *The `wannier90` input file*

1. Run `pwscf` to obtain the ground state of silicon

    ```bash title="Terminal"
    pw.x < silicon.scf > scf.out
    ```

2. Run `pwscf` to obtain the Bloch states on a uniform
    k-point grid. Note that we request the lower 4 (valence) bands

    ```bash title="Terminal"
    pw.x < silicon.nscf > nscf.out
    ```

3. Run `wannier90` to generate a list of the required overlaps (written
    into the `silicon.nnkp` file).

    ```bash title="Terminal"
    wannier90.x -pp silicon
    ```

4. Run `pw2wannier90` to compute the overlap between Bloch states and
    the projections for the starting guess (written in the `silicon.mmn`
    and `silicon.amn` files).

    ```bash title="Terminal"
    pw2wannier90.x < silicon.pw2wan > pw2wan.out
    ```

5. Run `wannier90` to compute the MLWFs.

    ```bash title="Terminal"
    wannier90.x silicon
    ```

6. Plot the bandstructure by adding the following commands to the input
    file `silicon.win`

    ```vi title="Input file"
    restart = plot
    
    bands_plot = true
    ```

    and re-running `wannier90`. The files `silicon_band.dat` and
    `silicon_band.gnu` are created. To plot the bandstructure using
    gnuplot

    ```bash title="Terminal"
    gnuplot
    ```

    ```gnuplot title="Gnuplot shell"
    load 'silicon_band.gnu'
    ```

    The k-point path for the bandstructure interpolation is set in the
    `kpoint_path` block. Try plotting along different paths.

## Stengel-Spaldin functional

- Outline: *Obtain MLWFs for the valence and low-lying conduction-band
    states of Si by minimizing the Stengel-Spaldin functional.
    Plot the interpolated bandstructure.*

- Directory: [`tutorial/tutorial36/`](https://github.com/wannier-developers/wannier90/tree/develop/tutorials/tutorial36)

- Input Files

    - `silicon.scf` *The `pwscf` input file for ground
        state calculation*

    - `silicon.nscf` *The `pwscf` input file to obtain
        Bloch states on a uniform grid*

    - `silicon.pw2wan` *Input file for `pw2wannier90`*

    - `silicon.win` *The `wannier90` input file*

- To run, the procedure is the same as the MV functional case.

## Further ideas

- Compare the spreads obtained by both functionals.
  Increase the grid and find out that the two functionals converges to the same behavior.

- Compare the Wannier-interpolated bandstructure with the both functionals.
  Find out that the interpolated bandstructures are almost the same.
