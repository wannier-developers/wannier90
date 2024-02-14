# 25: Gallium Arsenide &#151; Nonlinear shift current {#gallium-arsenide-nonlinear-shift-current .unnumbered}

-   Outline: *Calculate the nonlinear shift current of inversion
    asymmetric fcc Gallium Arsenide. In preparation for this example it
    may be useful to read Ref. [@ibanez-azpiroz_ab_2018]*

-   Directory: `tutorial/tutorial25/` *Files can be downlowaded from [here]{}*

-   Input files:

    -   `GaAs.scf` *The `pwscf` input file for ground state calculation*

    -   `GaAs.nscf` *The `pwscf` input file to obtain Bloch states on a
        uniform grid*

    -   `GaAs.pw2wan` *The input file for* `pw2wannier90`

    -   `GaAs.win` *The* `wannier90` *and* `postw90` *input file*

    1.  Run `pwscf` to obtain the ground state of Gallium Arsenide

        `pw.x < GaAs.scf > scf.out`

    2.  Run `pwscf` to obtain the ground state of Gallium Arsenide

        `pw.x < GaAs.nscf > nscf.out`

    3.  Run `Wannier90` to generate a list of the required overlaps
        (written into the `GaAs.nnkp` file)

        `wannier90.x -pp GaAs`

    4.  Run `pw2wannier90` to compute:

        -   The overlaps $\langle u_{n\bf{k}}|u_{n\bf{k+b}}\rangle$
            between spinor Bloch states (written in the `GaAs.mmn` file)

        -   The projections for the starting guess (written in the
            `GaAs.amn` file)

        `pw2wannier90.x < GaAs.pw2wan > pw2wan.out`

    5.  Run `wannier90` to compute MLWFs

        `wannier90.x GaAs`

    6.  Run `postw90` to compute nonlinear shift current

        `postw90.x GaAs` (serial execution)

        `mpirun -np 8 postw90.x GaAs` (example of parallel execution
        with 8 MPI processes)

## Shift current $\sigma^{abc}$ {#shift-current-sigmaabc .unnumbered}

The shift current tensor of GaAs has only one independent component that
is finite, namely $\sigma^{xyz}$. For its computation, set

    berry = true
    berry_task = sc

Like the linear optical conductivity, the shift current is a
frequency-dependent quantity. The frequency window and step is
controlled by `kubo_freq_min`, `kubo_freq_max` and `kubo_freq_step`, as
explained in the users guide.

The shift current requires an integral over the Brillouin zone. The
interpolated k-mesh is controlled by `berry_kmesh`, which has been set
to

    berry_kmesh = 100 100 100

We also need to input the value of the Fermi level in eV:

    fermi_energy = [insert your value here]

Due to the sum over intermediate states involved in the calculation of
the shift current, one needs to consider a small broadening parameter to
avoid numerical problems due to possible degeneracies (see parameter
$\eta$ in Eq. (36) of Ref. [@ibanez-azpiroz_ab_2018] and related
discussion). This parameter is controlled by `sc_eta`. It is normally
found that values between 0.01 eV and 0.1 eV yield an stable spectrum.
The default value is set to $0.04$ eV.

Finally, `sc_phase_conv` controls the phase convention used for the
Bloch sums. `sc_phase_conv=1` uses the so-called tight-binding
convention, whereby the Wannier centres are included into the phase,
while `sc_phase_conv=2` leaves the Wannier centres out of the phase.
These two possible conventions are explained in Ref. [@pythtb]. Note
that the overall shift-current spectrum does not depend on the chosen
convention, but the individual terms that compose it do.

On output, the program generates a set of 18 files named
`SEED-sc_***.dat`, which correspond to the different tensor components
of the shift current (note that the 9 remaining components until
totaling $3\times3\times3=27$ can be obtained from the 18 outputed by
taking into account that $\sigma^{abc}$ is symmetric under
$b\leftrightarrow c$ index exchange). For plotting the only finite
shift-current component of GaAs $\sigma^{xyz}$ (units of A/V$^{2}$) as
in the upper panel of Fig. 3 in Ref. [@ibanez-azpiroz_ab_2018],

    myshell> gnuplot
    gnuplot> plot 'GaAs-sc_xyz.dat' u 1:2 w l


