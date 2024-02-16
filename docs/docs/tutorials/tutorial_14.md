# 14: Linear Sodium Chain &#151; Transport properties {#linear-sodium-chain-transport-properties .unnumbered}

-   Outline: *Compare the quantum conductance of a periodic linear chain
    of Sodium atoms with that of a defected chain*

-   Directories: `tutorials/tutorial14/periodic`<br>
    &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;  `tutorials/tutorial14/defected`<br> 
    *Files can be downloaded from [here](https://github.com/wannier-developers/wannier90/tree/develop/tutorials/tutorial14)*

-   Input Files

    -    `Na_chain.scf` *The `pwscf` input file for ground
        state calculation*

    -    `Na_chain.nscf` *The `pwscf` input file to obtain
        Bloch states for the conduction states*

    -    `Na_chain.pw2wan` *Input file for `pw2wannier90`*

    -    `Na_chain.win` *The `wannier90` input file*

The periodic system contains two unit cells evenly distributed along the
supercell. Transport calculations are performed using
` transport_mode = bulk` and so the resulting quantum conductance
represents that of an infinite periodic chain.

The part of the `wannier90` input file that controls the transport part
of the calculation looks like:


```vi title="Input file"
transport = true
transport_mode = bulk
tran_read_ht = false
one_dim_axis = x
fermi_energy = -2.7401
tran_win_min = -5.0
tran_win_max = 5.0
tran_energy_step = 0.01
translation_centre_frac = 0.5 0.5 0.5
tran_num_bb = 2
```

The defected system uses a 13 atom supercell with the central atom
position altered to break symmetry. Setting `transport_mode = lcr` with
tell `wannier90` to treat the system as an infinite system with the
defect at its centre. The supercell is chosen so that is conforms to the
2c2 geometry (see User Guide for details). Each principal layer is 2
atoms long so that the conductor region contains the defected atom plus
a single atom on either side.

The transport section of the input file contains these key differences:

```vi title="Input file"
transport_mode = lcr
tran_num_ll = 2
tran_num_cell_ll = 2
```

Descriptions of these and other keywords related to the calculation of
transport properties can be found in the User Guide.

1.  Run `pwscf` and `wannier90` for the periodic system.

2.  Run `pwscf` and `wannier90` for the defected system.

3.  The quantum conductance is written to the files
    `periodic/Na_chain_qc.dat` and , respectively. Compare the quantum
    conductance of the periodic (bulk) calculation with the defected
    (LCR) calculation. Your plot should look like
    [this](#fig:Na_qc){reference-type="ref" reference="fig:Na_qc"}.


<figure markdown="span" id="fig:Na_qc">
![Image title](./Na_qc.webp){ width="500" }
<figcaption> Quantum conductance of periodic Sodium chain (black)
compared to that of the defected Sodium chain (red).</figcaption>
</figure>


