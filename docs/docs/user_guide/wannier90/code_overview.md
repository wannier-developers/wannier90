# Code Overview

`wannier90` can operate in two modes:

1.  *Post-processing mode:* read in the overlaps and projections from
    file as computed inside a first-principles code. We expect this to
    be the most common route to using `wannier90`, and is described in
    Ch. [5](#ch:wann-pp){reference-type="ref" reference="ch:wann-pp"};

2.  *Library mode:* as a set of library routines to be called from
    within a first-principles code that passes the overlaps and
    projections to the `wannier90` library routines and in return gets
    the unitary transformation corresponding to MLWF. This route should
    be used if the MLWF are needed within the first-principles code, for
    example in post-LDA methods such as LDA+U or SIC, and is described
    in Ch. [6](#ch:wann-lib){reference-type="ref"
    reference="ch:wann-lib"}.

<figure id="structure">
<div class="center">
<embed src="overview.pdf" style="width:5in" />
</div>
<figcaption>Schematic overview of the module structure of
<code>wannier90</code>. Modules may only use data and subroutines from
lower modules.</figcaption>
</figure>
