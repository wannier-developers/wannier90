Wannier90 autodoc
=================

The auto-generated documentation is built using FORD
(https://github.com/Fortran-FOSS-Programmers/ford), whose documentation is at
https://forddocs.readthedocs.io.

The published result is at https://wannier-developers.github.io/wannier90/. It is
built and deployed automatically by the GitHub Actions workflows in
.github/workflows/ (step_autodoc.yaml and autodoc.yaml), so it does not normally
need to be built by hand. The manual route below remains useful for previewing
changes locally.

Installing FORD
---------------

FORD must be installed in an environment SEPARATE from the one used for the
MkDocs user documentation in docs/requirements.txt. FORD pins `markdown ~= 3.4`
while mkdocs-material 9.7.7 requires `markdown >= 3.6`, so installing both
together fails with a pip resolver conflict. Use a dedicated virtualenv:

    python3 -m venv .venv-ford
    .venv-ford/bin/pip install ford

Graphviz must also be installed on the system, otherwise `graph: true` in
project.md cannot generate the call and dependency graphs.

Building locally
----------------

To create the documentation, run

    ford project.md

The main page of the documentation will be ./build/index.html

For working on the autodoc, you can use

    ford fast.md

which will not create the graphs and search indices, and thus be faster. The main page will then be ./build_fast/index.html
