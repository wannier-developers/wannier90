# Wannier90: A Guide for Contributors

As the community of developers grows, it is important to maintain a degree of consistency in the distribution in terms of coding style, documentation, tutorial examples and tests cases for all aspects of the functionality of the code.

## Coding style
* Code should conform to the Fortran 2003 standard, but avoid using new object-oriented programming (OOP) aspects, because some compilers do not yet support them fully
* Any MPI code should conform to version 2.x of the MPI standard
* `USE` statements should always have the `ONLY` option with an explicit list of variables
* Static variables should have the `SAVE` attribute
* Type conversions should be explicit
* Each source file should be a module and should:
  - contain the `IMPLICIT NONE` statement
  - depend on as few other modules as possible
  - be declared `PRIVATE` and have a minimal set of public entities
  - start with the GPL header that is common to all module files
  - be properly annotated with Fortran Documenter (FORD) comments
* Subroutines that are public should be named with a common prefix associated with the module in which they sit. E.g., the names of public subroutines in `wannierise_mod.F90` start with `wann_`, etc.
* Variable, parameter, function, subroutine and module names should be informative and, in so far as possible, maintain consistency of style with existing names

## Documentation and tutorial examples
There are a number of different types of documentation associated with the distribution and all relevant ones should be updated and extended to reflect modifications and additions to the code.
* `CHANGE.log` in the root directory of the distribution describes the updates with respect to the most recent release
* `/doc/user_guide/` is the main User Guide for the code, where new variables, input parameters, functionality and file formats should be described
* `/doc/tutorial/` is the tutorial guide for the examples in /examples/ and should be updated whenever a new example is added
* `/examples/README` provides a very brief description of each example and the associated functionality that it covers
FORD annotations should be included in all code that is developed

## Test suite
A set of tests is provided with Wannier90, in the folder `test-suite`.
If you add new functionality, **it is required that you also add the respective tests** to make sure the functionality always behaves as expected. 

Try to add only tests that run within a few seconds. In most cases, this is possible (tests do not need to be converged, just need to check that the functionality is working as expected).

Before committing, please check that the code compiles and that the tests run for you. To know how to write a test, and how to run them, read the README file inside the `test-suite` folder.

Also, when you create a pull request, Travis-CI will run the same tests and show a green tick or a red cross depending on whether all the tests (that do not require the interface) pass. This typically takes just a few minutes, so after you create a pull request please check that all tests have passed.

# GitHub repository and management of pull requests

The Git repository for the development of the Wannier90 code is hosted on GitHub:
https://github.com/wannier-developers/wannier90

We adopt the branching model described here:
http://nvie.com/posts/a-successful-git-branching-model/

The main development branch is called `develop`. Typically, a developer:
* creates a fork of the project
* creates a new branch, branching off from develop, for each different feature, and makes modifications on their fork/branch
 
  *Note*: using different branches is important. The reason is the following: if you create a pull request, until when it is accepted, if you continue committing in the same branch, the pull request gets updated. Now, if you work at the same time on two features, and one of the two is not yet acceptable, we cannot accept the pull request (so if they come in the same request, neither of them gets accepted). So please work in different branches, and do different pull requests for different features.
* When ready to merge changes back, a pull request is made on `develop`.

While working on a branch/fork, merge the develop branch back into it often in order to stay in synchronisation with the current development version of the code and to minimise potential merge conflicts. For the same reason, commit very often your work in small units (use different commits for different functionality, and also commit often in your fork while you work) - this also reduces the risk of merge conflicts.

*Note on the `master` branch*: The `master` branch points to the most recent release of the code. Commits and pull requests should *not* be done on `master`.

New official releases of the code will be tagged as `vX.Y` (major releases) or `vX.Y.Z` (minor releases).

## Necessary conditions for a pull request to be accepted
Pull requests made on the `develop` branch will go to the Wannier Developers' Group for consideration. Here is a checklist of necessary conditions that must be met for a pull request to be accepted:

* The coding style has been adopted
* `CHANGE.log` has been updated 
* `/doc/user_guide/` has been updated (e.g., if new input parameters, new functionality, or a new input/output files have been added)
* At least an example has been added to the tutorial set in `/examples`, and `/examples/README` and `/doc/tutorial/` have been updated accordingly (e.g., if new functionality has been added)
* A test case has been added to `/test-suite/` (e.g., if new functionality has been added)
* The code compiles and passes the set of tests in the test suite (which runs automatically when a pull request is made)

Whilst we expect to accept the majority of pull requests, it is possible that we will not accept all pull requests. In such cases we will always endeavour to explain to the developer our reasons for not doing so.

## List of contributors
Prior to each release of the code, the Wannier Developers’ Group will update the list of contributors in the `README` file in the root of the distribution so that contributions are appropriately attributed and contributors receive due recognition for their contributions. 
