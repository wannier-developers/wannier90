
#include "wannier90.hh"
#include <complex>

void wannier_setup(void*& w90glob, double* kpts, int* mp_grid, int nk, int nb, int nw, double* cell, int& nnfd, int* nnkp) {
        w90glob = w90_create(); // allocate an instance of the library data block

        cset_option(w90glob, "kpoints", kpts, nk, 3);
        cset_option(w90glob, "mp_grid", mp_grid);
        cset_option(w90glob, "num_bands", nb);
        cset_option(w90glob, "num_kpts", nk);
        cset_option(w90glob, "num_wann", nw);
        cset_option(w90glob, "unit_cell_cart", cell, 3, 3); // because cell is not [][3] now

        int ierr;
        cinput_setopt(w90glob, "gaas", ierr); // process necessary library options
        if (ierr != 0 ) exit(ierr);
        cinput_reader(w90glob, ierr); // process any other options
        if (ierr != 0 ) exit(ierr);

        cget_nn(w90glob, nnfd);  // return number of NN in FD scheme
        cget_nnkp(w90glob, nnkp); // return indexes of NN k-points in FD scheme
}

void wannier_run(void* w90glob, std::complex<double>* amat, double* eval, std::complex<double>* morigl, std::complex<double>* umat, double* wannier_ctr, double* wannier_spr) {

        cset_eigval(w90glob, eval); // contains eigenvalues
        cset_m_local(w90glob,morigl); // m matrix
        cset_u_opt(w90glob, amat); // initial projections
        cset_u_matrix(w90glob, umat); // results returned here

        int ierr;
        cdisentangle(w90glob, ierr);
        if (ierr != 0 ) exit(ierr);
        cwannierise(w90glob, ierr);
        if (ierr != 0 ) exit(ierr);
        cget_centres(w90glob, wannier_ctr);
        cget_spreads(w90glob, wannier_spr);
        w90_delete(w90glob);
}
