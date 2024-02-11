#include "wannier.hh"
#include <ISO_Fortran_binding.h>
#include <cassert>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <string>

using namespace std;

const int nb = 12;
const int nw = 8;
const int nk = 64;
const int nn = 8;

void reada(string fn, complex<double> amat[nk][nw][nb]) {
        string junk;
        ifstream infile(fn.c_str(), ios::in);
        getline(infile, junk);
        int nbl, nwl, nkl;
        infile >> nbl >> nkl >> nwl; // you could compare these to expectations

        double a, b;
        for (int j = 0; j < nk; ++j) {
                for (int k = 0; k < nw; ++k) {
                        for (int i = 0; i < nb; ++i) {
                                infile >> nbl >> nwl >> nkl >> a >> b;
                                assert(i == nbl - 1);
                                assert(j == nkl - 1);
                                assert(k == nwl - 1);
                                amat[j][k][i] = complex<double>(a, b);
                        }
                }
        }
}

void readm2(string fn, complex<double> mmat[nk][nk][nb][nb]) {
        string junk;
        ifstream infile(fn.c_str(), ios::in);
        getline(infile, junk);
        int nbl, nnl, nkl, nk2, nnm, nnn; // nnl is the number of neighbours in FD scheme
        infile >> nbl >> nkl >> nnl;      // you could compare these to expectations
        double a, b;
        for (int j = 0; j < nk; ++j) {
                for (int n = 0; n < nn; ++n) {
                        infile >> nkl >> nk2 >> nnl >> nnm >> nnn;
                        assert(nkl == j + 1);
                        int tgk = nk2 - 1;
                        assert(tgk >= 0);
                        assert(tgk < nk);
                        for (int k = 0; k < nb; ++k) {
                                for (int i = 0; i < nb; ++i) {
                                        infile >> a >> b;
                                        mmat[j][tgk][k][i] = complex<double>(a, b);
                                }
                        }
                }
        }
        infile.close();
}

void reade(string fn, double eval[nk][nb]) {
        string junk;
        ifstream infile(fn.c_str(), ios::in);
        int il, jl;
        for (int j = 0; j < nk; ++j) {
                for (int i = 0; i < nb; ++i) {
                        infile >> il >> jl >> eval[j][i];
                        assert(i == il - 1);
                        assert(j == jl - 1);
                }
        }
        infile.close();
}

int main(int argc, char* argv[]) {

        MPI_Init(&argc, &argv);

        int mpisize, mpirank;
        MPI_Comm comm = MPI_COMM_WORLD;
        MPI_Comm_size(comm, &mpisize);
        MPI_Comm_rank(comm, &mpirank);

        // necessary input data
        double kpt[64][3];
        double kptscrambled[64][3];
        int mp_grid[3] = { 4, 4, 4 };
        double cell[3][3] = { { -2.8258062938705995, 0, 2.8258062938705995 }, { 0, 2.8258062938705995, 2.8258062938705995 }, { -2.8258062938705995, 2.8258062938705995, 0 } };
        int excl[5] = { 1, 2, 3, 4, 5 };

        complex<double> amat[nk][nw][nb];
        complex<double> morigl[nk][nn][nb][nb];
        complex<double> morig[nk][nn][nb][nb];
        complex<double> mbig[nk][nk][nb][nb];
        complex<double> umat[nk][nw][nw];
        complex<double> uopt[nk][nw][nb];
        double eval[nk][nb];

        complex<double> samat[nk][nw][nb];
        complex<double> sumat[nk][nw][nw];
        complex<double> suopt[nk][nw][nb];
        double seval[nk][nb];
        int nnkp[nn][nk];

        reada("gaas.amn", amat);
        reade("gaas.eig", eval);
        readm2("gaas.mmn", mbig);

        // kpoint distribution
        int distk[nk];
        int nkl = nk / mpisize;
        if (nk % mpisize > 0) nkl++;
        for (int i = 0; i < nk; ++i) distk[i] = i / nkl;

        // kpoint vectors in w90 order
        int ctr = 0;
        for (int ika = 0; ika < mp_grid[0]; ++ika) {
                for (int ikb = 0; ikb < mp_grid[1]; ++ikb) {
                        for (int ikc = 0; ikc < mp_grid[2]; ++ikc) {
                                kpt[ctr][0] = double(ika) / double(mp_grid[0]);
                                kpt[ctr][1] = double(ikb) / double(mp_grid[1]);
                                kpt[ctr][2] = double(ikc) / double(mp_grid[2]);
                                ctr++;
                        }
                }
        }

        // scramble everything
        int random[64];
        for (int i = 0; i < nk; i++) random[i] = i;
        bool lrand = true;
        if (lrand) {
                for (int i = 0; i < nk; i++) {
                        int t = rand() % nk;
                        int x = random[i];
                        random[i] = random[t];
                        random[t] = x;
                }
        }
        for (int i = 0; i < nk; i++) {
                bool found = false;
                for (int j = 0; j < nk; j++) {
                        if (j == random[i] && !found) {
                                found = true;
                        } else if (j == random[i] && found) {
                                cerr << "nonuniq" << endl;
                                exit(9);
                        }
                }
                if (!found) {
                        cerr << "missing" << endl;
                        exit(9);
                }
        }

        for (int i = 0; i < nk; ++i) {
                kptscrambled[i][0] = kpt[random[i]][0];
                kptscrambled[i][1] = kpt[random[i]][1];
                kptscrambled[i][2] = kpt[random[i]][2];
                for (int j = 0; j < nw; ++j) {
                        for (int k = 0; k < nb; ++k) {
                                samat[i][j][k] = amat[random[i]][j][k];
                                suopt[i][j][k] = uopt[random[i]][j][k];
                        }
                }
                for (int j = 0; j < nw; ++j) {
                        for (int k = 0; k < nw; ++k) {
                                sumat[i][j][k] = umat[random[i]][j][k];
                        }
                }
                for (int k = 0; k < nb; ++k) {
                        seval[i][k] = eval[random[i]][k];
                }
        }

        // setup
        void* w90glob;
        int nntmp; // because nn above is const
        wannier_setup(w90glob, &distk[0], &kptscrambled[0][0], &mp_grid[0], nk, nb, nw, &cell[0][0], &excl[0], nntmp, &nnkp[0][0], comm);
        assert(nn == nntmp);

        for (int ik = 0; ik < nk; ++ik) {
                for (int in = 0; in < nn; ++in) {
                        for (int ib = 0; ib < nb; ++ib) {
                                for (int jb = 0; jb < nb; ++jb) {
                                        int sik = random[ik];                           // shuffled k index
                                        int sjk = random[nnkp[in][ik] - 1];             // shuffled neighbour k index
                                        morig[ik][in][ib][jb] = mbig[sik][sjk][ib][jb]; // mbig <u_k|u_k'>
                                }
                        }
                }
        }

        // copy appropriate k points to local array (decompose m)
        ctr = 0;
        for (int ik = 0; ik < nk; ++ik) {
                if (distk[ik] == mpirank) {
                        for (int in = 0; in < nn; ++in) {
                                for (int ib = 0; ib < nb; ++ib) {
                                        for (int jb = 0; jb < nb; ++jb) {
                                                morigl[ctr][in][ib][jb] = morig[ik][in][ib][jb];
                                        }
                                }
                        }
                        ctr++;
                }
        }

        // basic arrays for the resulting spreads for printout
        double wannier_ctr[nw][3];
        double wannier_spr[nw];

        wannier_run(w90glob, &samat[0][0][0], &seval[0][0], &morigl[0][0][0][0], &sumat[0][0][0], &wannier_ctr[0][0], &wannier_spr[0]);

        if (mpirank == 0) {
                cout << fixed << setprecision(10);
                for (int iw = 0; iw < nw; ++iw) {
                        cout << setw(20) << wannier_ctr[iw][0];
                        cout << setw(20) << wannier_ctr[iw][1];
                        cout << setw(20) << wannier_ctr[iw][2];
                        cout << setw(20) << wannier_spr[iw];
                        cout << endl;
                }
        }

        MPI_Finalize();
        return 0;
}
