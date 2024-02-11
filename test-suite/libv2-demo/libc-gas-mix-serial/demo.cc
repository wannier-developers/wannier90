#include <ISO_Fortran_binding.h>
#include <cassert>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "wannier.hh"

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

void readm(string fn, complex<double> mmat[nk][nk][nb][nb]) {
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

        // necessary input data
        double kpt[64][3];
        int mp_grid[3] = { 4, 4, 4 };
        double cell[3][3] = { { -2.8258062938705995, 0, 2.8258062938705995 }, { 0, 2.8258062938705995, 2.8258062938705995 }, { -2.8258062938705995, 2.8258062938705995, 0 } };

        complex<double> amat[nk][nw][nb];
        complex<double> morigl[nk][nn][nb][nb];
        complex<double> mbig[nk][nk][nb][nb];
        complex<double> umat[nk][nw][nw]; // normally nk,nw,nw
        double eval[nk][nb];
        int nnkp[nn][nk];

        reada("gaas.amn", amat);
        reade("gaas.eig", eval);
        readm("gaas.mmn", mbig);

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


        //setup
        void* w90glob;
        int nntmp; // because nn above is const
        wannier_setup(w90glob, &kpt[0][0], &mp_grid[0], nk, nb, nw, &cell[0][0], nntmp, &nnkp[0][0]);
        assert(nn == nntmp);

        // copy appropriate k points to local array (decompose m)
        for (int ik = 0; ik < nk; ++ik) {
                for (int in = 0; in < nn; ++in) {
                        for (int ib = 0; ib < nb; ++ib) {
                                for (int jb = 0; jb < nb; ++jb) {
                                        int jk = nnkp[in][ik] - 1;
                                        morigl[ik][in][ib][jb] = mbig[ik][jk][ib][jb];
                                }
                        }
                }
        }

        double wannier_ctr[nw][3];
        double wannier_spr[nw];

        wannier_run(w90glob, &amat[0][0][0], &eval[0][0], &morigl[0][0][0][0], &umat[0][0][0], &wannier_ctr[0][0], &wannier_spr[0]);

        cout << fixed << setprecision(10);
        for (int iw = 0; iw < nw; ++iw) {
                cout << setw(20) << wannier_ctr[iw][0];
                cout << setw(20) << wannier_ctr[iw][1];
                cout << setw(20) << wannier_ctr[iw][2];
                cout << setw(20) << wannier_spr[iw];
                cout << endl;
        }
        return 0;
}
