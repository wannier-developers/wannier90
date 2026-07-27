//
// Copyright (C) 2026 Wannier Developer Group
//
// This library is free software; you can redistribute it
// and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software
// Foundation; either version 2.1 of the License, or (at your
// option) any later version.
//
// This library is distributed in the hope that it will be
// useful,but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.  See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General
// Public License along with this library; if not, see
// <https://www.gnu.org/licenses/>.
//
// The webpage of the Wannier90 code is
// <https://www.wannier.org>.
//
// The Wannier90 code is hosted on GitHub
// <https://github.com/wannier-developers/wannier90>
//

extern "C" {
        #include "wannier90.h"
}
#include <cassert>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>

using namespace std;

// functions for reading .mmn, .amn, .eig and certain values from .win file
int nband(string);
int nwann(string);
void kpts(double[][3], int, string);
void mpgrid(int[3], string);
void rlat(double[][3], string);
void readm(string, int**, int***, int, int, int, complex<double>*);
void reada(string, int, int, int, complex<double>*);
void reade(string, int, int, double*);

int main(int argc, char* argv[]) {

#ifdef MPI
        // before argv parsing: MPI_Init may consume launcher-added arguments
        MPI_Init(&argc, &argv);
#endif
        // check command-line argument and input file existence
        if ((argc != 2) || (!filesystem::exists(argv[1]))) {
                cerr << "usage: " << argv[0] << " xxx.win" << endl;
                exit(1);
        }

        // .win file (argument 1) as a string
        ifstream winfile(argv[1], ios::in);
        stringstream filestream;
        filestream << winfile.rdbuf();
        string filestring { filestream.str() };
        winfile.close();

        // clean up win format
        regex regexp { R"([#!].*)" }; // remove comments
        filestring = regex_replace(filestring, regexp, string(""));
        regexp = { R"(:|=)" }; // remove define/assign
        filestring = regex_replace(filestring, regexp, string(""));
        regexp = { R"(^ *)" }; // remove define/assign
        filestring = regex_replace(filestring, regexp, string(""));

        regexp = { R"(^(.*)\.win$)" };
        smatch match;
        string fn { argv[1] };
        regex_search(fn, match, regexp);
        assert(!match.empty());
        string root = match[1];

        double uc[3][3];
        rlat(uc, filestring);

        int nkabc[3];
        mpgrid(nkabc, filestring);
        int nk = nkabc[0] * nkabc[1] * nkabc[2];
        double(*kpt)[3] = new double[nk][3];

        kpts(kpt, nk, filestring);

        int nw = nwann(filestring);
        int nb = nband(filestring);
        if (nb == 0) nb = nw; // not specified in input



        ///////////// LIBRARY

        w90_data w90glob;
        w90_create(&w90glob);
#ifdef MPI
        w90_set_comm(w90glob, MPI_COMM_WORLD);
#endif

        w90_set_option_double2d(w90glob, "kpoints", &kpt[0][0], nk, 3);
        w90_set_option_int1d(w90glob, "mp_grid", nkabc, 3);
        w90_set_option_int(w90glob, "num_bands", nb);
        w90_set_option_int(w90glob, "num_wann", nw);
        w90_set_option_double2d(w90glob, "unit_cell_cart", &uc[0][0], 3, 3);
        w90_set_option_logical(w90glob, "dump_inputs", true);

        int ierr;
        w90_input_setopt(w90glob, root.c_str(), &ierr); // process necessary library options
        assert(ierr == 0);

        w90_input_reader(w90glob, &ierr); // process any other options
        assert(ierr == 0);

        int nnfd;
        w90_get_nn(w90glob, &nnfd, &ierr); // return number of NN in FD scheme
        assert(ierr == 0);
        cout << nnfd <<endl;

        // prepare nnkp and gkpb arrans
        // horrible contortions to make multidimensional arrays...
        int* nndata = new int[nk * nnfd];
        int** nnkp = new int*[nnfd];
        for (int i = 0; i < nnfd; ++i) {
                nnkp[i] = &nndata[i * nk];
        }

        int* gkpbdata = new int[nk * nnfd * 3];
        int*** gkpb = new int**[nnfd];
        for (int j = 0; j < nnfd; ++j) {
                gkpb[j] = new int*[nk * nnfd];

                for (int i = 0; i < nk; ++i) {
                        gkpb[j][i] = new int[3];
                        gkpb[j][i] = &gkpbdata[j * nk * 3 + i * 3];
                }
        }

        w90_get_nnkp(w90glob, &nnkp[0][0], &ierr); // return indexes of NN k-points in FD scheme
        assert(ierr == 0);
        w90_get_gkpb(w90glob, &gkpb[0][0][0], &ierr);
        assert(ierr == 0);

        // printout nnkp data for testing
        for (int i = 0; i < nk; ++i) {
                for (int j = 0; j < nnfd; ++j) {
                        cout << "NNKP: " << setw(4) << i << setw(4) << nnkp[j][i]
                             << setw(4) << gkpb[j][i][0]
                             << setw(4) << gkpb[j][i][1]
                             << setw(4) << gkpb[j][i][2] << endl;
                }
        }

        // prepare main data arrays
        complex<double>* adata = new complex<double>[nb * nw * nk];
        complex<double>* mdata = new complex<double>[nb * nb * nk * nnfd];
        complex<double>* umat = new complex<double>[nw * nw * nk];
        double* edata = new double[nb * nk];

        w90_set_m_local(w90glob, reinterpret_cast<_Complex double*>(mdata)); // m matrix
        w90_set_u_matrix(w90glob, reinterpret_cast<_Complex double*>(umat)); // results returned here
        w90_set_u_opt(w90glob, reinterpret_cast<_Complex double*>(adata));   // initial projections
        w90_set_eigval(w90glob, edata);  // contains eigenvalues

        fn = root + ".mmn";
        readm(fn, nnkp, gkpb, nk, nb, nnfd, mdata);

        fn = root + ".amn";
        reada(fn, nk, nb, nw, adata);

        fn = root + ".eig";
        if (filesystem::exists(fn)) reade(fn, nk, nb, edata);

        w90_disentangle(w90glob, &ierr);
        w90_project_overlap(w90glob, &ierr);
        assert(ierr == 0);

        w90_wannierise(w90glob, &ierr);
        assert(ierr == 0);

        double wannier_ctr[nw][3];
        double wannier_spr[nw];
        w90_get_centres(w90glob, &wannier_ctr[0][0]);
        w90_get_spreads(w90glob, wannier_spr);
        w90_delete(&w90glob);

#ifdef MPI
        MPI_Finalize();
#endif
        return 0;
}

int nwann(string filestring) {
        regex regexp { R"(num_wann\s*(\d+)\s*\n)" };
        smatch match;
        regex_search(filestring, match, regexp);
        assert(!match.empty());
        return stoi(match[1]);
}
int nband(string filestring) {
        regex regexp { R"(num_bands\s*(\d+)\s*\n)" };
        smatch match;
        regex_search(filestring, match, regexp);
        if (!match.empty()) {
                return stoi(match[1]);
        } else {
                return 0;
        }
}

void rlat(double rlat[][3], string filestring) {
        // unit cell specification
        regex latticeregex(R"((begin unit_cell_cart\s*)\n\s*(bohr|ang)?\s*\n?\s*(-?\d*\.\d+)\s*\n?\s*(-?\d*\.\d+)\s*\n?\s*(-?\d*\.\d+)\s*\n?\s*(-?\d*\.\d+)\s*\n?\s*(-?\d*\.\d+)\s*\n?\s*(-?\d*\.\d+)\s*\n?\s*(-?\d*\.\d+)\s*\n?\s*(-?\d*\.\d+)\s*\n?\s*(-?\d*\.\d+)\s*\n?\s*(end unit_cell_cart))", std::regex_constants::icase);
        smatch latticematch;
        regex_search(filestring, latticematch, latticeregex);
        assert(!latticematch.empty());
        double fac = 1.0;
        if (latticematch[2] == "bohr") fac = 0.52917720859;
        rlat[0][0] = stod(latticematch[3]) * fac;
        rlat[0][1] = stod(latticematch[4]) * fac;
        rlat[0][2] = stod(latticematch[5]) * fac;
        rlat[1][0] = stod(latticematch[6]) * fac;
        rlat[1][1] = stod(latticematch[7]) * fac;
        rlat[1][2] = stod(latticematch[8]) * fac;
        rlat[2][0] = stod(latticematch[9]) * fac;
        rlat[2][1] = stod(latticematch[10]) * fac;
        rlat[2][2] = stod(latticematch[11]) * fac;

        /*
        rlat[0][0] = stod(latticematch[3]) * fac;
        rlat[1][0] = stod(latticematch[4]) * fac;
        rlat[2][0] = stod(latticematch[5]) * fac;
        rlat[0][1] = stod(latticematch[6]) * fac;
        rlat[1][1] = stod(latticematch[7]) * fac;
        rlat[2][1] = stod(latticematch[8]) * fac;
        rlat[0][2] = stod(latticematch[9]) * fac;
        rlat[1][2] = stod(latticematch[10]) * fac;
        rlat[2][2] = stod(latticematch[11]) * fac;
        */
}

void mpgrid(int nkabc[3], string filestring) {
        regex nkregex { R"(mp_grid\s*(\d+)\s(\d+)\s(\d+)\s*\n)" };
        smatch nkmatch;
        regex_search(filestring, nkmatch, nkregex);
        assert(!nkmatch.empty());
        nkabc[0] = stoi(nkmatch[1]);
        nkabc[1] = stoi(nkmatch[2]);
        nkabc[2] = stoi(nkmatch[3]);
        assert(nkabc[0] > 0);
        assert(nkabc[1] > 0);
        assert(nkabc[2] > 0);
}

void kpts(double kpts[][3], int nkexpected, string filestring) {
        regex kptregex(R"(begin kpoints\s*\n((.|\n)*)\s*\nend kpoints)", std::regex_constants::icase);
        smatch kptmatch;
        regex_search(filestring, kptmatch, kptregex);
        assert(!kptmatch.empty());
        stringstream kstream(kptmatch[1]);
        int nkfound { 0 };
        string line;
        while (getline(kstream, line)) {
                if (line.empty()) continue;
                stringstream linestream(line);
                double kx, ky, kz;
                if (!(linestream >> kx >> ky >> kz)) continue;
                assert(nkfound < nkexpected);
                kpts[nkfound][0] = kx;
                kpts[nkfound][1] = ky;
                kpts[nkfound][2] = kz;
                nkfound++;
        }
        assert(nkfound == nkexpected);
}

void readm(string filename, int** nnkp, int*** gkpb, int nkexpect, int nbexpect, int nnexpect, complex<double>* m) {
        string header;
        ifstream mfile(filename, ios::in);
        getline(mfile, header);
        int nnfile, nkfile, nbfile;
        mfile >> nbfile >> nkfile >> nnfile;
        assert(nnfile == nnexpect);
        assert(nkfile == nkexpect);
        assert(nbfile == nbexpect);

        for (int ik = 0; ik < nkexpect; ++ik) {
                for (int jk = 0; jk < nnexpect; ++jk) {
                        int tk, tkp, tx, ty, tz;
                        mfile >> tk >> tkp >> tx >> ty >> tz;
                        if (mfile) {

                                // some of the examples have reordered m matrix files
                                int kk;
                                for (kk = 0; kk < nnexpect; ++kk) {
                                        if (nnkp[kk][ik] == tkp && tx == gkpb[kk][ik][0] && ty == gkpb[kk][ik][1] && tz == gkpb[kk][ik][2]) break;
                                }

                                assert(tk == ik + 1);
                                assert(tkp == nnkp[kk][ik]);
                                assert(tx == gkpb[kk][ik][0]);
                                assert(ty == gkpb[kk][ik][1]);
                                assert(tz == gkpb[kk][ik][2]);

                                int ctr = 0;
                                int off = (ik * nnexpect + kk) * nbexpect * nbexpect;
                                for (int ib = 0; ib < nbexpect; ++ib) {
                                        for (int jb = 0; jb < nbexpect; ++jb) {
                                                double r, c;
                                                mfile >> r >> c;
                                                m[off + ctr++] = complex<double>(r, c);
                                        }
                                }
                        }
                }
        }
        mfile.close();
}
void reada(string filename, int nkexpect, int nbexpect, int nwexpect, complex<double>* a) {
        string header;
        ifstream afile(filename, ios::in);
        getline(afile, header);
        int nwfile, nkfile, nbfile;
        afile >> nbfile >> nkfile >> nwfile;
        assert(nwfile == nwexpect);
        assert(nkfile == nkexpect);
        assert(nbfile == nbexpect);

        int ctr { 0 };
        for (int ik = 0; ik < nkexpect; ++ik) {
                for (int iw = 0; iw < nwexpect; ++iw) {
                        for (int ib = 0; ib < nbexpect; ++ib) {
                                double r, c;
                                int tb, tk, tw;
                                afile >> tb >> tw >> tk >> r >> c;
                                if (afile) {
                                        assert(tb == ib + 1);
                                        assert(tk == ik + 1);
                                        assert(tw == iw + 1);
                                        a[ctr++] = complex<double>(r, c);
                                }
                        }
                }
        }
        afile.close();
}
void reade(string filename, int nkexpect, int nbexpect, double* e) {
        string header;
        ifstream efile(filename, ios::in);
        int ctr { 0 };
        for (int ik = 0; ik < nkexpect; ++ik) {
                for (int ib = 0; ib < nbexpect; ++ib) {
                        double r;
                        int tb, tk;
                        efile >> tb >> tk >> r;
                        if (efile) {
                                assert(tb == ib + 1);
                                assert(tk == ik + 1);
                                e[ctr++] = r;
                        }
                }
        }
        efile.close();
}
