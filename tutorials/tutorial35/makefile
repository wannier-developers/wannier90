outfiles := bands.out \
bandsx.out \
silicon.bands.dat \
silicon.bands.dat.gnu \
silicon.bands.dat.proj \
silicon.bands.dat.proj.projwfc_up \
silicon.bands.dat.rap \
silicon.pdos_atm\#1\(C\)_wfc\#1\(s\) \
silicon.pdos_atm\#1\(C\)_wfc\#2\(p\) \
silicon.pdos_atm\#2\(C\)_wfc\#1\(s\) \
silicon.pdos_atm\#2\(C\)_wfc\#2\(p\) \
silicon.pdos_tot \
silicon.plotband.gnu \
silicon.plotband.pdf \
silicon.plotband.ps \
silicon.projbands.gnu \
silicon.projbands.gnu_projected.pdf \
silicon.projbands.gnu_projected.ps \
out/ \
plotband.out \
projwfc.out \
scf.out \
CRASH \
silicon.plotband.gnu.dat \
silicon.werr \
silicon.wout \
nscf.out \
pw2wan.out \
silicon.amn \
silicon.eig \
silicon.mmn \
silicon.node_00000.werr \
silicon_bandsdiff.pdf

infiles := silicon.bands \
silicon.bandsx \
silicon.nscf \
silicon.plotband \
silicon.projwfc \
silicon.pw2wan \
silicon.scf \
silicon.win

default:
	@echo "Usage: make clean"

clean:
	echo "Cleaning up..."
	rm -r $(outfiles)

# Usually, it is mpirun. If not set, run in serial.
# MPI_CMD = mpirun -n 8

test:
	@echo "Testing all calculations"
	$(MPI_CMD) pw.x < silicon.scf > scf.out
	$(MPI_CMD) pw.x < silicon.bands > bands.out
	$(MPI_CMD) bands.x < silicon.bandsx > bandsx.out
	$(MPI_CMD) pw.x < silicon.nscf > nscf.out
	$(MPI_CMD) wannier90.x -pp silicon
	$(MPI_CMD) pw2wannier90.x < silicon.pw2wan > pw2wan.out
	$(MPI_CMD) wannier90.x silicon
	./silicon_bandsdiff.gnu
	@echo "All calculations finished successfully"
