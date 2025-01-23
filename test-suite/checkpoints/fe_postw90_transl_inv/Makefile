all: $(patsubst %.win, %.chk.fmt.bz2, $(wildcard *.win))

%.chk.fmt.bz2: %.win
	$(eval SEEDNAME:=$(patsubst %.win, %, $<))
	echo $(SEEDNAME)
	../../../wannier90.x $< && ../../../w90chk2chk.x -u2f $(SEEDNAME) && bzip2 $(SEEDNAME).chk.fmt && rm $(SEEDNAME).chk

.PHONY: all