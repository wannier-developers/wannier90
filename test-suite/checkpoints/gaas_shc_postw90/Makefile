all: $(patsubst %.win, %.chk.fmt.bz2, $(wildcard *.win))

%.chk.fmt.bz2: %.win %.mmn %.amn
	$(eval SEEDNAME:=$(patsubst %.win, %, $<))
	echo $(SEEDNAME)
	../../../wannier90.x $< && ../../../w90chk2chk.x -u2f $(SEEDNAME) && bzip2 $(SEEDNAME).chk.fmt && rm $(SEEDNAME).chk

%.mmn: %.mmn.bz2
	$(eval SEEDNAME:=$(patsubst %.mmn.bz2, %, $<))
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).mmn

%.amn: %.amn.bz2
	$(eval SEEDNAME:=$(patsubst %.amn.bz2, %, $<))
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).amn

.PHONY: all
