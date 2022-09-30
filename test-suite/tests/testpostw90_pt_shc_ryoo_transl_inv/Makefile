all: $(patsubst %.chk.fmt.bz2, %.chk, $(wildcard *.chk.fmt.bz2)) $(patsubst %.mmn.bz2, %.mmn, $(wildcard *.mmn.bz2)) $(patsubst %.amn.bz2, %.amn, $(wildcard *.amn.bz2)) $(patsubst %.spn.bz2, %.spn, $(wildcard *.spn.bz2)) $(patsubst %.sHu.bz2, %.sHu, $(wildcard *.sHu.bz2)) $(patsubst %.sIu.bz2, %.sIu, $(wildcard *.sIu.bz2))

%.chk: %.chk.fmt.bz2
	$(eval SEEDNAME:=$(patsubst %.chk.fmt.bz2, %, $<))
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).chk.fmt && ../../../w90chk2chk.x -f2u $(SEEDNAME) && rm $(SEEDNAME).chk.fmt

%.mmn: %.mmn.bz2
	$(eval SEEDNAME:=$(patsubst %.mmn.bz2, %, $<))
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).mmn

%.amn: %.amn.bz2
	$(eval SEEDNAME:=$(patsubst %.amn.bz2, %, $<))
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).amn

%.spn: %.spn.bz2
	$(eval SEEDNAME:=$(patsubst %.spn.bz2, %, $<))
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).spn

%.sHu: %.sHu.bz2
	$(eval SEEDNAME:=$(patsubst %.sHu.bz2, %, $<))
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).sHu

%.sIu: %.sIu.bz2
	$(eval SEEDNAME:=$(patsubst %.sIu.bz2, %, $<))
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).sIu

.PHONY: all