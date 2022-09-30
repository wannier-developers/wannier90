all: $(patsubst %.chk.fmt.bz2, %.chk, $(wildcard *.chk.fmt.bz2)) $(patsubst %.mmn.bz2, %.mmn, $(wildcard *.mmn.bz2)) $(patsubst %.uHu.bz2, %.uHu, $(wildcard *.uHu.bz2))

%.chk: %.chk.fmt.bz2
	$(eval SEEDNAME:=$(patsubst %.chk.fmt.bz2, %, $<))
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).chk.fmt && ../../../w90chk2chk.x -f2u $(SEEDNAME) && rm $(SEEDNAME).chk.fmt

%.uHu: %.uHu.bz2
	$(eval SEEDNAME:=$(patsubst %.uHu.bz2, %, $<))
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).uHu 

%.mmn: %.mmn.bz2
	$(eval SEEDNAME:=$(patsubst %.mmn.bz2, %, $<))
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).mmn


.PHONY: all
