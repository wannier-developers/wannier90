all: $(patsubst %.mmn.bz2, %.mmn, $(wildcard *.mmn.bz2)) $(patsubst %.amn.bz2, %.amn, $(wildcard *.amn.bz2))

%.mmn: %.mmn.bz2
	$(eval SEEDNAME:=$(patsubst %.mmn.bz2, %, $<))
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).mmn

%.amn: %.amn.bz2
	$(eval SEEDNAME:=$(patsubst %.amn.bz2, %, $<))
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).amn


.PHONY: all
