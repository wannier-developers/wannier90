all: $(patsubst %.chk.fmt.bz2, %.chk, $(wildcard *.chk.fmt.bz2))

%.chk: %.chk.fmt.bz2
	$(eval SEEDNAME:=$(patsubst %.chk.fmt.bz2, %, $<))
	echo $(SEEDNAME)
	cat $< | bunzip2 - > $(SEEDNAME).chk.fmt && ../../../w90chk2chk.x -f2u $(SEEDNAME) && rm $(SEEDNAME).chk.fmt

.PHONY: all