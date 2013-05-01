# Intel compiler (e.g., yana)
FC=ifort
FLAGS= -O3 -llapack
LINK=$(FC)

# gfortran compiler (e.g., yana)
#FC=gfortran
#CC=gcc
#CFLAGS = -Wall
#FLAGS= -w -O3 -ffixed-line-length-none
#LINK=$(FC)

# Pathscale compiler (e.g., jacquard)
#FC=pathf90
#FLAGS= -w -O3 -llapack
#LINK=$(FC)


# ------------ BASIC PARAMS ------------
SRC=src
OBJ=obj
COMPTARGET=wan2pov
OBJECTS= \
	$(OBJ)/main.o
GENMODS= \
	$(OBJ)/general.o \
	$(OBJ)/io.o 


MODULES=$(OBJ)/write_df3.o\
	$(OBJ)/driver.o
#MODULES= \
#	$(OBJ)/wannier.o \
#	$(OBJ)/twowannier.o \
#	$(OBJ)/driver.o


$(COMPTARGET): $(MODULES) $(OBJECTS)
	@echo ".linking object files to final binary"
	@rm -f $(COMPTARGET)
	@cd $(OBJ) ; $(LINK) $(FLAGS) -o $(COMPTARGET) `ls *.o`
	@mv $(OBJ)/$(COMPTARGET) .

 

# ------------ MODULES ------------
#  nested string substitution to turn 'f90' into 'o' and 'obj' into 'src'
$(OBJ)/general.o: $(SRC)/general.f90
	@echo ".compiling file $(subst $(OBJ),$(SRC),$(subst .o,.f90,$@))"
	@rm -f $@
	@cp $(subst $(OBJ),$(SRC),$(subst .o,.f90,$@)) $(OBJ)/
	@cd $(OBJ) ; $(FC) $(FLAGS) -c $(subst $(OBJ)/,./,$(subst .o,.f90,$@))

$(OBJ)/io.o: $(SRC)/io.f90
	@echo ".compiling file $(subst $(OBJ),$(SRC),$(subst .o,.f90,$@))"
	@rm -f $@
	@cp $(subst $(OBJ),$(SRC),$(subst .o,.f90,$@)) $(OBJ)/
	@cd $(OBJ) ; $(FC) $(FLAGS) -c $(subst $(OBJ)/,./,$(subst .o,.f90,$@))

$(OBJ)/write_df3.o: $(GENMODS) $(SRC)/write_df3.c
	@echo ".compiling file $(subst $(OBJ),$(SRC),$(subst .o,.c,$@))"
	@rm -f $@
	@cp $(subst $(OBJ),$(SRC),$(subst .o,.c,$@)) $(OBJ)/
	@cd $(OBJ) ; $(CC) $(CFLAGS) -c $(subst $(OBJ)/,./,$(subst .o,.c,$@))

$(OBJ)/driver.o: $(GENMODS) $(SRC)/driver.f90
	@echo ".compiling file $(subst $(OBJ),$(SRC),$(subst .o,.f90,$@))"
	@rm -f $@
	@cp $(subst $(OBJ),$(SRC),$(subst .o,.f90,$@)) $(OBJ)/
	@cd $(OBJ) ; $(FC) $(FLAGS) -c $(subst $(OBJ)/,./,$(subst .o,.f90,$@))


# ------------ OBJECTS ------------
$(OBJ)/main.o: $(GENMODS) $(MODULES) $(SRC)/main.f90
	@echo ".compiling file $(subst $(OBJ),$(SRC),$(subst .o,.f90,$@))"
	@rm -f $@
	@cp $(subst $(OBJ),$(SRC),$(subst .o,.f90,$@)) $(OBJ)/
	@cd $(OBJ) ; $(FC) $(FLAGS) -c $(subst $(OBJ)/,./,$(subst .o,.f90,$@))




clean:
	@rm -fr obj $(COMPTARGET)
	@mkdir $(OBJ)
