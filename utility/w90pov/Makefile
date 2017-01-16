include ../../make.inc

SRC=src
OBJ=obj
COMPTARGET=w90pov
OBJECTS= \
	$(OBJ)/main.o
GENMODS= \
	$(OBJ)/general.o \
	$(OBJ)/io.o 

MODULES=$(OBJ)/write_df3.o\
	$(OBJ)/driver.o

$(COMPTARGET): $(MODULES) $(OBJECTS) $(GENMODS)
	$(F90) $(FCOPTS) -o $(COMPTARGET) $^

$(OBJ)/write_df3.o: $(GENMODS)
$(OBJ)/driver.o: $(GENMODS)
$(OBJ)/main.o: $(GENMODS) $(MODULES)
$(OBJ)/io.o: $(OBJ)/general.o

$(OBJ)/%.o : $(SRC)/%.f90
	$(F90) $(FCOPTS) -c $< -o $@

$(OBJ)/%.o : $(SRC)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -fr $(COMPTARGET) $(OBJ)/*.o *.mod
