#
include ../Makefile_definitions

LIB = $(LIB_DIR)libmy_lapack.a

OBJ =  dgeequ.o  dgetrf.o  dlamch.o  ieeeck.o  lsame.o\
       dgetf2.o  dgetrs.o  dlaswp.o  ilaenv.o

all : $(LIB)
	rm -f *.o

$(LIB) : $(LIB)($(OBJ))

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f
	$(F90) -c $(FG)  $<
	ar ruv $(LIB) $*.o

#
# Specific dependency of source files which uses external
# fortran moudles. (The order of compilation matters.)


clean :
	rm -f *.o

clean_lib :
	rm -f $(NIB)
