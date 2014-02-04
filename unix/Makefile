#
#Include system definitions

include ../Makefile_definitions


# Library in which object files will be stored.

LIB=$(LIB_DIR)libsubs.a

#
# Entries in a library
#
OBJ= $(LIB)(dir_acc_pars_gen.o)\
     $(LIB)(date_time.o)\
     $(LIB)(verbose_info.o)\
     $(LIB)(wr_pwd.o)\
     $(LIB)(set_line_buffering.o)\
     $(LIB)(tune.o)\
     $(LIB)(gen_asci_open.o)\
     $(LIB)(f01btf.o)\
     $(LIB)(f04ayf.o)\
     $(LIB)(term_io.o)

all : $(LIB)
	rm -f *.o
$(LIB) : $(OBJ)	

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f 
	$(F90) -c $(FG)  $<
	ar ruv $(LIB) $*.o

# 
# Specific dependency of source files which uses external 
# fortran moudles. (The order of compilation matters.)
#

clean :
	rm -f *.o

clean_lib :
	rm -f $(LIB)
