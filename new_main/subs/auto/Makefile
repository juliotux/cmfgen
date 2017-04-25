# Include local system definitions

include ../../../Makefile_definitions

# Library where object files will/are stored

LIB=$(LIB_DIR)libdev.a
LIB_MOD=$(LIB_DIR)libmod.a

#
# Entries for library
#
OBJ = $(LIB)(rd_auto_v1.o)\
      $(LIB)(add_auto_rates.o)\
      $(LIB)(steq_auto_v2.o)

all : $(LIB)
	rm -f *.o
$(LIB) : $(OBJ)	

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f 
	$(F90) -c $(FG)  $< 
	ar ruv $(LIB) $*.o

# Specific dependency of source files which uses external 
# fortran moudles. (The order of compilation matters.)

ST = $(LIB_MOD)(steq_data_mod.o)

$(LIB)(steq_auto_v2.o) : steq_auto_v2.f  $(ST)
	$(F90) -c $(FG) steq_auto_v2.f
	ar ruv $(LIB) steq_auto_v2.o

$(ST):
	(cd ../mod_subs ; make)

clean :
	rm -f *.o 

clean_lib :
	rm $(LIB)
