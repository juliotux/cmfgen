# Include local system definitions

include ../../../Makefile_definitions

# Library where object files will/are stored

LIB=$(LIB_DIR)libdev.a
LIB_MOD=$(LIB_DIR)libmod.a

#
# Entries for library
#
OBJ = $(LIB)(steq_ba_two_phot_rate_v3.o)\
      $(LIB)(get_j_for_two_phot.o)\
      $(LIB)(two_phot_rate.o)\
      $(LIB)(set_two_phot_v2.o)\
      $(LIB)(set_two_phot_v3.o)\
      $(LIB)(two_phot_opac.o)\
      $(LIB)(two_phot_opac_v3.o)\
      $(LIB)(two_phot_var_opac.o)

all : $(LIB)
	rm -f *.o
$(LIB) : $(OBJ)	

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f $(TWO)
	$(F90) -c  $(FG) $< 
	ar ruv $(LIB) $*.o

# Specific dependency of source files which uses external 
# fortran moudles. (The order of compilation matters.)


TWO = $(LIB_MOD)(two_phot_mod.o)

ST = $(LIB_MOD)(steq_data_mod.o)

$(LIB)(steq_ba_two_phot_rate_v3.o) : steq_ba_two_phot_rate_v3.f $(TWO) $(ST) 
	$(F90) -c  $(FG) steq_ba_two_phot_rate_v3.f
	ar ruv $(LIB) steq_ba_two_phot_rate_v3.o

$(TWO):
	(cd ../../mod_subs ; make)

$(ST):
	(cd ../../mod_subs ; make)

clean :
	rm -f *.o 

clean_lib :
	rm $(LIB)
