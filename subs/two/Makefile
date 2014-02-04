# Include local system definitions

include ../../Makefile_definitions

# Library where object files will/are stored

LIB=$(LIB_DIR)libsubs.a
LIB_MOD=$(LIB_DIR)libmod.a

#
# Entries for library
#
OBJ = $(LIB)(steq_ba_two_phot.o)\
      $(LIB)(two_phot_rate.o)\
      $(LIB)(set_two_phot.o)\
      $(LIB)(two_phot_opac.o)\
      $(LIB)(two_phot_var_opac.o)

all : $(LIB)
	cp $(INSTALL_DIR)perm_mod_file jnk.mod
	mv -f *.mod $(MOD_DIR)
	rm -f *.o
$(LIB) : $(OBJ)	

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f $(LIB_MOD)(two_phot_mod.o) 
	$(F90) -c $(FG)  $< 
	ar ruv $(LIB) $*.o

# Specific dependency of source files which uses external 
# fortran moudles. (The order of compilation matters.)


$(LIB_MOD)(two_phot_mod.o):
	(cd ../new_main/mod_subs; make)

clean :
	rm -f *.o 

clean_lib :
	rm $(LIB)
