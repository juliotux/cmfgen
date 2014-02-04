# Include local system definitions

include ../../Makefile_definitions

# Library where object files will/are stored

LIB=$(LIB_DIR)libsubs.a

#
# Entries in a library
#
OBJ = $(LIB)(chg_exch_mod.o)\
      $(LIB)(set_chg_exch_v2.o)\
      $(LIB)(verify_chg_exch.o)\
      $(LIB)(eval_chg_rates.o)\
      $(LIB)(steq_ba_chg_exch.o)\
      $(LIB)(wr_chg_cool.o)

all : $(LIB)
	cp -f $(INSTALL_DIR)perm_mod_file jnk.mod
	mv -f *.mod $(MOD_DIR)
	rm -f *.o
$(LIB) : $(OBJ)	

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f $(LIB)(chg_exch_mod.o) 
	$(F90) -c $(FG)  $< 
	ar ruv $(LIB) $*.o

# Specific dependency of source files which uses external 
# fortran moudles. (The order of compilation matters.)

$(LIB)(chg_exch_mod.o) : chg_exch_mod.f 
	$(F90) -c $(FG) chg_exch_mod.f 
	ar ruv $(LIB) chg_exch_mod.o

#$(LIB)(.o) : usr_option.f $(LIB)(wr_string.o) 

clean :
	rm -f *.o 

clean_lib :
	rm $(LIB)
