# TEST MAKEFILE

#Include file with system dependent definitions

include ../Makefile_definitions

OPTIONS=$(FG)


# These names are for the fortran linker

LN1=tools
LN2=pgplt
LN3=subs
LN4=newsubs
LN5=dev
LN6=mod

# These names refer to the full file name

FN1=$(LIB_DIR)lib$(LN1).a
FN2=$(LIB_DIR)lib$(LN2).a
FN3=$(LIB_DIR)lib$(LN3).a
FN4=$(LIB_DIR)lib$(LN4).a
FN5=$(LIB_DIR)lib$(LN5).a
FN6=$(LIB_DIR)lib$(LN6).a

# We only make dispgen.exe if there has been a change to the main *.f files,
# or if the libraries have changed since the last link.
  
all : $(EXE_DIR)cmfgen_dev.exe

UPDATE=$(LIB_DIR)libtools.a(update_keyword.o)

# Give the dependencies of dispgen.exe and how it is created.

$(EXE_DIR)cmfgen_dev.exe : cmfgen.o cmfgen_sub.o \
                                  $(FN1) $(FN2) $(FN3) $(FN4) $(FN5) $(FN6)
	$(F90) cmfgen.o cmfgen_sub.o \
              -L$(LIB_DIR) -l$(LN6) -l$(LN5) -l$(LN4) -l$(LN3) -l$(LN2) -l$(LN5) -l$(LN3) -l$(LN1) -l$(LN6) \
               $(LOCLIB) $(LAPACK) $(BLAS) \
              -o $(EXE_DIR)cmfgen_dev.exe

USR=$(LIB_DIR)libtools.a(usr_option.o)

MODS= $(LIB_DIR)libmod.a(mod_cmfgen.o)\
      $(LIB_DIR)libmod.a(cmf_sob_mod.o)\
      $(LIB_DIR)libmod.a(control_variable_mod.o)\
      $(LIB_DIR)libmod.a(ang_qw_mod.o)\
      $(LIB_DIR)libmod.a(line_mod.o)\
      $(LIB_DIR)libmod.a(line_vec_mod.o)\
      $(LIB_DIR)libmod.a(mod_lev_dis_blk.o)\
      $(LIB_DIR)libmod.a(opac_mod.o)\
      $(LIB_DIR)libmod.a(radiation_mod.o)\
      $(LIB_DIR)libmod.a(steq_data_mod.o)\
      $(LIB_DIR)libmod.a(var_rad_mod.o)

cmfgen.o : cmfgen.f $(MODS) $(USR)
	$(F90) -c $(OPTIONS) cmfgen.f 

#cmfgen.o : cmfgen.f $(MODS)
#	$(F90) -c $(OPTIONS) cmfgen.f 

cmfgen_sub.o : cmfgen_sub.f $(MODS) $(UPDATE) GET_J_CHI.INC  
	$(F90) -c $(OPTIONS) cmfgen_sub.f 

$(USR):
	(cd ../tools ; make )

$(UPDATE):
	(cd ../tools ; make)

$(MODS): mod

mod:
	(cd mod_subs ; make )

# We now MAKE the required libraries. Note that this double
# dependence is rquired to ensure that the libraries are updated.
# Otherwise Makfile simple checks date of library.

$(FN1): lib1
$(FN2): lib2
$(FN3): lib3
$(FN4): lib4
$(FN5): lib5
$(FN6): lib6

lib1:
	(cd ../tools; make )
lib2:
	(cd ../pgplt; make )
lib3:
	(cd ../subs; make )
	(cd ../stark; make )
lib4:
	(cd ../newsubs; make )
lib5:
	(cd subs; make )
lib6:
	(cd mod_subs; make )


clean :
	rm -f *.o

clean_lib :
	rm -f $(LIB)
