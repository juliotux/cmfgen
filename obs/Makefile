#Include file with system dependent definitions

include ../Makefile_definitions

OPTION=$(FG)

#
# These names are for the fortran linker
#

LN1=tools
LN2=disp
LN3=subs
LN4=newsubs
LN5=pgplt
LN6=mod
LN7=dev

#
# These names refer to the full file name
#

FN1=$(LIB_DIR)lib$(LN1).a
FN2=$(LIB_DIR)lib$(LN2).a
FN3=$(LIB_DIR)lib$(LN3).a
FN4=$(LIB_DIR)lib$(LN4).a
FN5=$(LIB_DIR)lib$(LN5).a
FN6=$(LIB_DIR)lib$(LN6).a
FN7=$(LIB_DIR)lib$(LN7).a

#
# We only make the executables if there has been a change to the main *.f files,
# or if the libraries have changed since the last link.
#
 
OBJ1  =  cmf_flux_v5.o cmf_flux_sub_v5.o mod_cmf_obs.o mod_freq_obs.o \
        cmf_flux_cntrl_var_mod.o rd_cmf_flux_controls.o check_cmf_flux_param_consis.o \
        rd_model_file.o revise_obs_p.o  ins_line_v5.o obs_frame_sub_v9.o\
        ins_line_obs_v4.o rd_cont_j_obs.o set_forbid_zero.o
 
OBJ2  =  obs_frame.o obs_frame_sub_v2.o 
 
all: $(EXE_DIR)cmf_flux.exe $(EXE_DIR)obs_frame.exe

#
# Give the dependencies of cmf_flux.exe and how it is created.
# 

$(EXE_DIR)cmf_flux.exe :  $(OBJ1) $(FN1) $(FN2) $(FN3) $(FN4) $(FN5) $(FN6) $(FN7)
	ar -xo ../lib/libsubs.a griem_stark_mod.o
	$(F90) $(OBJ1) griem_stark_mod.o \
              -L$(LIB_DIR) -l$(LN7) -l$(LN4) -l$(LN2) -l$(LN3) -l$(LN5) -l$(LN1) -l$(LN6) \
               $(LOCLIB) \
               -o $(EXE_DIR)cmf_flux.exe

$(EXE_DIR)obs_frame.exe :  $(OBJ2) $(FN1) $(FN2) $(FN3) $(FN4)
	$(F90) $(OBJ2) \
              -L$(LIB_DIR) -l$(LN4) -l$(LN2) -l$(LN3) -l$(LN1) \
               $(LOCLIB) \
              -o  $(EXE_DIR)obs_frame.exe

# General rule to compile a fortran ".f" file.

%.o : %.f
	$(F90) -c $(OPTION)  $<

#
# I use F90noamp because of an issue with the PGF95 compiler on the Mac PRO.
#

cmf_flux_v5.o : cmf_flux_v5.f  mod_cmf_obs.o  cmf_flux_cntrl_var_mod.o $(USR)
	$(F90) -c $(OPTION) cmf_flux_v5.f 

#	$(F90nomp) -c $(OPTION) cmf_flux_v5.f 



cmf_flux_sub_v5.o : cmf_flux_sub_v5.f  mod_cmf_obs.o  mod_freq_obs.o cmf_flux_cntrl_var_mod.o $(GEN) \
   COMP_JCONT_V4.INC  EVAL_LTE_INC_V5.INC\
   OPACITIES_V5.INC   XOPAC_V4.INC
	$(F90) -c $(OPTION) cmf_flux_sub_v5.f 

#	$(F90) -c $(OPTION) cmf_flux_sub_v5.f 




USR = ../lib/libtools.a(usr_option.o)
GEN = ../lib/libmod.a(mod_lev_dis_blk.o)

rd_cmf_flux_controls.o : cmf_flux_cntrl_var_mod.o
check_cmf_flux_param_consis.o : cmf_flux_cntrl_var_mod.o

$(GEN):
	(cd ../subs; make )
$(USR):
	(cd ../tools; make )

 
$(FN1) : lib1
$(FN2) : lib2
$(FN3) : lib3
$(FN4) : lib4
$(FN5) : lib5
$(FN6) : lib6
$(FN7) : lib7

# We now MAKE the required libraries.

lib1:
	(cd ../tools; make)
lib2:
	(cd ../disp/subs; make)
lib3:
	(cd ../subs; make )
	(cd ../unix; make )
	(cd ../stark; make )
lib4:
	(cd ../newsubs; make)
lib5:
	(cd ../pgplt; make)
lib6:
	(cd ../new_main/mod_subs; make)
lib7:
	(cd ../new_main/subs; make)

clean :
	rm -f *.o
