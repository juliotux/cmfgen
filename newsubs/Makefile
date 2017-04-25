# Include local system definitions

include ../Makefile_definitions

# Library where object files will/are stored
         
LIB=$(LIB_DIR)libnewsubs.a
LIB_MOD=$(LIB_DIR)libmod.a

#
# Entries in a library
#
OBJ= $(LIB)(ave_level_energy.o)\
     $(LIB)(colcool_sl_v4.o)\
     $(LIB)(colcool_sl_v5.o)\
     $(LIB)(cut_points_v3.o)\
     $(LIB)(do_hmi_ff.o)\
     $(LIB)(ff_res_gaunt.o)\
     $(LIB)(genopaeta_v8.o)\
     $(LIB)(genopaeta_v10.o)\
     $(LIB)(genosc_v8.o)\
     $(LIB)(genosc_v9.o)\
     $(LIB)(get_scl_xray_fluxes_v1.o)\
     $(LIB)(newpertjd.o)\
     $(LIB)(omega_gen_v3.o)\
     $(LIB)(phot_data_mod.o)\
     $(LIB)(quad_multi_v9.o)\
     $(LIB)(rayleigh_scat.o)\
     $(LIB)(rdphot_gen_v2.o)\
     $(LIB)(rd_hyd_bf_data.o)\
     $(LIB)(rd_phot_die_v1.o)\
     $(LIB)(rd_rv_file_v2.o)\
     $(LIB)(rd_mod_density.o)\
     $(LIB)(rd_store_log.o)\
     $(LIB)(rd_xray_spec.f)\
     $(LIB)(set_edge_freq_v3.o)\
     $(LIB)(set_edge_freq_v4.o)\
     $(LIB)(set_x_freq_v2.o)\
     $(LIB)(sm_phot_v3.o)\
     $(LIB)(subcol_multi_v4.o)\
     $(LIB)(subcol_multi_v6.o)\
     $(LIB)(sub_phot_gen.o)\
     $(LIB)(xcross_v2.o)\
     $(LIB)(xray_data_mod.o)\
     $(LIB)(var_op_v6.o)

all : $(LIB)
	cp -f $(INSTALL_DIR)perm_mod_file jnk.mod
	mv -f *.mod $(MOD_DIR)
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
# NB: change_folor.f, change_pen.f, gramon_pgplot.f, vectorpen.f
#     depend on ../tools/gen_in

$(LIB)(ff_res_gaunt.o)    : ff_res_gaunt.f    $(LIB)(phot_data_mod.o)
$(LIB)(rdphot_gen_v2.o)   : rdphot_gen_v2.f   $(LIB)(phot_data_mod.o)

$(LIB)(sub_phot_gen.o)    : sub_phot_gen.f    $(LIB)(phot_data_mod.o) $(LIB)(rd_hyd_bf_data.o) $(LIB)(xray_data_mod.o)

GEN = $(LIB_MOD)(mod_lev_dis_blk.o)

$(LIB)(genopaeta_v10.o)  : genopaeta_v10.f    $(GEN)
$(LIB)(genopaeta_v8.o)   : genopaeta_v8.f     $(GEN)
$(LIB)(quad_multi_v8.o)  : quad_multi_v8.f    $(GEN)

$(LIB)(xcross_v2.o)     : xcross_v2.f       $(LIB)(xray_data_mod.o)
$(LIB)(set_x_freq_v2.o) : set_x_freq_v2.f   $(LIB)(xray_data_mod.o)

$(LIB)(get_scl_xray_fluxes_v1.o)     : get_scl_xray_fluxes_v1.f       $(LIB)(rd_xray_spec.o)

$(GEN):
	(cd ../subs; make )

clean :
	rm -f *.o

clean_lib :
	rm -f $(LIB)
