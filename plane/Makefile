#
# Include local system definitions

include ../Makefile_definitions

OPTION=$(FG)

# Library where object files will/are stored

LIB=$(LIB_DIR)libsubs.a
LIB_MOD=$(LIB_DIR)libmod.a
#
# Entries for library
#

MOD=  $(LIB)(mod_ray_mom_store.o)\
      $(LIB)(mod_var_mom_j_cmf.o)\
      $(LIB)(mod_var_jrel_v2.o)\
      $(LIB)(mod_space_grid_v2.o)

OBJ = $(LIB)(comp_add_edd_facs.o)\
      $(LIB)(characteristics_v2.o)\
      $(LIB)(cmf_formal_rel_v4.o)\
      $(LIB)(derivr.o)\
      $(LIB)(derivs.o)\
      $(LIB)(drhsdchi_pp.o)\
      $(LIB)(define_grid_v2.o)\
      $(LIB)(expn.o)\
      $(LIB)(edd_j_hub_var_v2.o)\
      $(LIB)(edd_jrel_var_v3.o)\
      $(LIB)(fcomp_pp_v2.o)\
      $(LIB)(fgrey_norel_v1.o)\
      $(LIB)(fg_j_cmf_v12.o)\
      $(LIB)(fg_j_cmf_v13.o)\
      $(LIB)(gauleg.o)\
      $(LIB)(get_ibound.o)\
      $(LIB)(get_jh_at_prev_time_step.o)\
      $(LIB)(get_moms_non_rel.o)\
      $(LIB)(get_moms_rel.o)\
      $(LIB)(get_rsqh_rel.o)\
      $(LIB)(jgrey_hub_ddt_v2.o)\
      $(LIB)(jgrey_hub_ddt_v4.o)\
      $(LIB)(j_weight.o)\
      $(LIB)(h_weight.o)\
      $(LIB)(k_weight.o)\
      $(LIB)(n_weight.o)\
      $(LIB)(mom_j_cmf_v8.o)\
      $(LIB)(mom_j_cmf_v11.o)\
      $(LIB)(mom_j_ddt_v4.o)\
      $(LIB)(mom_j_pp_v1.o)\
      $(LIB)(mom_jrel_v8.o)\
      $(LIB)(mom_jrel_grey_v1.o)\
      $(LIB)(mom_jrel_grey_v2.o)\
      $(LIB)(optdepth_v4.o)\
      $(LIB)(out_jh.o)\
      $(LIB)(pp_form_cmf_v2.o)\
      $(LIB)(pp_var_mom_cmf_v1.o)\
      $(LIB)(pp_edd_var_cmf_v1.o)\
      $(LIB)(pp_mom_cmf_v1.o)\
      $(LIB)(rel_variables.o)\
      $(LIB)(runge_kutta.o)\
      $(LIB)(solve_cmf_formal_v2.o)\
      $(LIB)(var_jrel_v4.o)\
      $(LIB)(var_mom_j_ddt_v4.o)\
      $(LIB)(var_mom_pp_v1.o)\
      $(LIB)(var_mom_j_cmf_v9.o)\
      $(LIB)(var_mom_j_cmf_v12.o)\
      $(LIB)(var_mom_jext_cmf_v3.o)\
      $(LIB)(velocity_law.o)

all : $(LIB)
	cp -f $(INSTALL_DIR)perm_mod_file jnk.mod
	mv -f *.mod $(MOD_DIR)
	rm -f *.o

$(LIB) : $(MOD) $(OBJ)	

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f 
	$(F90) -c $(OPTION) $<
	ar ruv $(LIB) $*.o

%.o : %.f 
	$(F90) -c $(OPTION) $<
	ar ruv $(LIB) $*.o

$(LIB)(characteristics_v2.o) : characteristics_v2.f $(LIB)(mod_space_grid_v2.o)
	$(F90) -c $(OPTION) characteristics_v2.f
	ar ruv $(LIB) characteristics_v2.o

$(LIB)(cmf_formal_rel_v2.o) : cmf_formal_rel_v2.f $(LIB)(mod_space_grid_v2.o) $(LIB)(mod_ray_mom_store.o)
	$(F90) -c $(OPTION) cmf_formal_rel_v2.f
	ar ruv $(LIB) cmf_formal_rel_v2.o

$(LIB)(define_grid_v2.o) : define_grid_v2.f $(LIB)(mod_space_grid_v2.o)
	$(F90) -c $(OPTION) define_grid_v2.f
	ar ruv $(LIB) define_grid_v2.o

$(LIB)(drhsdchi_pp.o) : drhsdchi_pp.f $(LIB)(mod_trap_derivatives.o)
	$(F90) -c $(OPTION) drhsdchi_pp.f
	ar ruv $(LIB) drhsdchi_pp.o

$(LIB)(edd_jrel_var_v2.o) : edd_jrel_var_v2.f $(LIB)(mod_var_jrel_v2.o)
	$(F90) -c $(OPTION) edd_jrel_var_v2.f
	ar ruv $(LIB) edd_jrel_var_v2.o

$(LIB)(fg_j_cmf_v12.o) : fg_j_cmf_v12.f $(LIB)(mod_ray_mom_store.o)
	$(F90) -c $(OPTION) fg_j_cmf_v12.f
	ar ruv $(LIB) fg_j_cmf_v12.o

$(LIB)(get_moms_rel.o) : get_moms_rel.f $(LIB)(mod_ray_mom_store.o)
	$(F90) -c $(OPTION) get_moms_rel.f
	ar ruv $(LIB) get_moms_rel.o

$(LIB)(get_moms_non_rel.o) : get_moms_non_rel.f $(LIB)(mod_ray_mom_store.o)
	$(F90) -c $(OPTION) get_moms_non_rel.f
	ar ruv $(LIB) get_moms_non_rel.o

$(LIB)(mom_j_cmf_v8.o) : mom_j_cmf_v8.f $(LIB)(mod_ray_mom_store.o)
	$(F90) -c $(OPTION) mom_j_cmf_v8.f
	ar ruv $(LIB) mom_j_cmf_v8.o

$(LIB)(mom_j_ddt_v2.o) : mom_j_ddt_v2.f $(LIB)(mod_ray_mom_store.o)
	$(F90) -c $(OPTION) mom_j_ddt_v2.f
	ar ruv $(LIB) mom_j_ddt_v2.o

$(LIB)(mom_jrel_v7.o) : mom_jrel_v7.f $(LIB)(mod_ray_mom_store.o)
	$(F90) -c $(OPTION) mom_jrel_v7.f
	ar ruv $(LIB) mom_jrel_v7.o

$(LIB)(optdepth_v2.o) : optdepth_v2.f $(LIB)(mod_space_grid_v2.o)
	$(F90) -c $(OPTION) optdepth_v2.f
	ar ruv $(LIB) optdepth_v2.o

$(LIB)(solve_cmf_formal_v2.o) : solve_cmf_formal_v2.f $(LIB)(mod_space_grid_v2.o)
	$(F90) -c $(OPTION) solve_cmf_formal_v2.f
	ar ruv $(LIB) solve_cmf_formal_v2.o

$(LIB)(var_mom_j_cmf_v9.o) : var_mom_j_cmf_v9.f $(LIB)(mod_ray_mom_store.o) $(LIB)(mod_var_mom_j_cmf.o)
	$(F90) -c $(OPTION) var_mom_j_cmf_v9.f
	ar ruv $(LIB) var_mom_j_cmf_v9.o

$(LIB)(var_mom_j_ddt_v2.o) : var_mom_j_ddt_v2.f $(LIB)(mod_ray_mom_store.o)
	$(F90) -c $(OPTION) var_mom_j_ddt_v2.f
	ar ruv $(LIB) var_mom_j_ddt_v2.o

$(LIB)(var_jrel_v3.o) : var_jrel_v3.f $(LIB)(mod_var_jrel_v2.o) $(LIB)(mod_ray_mom_store.o) $(LIB)(mod_trap_derivatives.o) 
	$(F90) -c $(OPTION) var_jrel_v3.f
	ar ruv $(LIB) var_jrel_v3.o

$(LIB)(var_mom_jext_cmf_v3.o) : var_mom_jext_cmf_v3.f $(LIB)(mod_ray_mom_store.o)
	$(F90) -c $(OPTION) var_mom_jext_cmf_v3.f
	ar ruv $(LIB) var_mom_jext_cmf_v3.o

$(LIB)(mod_trap_derivatives.o) :
	(cd ../subs ; make)

clean :
	rm -f *.o 

clean_lib :
	rm -f *.a
