# Include local system definitions

include ../../Makefile_definitions

OPTION=$(FG)

# Library where object files will/are stored

LIB=$(LIB_DIR)libdev.a
LIB_MOD=$(LIB_DIR)libmod.a

# Entries for library

OBJ = $(LIB)(add_par_to_full_v2.o)\
      $(LIB)(adjust_cmat_to_lte.o)\
      $(LIB)(adjust_sn_r_grid.o)\
      $(LIB)(adjust_pops.o)\
      $(LIB)(auto_clump_rev.o)\
      $(LIB)(ba_update_v7.o)\
      $(LIB)(check_hydro_def.o)\
      $(LIB)(check_ions_present.o)\
      $(LIB)(check_tmin.o)\
      $(LIB)(cmf_blkband_v3.o)\
      $(LIB)(create_iv_links_v2.o)\
      $(LIB)(determine_nse.o)\
      $(LIB)(det_lte_ed.o)\
      $(LIB)(ddt_work.o)\
      $(LIB)(do_species_decays_v2.o)\
      $(LIB)(do_full_r_grid_v1.o)\
      $(LIB)(do_vel_regrid_v2.o)\
      $(LIB)(do_tau_regrid_v2.o)\
      $(LIB)(do_gam_abs_approx_v2.o)\
      $(LIB)(do_cmf_hydro_v2.o)\
      $(LIB)(do_ng_accel_band_v2.o)\
      $(LIB)(eval_adiabatic_v3.o)\
      $(LIB)(evalse_qwvj_v7.o)\
      $(LIB)(evalse_lowt_v1.o)\
      $(LIB)(evalse_x_qwvj_v4.o)\
      $(LIB)(eval_temp_ddt_v2.o)\
      $(LIB)(fiddle_pop_corrections_v2.o)\
      $(LIB)(fixpop_in_ba_v2.o)\
      $(LIB)(fixpop_t.o)\
      $(LIB)(fixpop_in_ba_v3.o)\
      $(LIB)(generate_full_matrix_v3.o)\
      $(LIB)(get_decay_lum_v1.o)\
      $(LIB)(get_non_loc_gamma_energy_v2.o)\
      $(LIB)(get_pops_at_prev_time_step_v5.o)\
      $(LIB)(old_grid_module.o)\
      $(LIB)(prevent_low_t.o)\
      $(LIB)(read_ba_data_v3.o)\
      $(LIB)(read_bcd_mat.o)\
      $(LIB)(read_time_model_v1.o)\
      $(LIB)(read_seq_time_file_v1.o)\
      $(LIB)(rv_sn_model_sniin.o)\
      $(LIB)(set_ba_storage.o)\
      $(LIB)(set_depth_consistency.o)\
      $(LIB)(set_rv_hydro_model_v3.o)\
      $(LIB)(set_imp_vec.o)\
      $(LIB)(set_pops_to_pure_lte.o)\
      $(LIB)(set_sn_r_grid.o)\
      $(LIB)(set_new_grid.o)\
      $(LIB)(solveba_v13.o)\
      $(LIB)(steq_multi_v9.o)\
      $(LIB)(steqne_v4.o)\
      $(LIB)(steq_advec_v4.o)\
      $(LIB)(steq_co_mov_deriv_v3.o)\
      $(LIB)(store_old_grid.o)\
      $(LIB)(store_ba_data_v3.o)\
      $(LIB)(sum_steq_sol.o)\
      $(LIB)(tst_rd_seq.o)\
      $(LIB)(var_op_v10.o)\
      $(LIB)(var_x_opa_eta_v4.o)\
      $(LIB)(vsebyj_multi_v8.o)\
      $(LIB)(vsebyj_x_v7.o)\
      $(LIB)(write_cmat.o)\
      $(LIB)(write_seq_time_file_v1.o)\
      $(LIB)(update_ba_for_line.o)\
      $(LIB)(wr_asci_steq.o)

all : $(LIB) XNUC XAUTO XCHG XTWO XNON
	cp -f $(INSTALL_DIR)perm_mod_file jnk.mod
	mv -f *.mod $(MOD_DIR)
	rm -f *.o

$(LIB) : $(OBJ)

STEQ=$(LIB_MOD)(steq_data_mod.o)
CHG=$(LIB_MOD)(chg_exch_mod_v3.o)
MOD=$(LIB_MOD)(mod_cmfgen.o)
DIS=$(LIB_MOD)(mod_lev_dis_blk.o)
CONT=$(LIB_MOD)(control_variable_mod.o)
NUC=$(LIB_MOD)(nuc_iso_mod.o)
HYD=$(LIB_MOD)(cmf_hydro_module.o)

UP_MODS= $(LIB_DIR)libmod.a(mod_cmfgen.o)\
       $(LIB_DIR)libmod.a(control_variable_mod.o)\
       $(LIB_DIR)libmod.a(line_mod.o)\
       $(LIB_DIR)libmod.a(line_vec_mod.o)\
       $(LIB_DIR)libmod.a(opac_mod.o)\
       $(LIB_DIR)libmod.a(radiation_mod.o)\
       $(LIB_DIR)libmod.a(steq_data_mod.o)\
       $(LIB_DIR)libmod.a(var_rad_mod.o)


# Note: Not all files depend on mod_cmfgen and steq_data_mod, but assuming this dependency
# reduces signicantly the number of specific specifications.

$(LIB)(%.o) : %.f  $(STEQ) $(MOD) $(CONT)
	$(F90) -c $(OPTION) $<
	ar ruv $(LIB) $*.o

# General rules to make obects in a library and object files
# Specific dependency of source files which uses external 
# fortran moudles. (The order of compilation matters.)

$(LIB)(adjust_pops.o) : adjust_pops.f $(LIB)(old_grid_module.o) $(MOD) $(CONT)
	$(F90) -c $(OPTION) adjust_pops.f
	ar ruv $(LIB) adjust_pops.o

$(LIB)(determine_nse.o) : determine_nse.f $(STEQ) $(CHG) $(MOD)
	$(F90) -c $(OPTION) determine_nse.f
	ar ruv $(LIB) determine_nse.o

$(LIB)(do_cmf_hydro_v2.o) : do_cmf_hydro_v2.f $(LIB)(old_grid_module.o) $(LIB)(cmf_hydro_module.o)
	$(F90) -c $(OPTION) do_cmf_hydro_v2.f
	ar ruv $(LIB) do_cmf_hydro_v2.o

# Optimizaton seems not to work with PGF95 compiler for next routine

$(LIB)(eval_temp_ddt_v2.o) : eval_temp_ddt_v2.f $(STEQ) $(MOD)
	$(F90) -c $(OPTION)  eval_temp_ddt_v2.f
	ar ruv $(LIB) eval_temp_ddt_v2.o

$(LIB)(get_non_loc_gamma_energy.o) : get_non_loc_gamma_energy.f $(CONT) $(NUC)
	$(F90) -c $(OPTION) get_non_loc_gamma_energy.f
	ar ruv $(LIB) get_non_loc_gamma_energy.o

$(LIB)(store_old_grid.o) : store_old_grid.f $(MOD)  $(LIB)(old_grid_module.o)
	$(F90) -c $(OPTION) store_old_grid.f
	ar ruv $(LIB) store_old_grid.o

$(LIB)(update_ba_for_line.o) : update_ba_for_line.f $(UP_MODS)
	$(F90) -c $(OPTION)  update_ba_for_line.f
	ar ruv $(LIB) update_ba_for_line.o

$(LIB)(var_op_v10.o) : var_op_v10.f $(DIS)
	$(F90) -c $(OPTION)  var_op_v10.f
	ar ruv $(LIB) var_op_v10.o

$(STEQ):
	(cd ../mod_subs ; make)

$(MOD):
	(cd ../mod_subs ; make)

$(NUC):
	(cd ../mod_subs ; make)

$(UP_MOD):
	(cd ../mod_subs ; make)

$(DIS):
	(cd ../mod_subs ; make)

$(CHG):
	(cd chg ; make)

XNUC:
	(make -f Make_nuc)

XAUTO:
	(cd auto; make)
XCHG:
	(cd chg; make)
XTWO:
	(cd two; make)
XNON:
	(cd non_therm; make)

clean :
	rm -f *.o 

clean_lib :
	rm $(LIB)
