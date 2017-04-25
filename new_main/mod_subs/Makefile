# Include local system definitions

include ../../Makefile_definitions

# Library where object files will/are stored

LIB=$(LIB_DIR)libdev.a
LIB_MOD=$(LIB_DIR)libmod.a
LIB_TOOLS=$(LIB_DIR)libtools.a

# Entries for library

MOD = $(LIB_MOD)(mod_cmfgen.o)\
      $(LIB_MOD)(ang_qw_mod.o)\
      $(LIB_MOD)(chg_exch_mod_v3.o)\
      $(LIB_MOD)(cmf_sob_mod.o)\
      $(LIB_MOD)(control_variable_mod.o)\
      $(LIB_MOD)(line_mod.o)\
      $(LIB_MOD)(line_vec_mod.o)\
      $(LIB_MOD)(mod_lev_dis_blk.o)\
      $(LIB_MOD)(nuc_iso_mod.o)\
      $(LIB_MOD)(opac_mod.o)\
      $(LIB_MOD)(radiation_mod.o)\
      $(LIB_MOD)(steq_data_mod.o)\
      $(LIB_MOD)(two_phot_mod.o)\
      $(LIB_MOD)(var_rad_mod.o)

OBJ = $(LIB)(adjust_r_grid_v4.o)\
      $(LIB)(adjust_den_vecs.o)\
      $(LIB)(check_param_consistency.o)\
      $(LIB)(comp_j_blank.o)\
      $(LIB)(comp_grey_v2.o)\
      $(LIB)(comp_grey_v4.o)\
      $(LIB)(comp_opac.o)\
      $(LIB)(comp_var_opac.o)\
      $(LIB)(do_var_cont.o)\
      $(LIB)(eval_lte_v5.o)\
      $(LIB)(grey_t_iterate.o)\
      $(LIB)(init_opac_mod.o)\
      $(LIB)(init_line_opac_var_v2.o)\
      $(LIB)(rd_control_variables.o)\
      $(LIB)(rd_cont_j.o)\
      $(LIB)(set_line_opac.o)\
      $(LIB)(set_ang_qw.o)\
      $(LIB)(set_frequency_grid_v2.o)\
      $(LIB)(set_new_model_estimates.o)\
      $(LIB)(set_abund_clump.o)\
      $(LIB)(solve_for_pops.o)\
      $(LIB)(sub_cmf_line.o)\
      $(LIB)(sub_sob_line_v3.o)\
      $(LIB)(sup_to_full_v4.o)

all : $(LIB)
	cp -f $(INSTALL_DIR)perm_mod_file jnk.mod
	mv -f *.mod $(MOD_DIR)
	rm -f *.o

$(LIB) : $(MOD) $(OBJ)	

# General rules to make obects in a library and object files

$(LIB_MOD)(%.o) : %.f 
	$(F90) -c $(FG) $< 
	ar ruv $(LIB_MOD) $*.o

$(LIB)(%.o) : %.f 
	$(F90) -c $(FG) $< 
	ar ruv $(LIB) $*.o

# Specific dependency of source files which uses external 
# fortran moudles. (The order of compilation matters.)

$(LIB)(adjust_r_grid_v4.o) : adjust_r_grid_v4.f $(LIB_MOD)(mod_cmfgen.o) $(LIB_TOOLS)(update_keyword.o)
	$(F90) -c $(FG) adjust_r_grid_v4.f
	ar ruv $(LIB) adjust_r_grid_v4.o

$(LIB)(check_param_consistency.o) : check_param_consistency.f $(LIB_MOD)(control_variable_mod.o)
	$(F90) -c $(FG) check_param_consistency.f
	ar ruv $(LIB) check_param_consistency.o

$(LIB)(comp_grey_v2.o) : comp_grey_v2.f $(LIB_MOD)(control_variable_mod.o) $(LIB_MOD)(mod_cmfgen.o)\
                                  $(LIB_MOD)(ang_qw_mod.o)
	$(F90) -c $(FG) comp_grey_v2.f
	ar ruv $(LIB) comp_grey_v2.o

$(LIB)(comp_grey_v4.o) : comp_grey_v4.f $(LIB_MOD)(control_variable_mod.o) $(LIB_MOD)(mod_cmfgen.o)\
                                  $(LIB_MOD)(ang_qw_mod.o)
	$(F90) -c $(FG) comp_grey_v4.f
	ar ruv $(LIB) comp_grey_v4.o

$(LIB)(comp_j_blank.o) : comp_j_blank.f $(LIB)(opac_mod.o) $(LIB_MOD)(control_variable_mod.o) $(LIB_MOD)(radiation_mod.o) $(LIB_MOD)(mod_cmfgen.o) $(LIB_MOD)(ang_qw_mod.o)
	$(F90) -c $(FG) comp_j_blank.f
	ar ruv $(LIB) comp_j_blank.o

$(LIB)(comp_opac.o) : comp_opac.f $(LIB)(opac_mod.o) $(LIB_MOD)(control_variable_mod.o) $(LIB_MOD)(mod_cmfgen.o)
	$(F90) -c $(FG) comp_opac.f
	ar ruv $(LIB) comp_opac.o

$(LIB)(comp_var_opac.o) : comp_var_opac.f $(LIB)(opac_mod.o) $(LIB_MOD)(control_variable_mod.o) $(LIB_MOD)(mod_cmfgen.o)
	$(F90) -c $(FG) comp_var_opac.f
	ar ruv $(LIB) comp_var_opac.o

$(LIB)(do_var_cont.o) : do_var_cont.f $(LIB)(opac_mod.o) $(LIB_MOD)(control_variable_mod.o)\
             $(LIB_MOD)(line_mod.o) $(LIB_MOD)(var_rad_mod.o) $(LIB_MOD)(radiation_mod.o) $(LIB_MOD)(mod_cmfgen.o)\
             $(LIB_MOD)(ang_qw_mod.o)
	$(F90) -c $(FG) do_var_cont.f
	ar ruv $(LIB) do_var_cont.o
                                                                                                                      
$(LIB)(eval_lte_v5.o) : eval_lte_v5.f $(LIB_MOD)(mod_cmfgen.o)
	$(F90) -c $(FG) eval_lte_v5.f
	ar ruv $(LIB) eval_lte_v5.o

$(LIB)(grey_t_iterate.o) : grey_t_iterate.f $(LIB)(opac_mod.o) $(LIB_MOD)(control_variable_mod.o)\
             $(LIB_MOD)(line_mod.o) $(LIB_MOD)(line_vec_mod.o) $(LIB_MOD)(mod_cmfgen.o) $(LIB_MOD)(ang_qw_mod.o)
	$(F90) -c $(FG) grey_t_iterate.f
	ar ruv $(LIB) grey_t_iterate.o
                                                                                                                      
$(LIB)(init_line_opac_var_v2.o) : init_line_opac_var_v2.f $(LIB)(opac_mod.o) $(LIB_MOD)(control_variable_mod.o)\
             $(LIB_MOD)(line_mod.o) $(LIB_MOD)(mod_cmfgen.o) $(LIB_MOD)(var_rad_mod.o)
	$(F90) -c $(FG) init_line_opac_var_v2.f
	ar ruv $(LIB) init_line_opac_var_v2.o
                                                                                                                      
$(LIB)(rd_cont_j.o) : rd_cont_j.f $(LIB)(opac_mod.o) $(LIB_MOD)(radiation_mod.o) $(LIB_MOD)(mod_cmfgen.o)\
             $(LIB_MOD)(ang_qw_mod.o)
	$(F90) -c $(FG) rd_cont_j.f
	ar ruv $(LIB) rd_cont_j.o
                                                                                                                      
$(LIB)(rd_control_variables.o) : rd_control_variables.f $(LIB_MOD)(control_variable_mod.o) $(LIB_MOD)(mod_cmfgen.o)
	$(F90) -c $(FG) rd_control_variables.f
	ar ruv $(LIB) rd_control_variables.o

$(LIB)(set_abund_clump.o) : set_abund_clump.f $(LIB_MOD)(control_variable_mod.o) $(LIB_MOD)(mod_cmfgen.o)
	$(F90) -c $(FG) set_abund_clump.f
	ar ruv $(LIB) set_abund_clump.o
                                                                                                                      
$(LIB)(set_ang_qw.o) : set_ang_qw.f $(LIB_MOD)(ang_qw_mod.o)
	$(F90) -c $(FG) set_ang_qw.f
	ar ruv $(LIB) set_ang_qw.o

$(LIB)(set_frequency_grid_v2.o) : set_frequency_grid_v2.f $(LIB_MOD)(control_variable_mod.o)\
              $(LIB_MOD)(line_vec_mod.o) $(LIB_MOD)(mod_cmfgen.o) 
	$(F90) -c $(FG) set_frequency_grid_v2.f
	ar ruv $(LIB) set_frequency_grid_v2.o
                                                                                                                      
$(LIB)(set_line_opac.o) : set_line_opac.f $(LIB)(opac_mod.o) $(LIB_MOD)(control_variable_mod.o)\
              $(LIB_MOD)(line_mod.o) $(LIB_MOD)(mod_cmfgen.o) $(LIB_MOD)(line_vec_mod.o)
	$(F90) -c $(FG) set_line_opac.f
	ar ruv $(LIB) set_line_opac.o
                                                                                                                      
$(LIB)(set_new_model_estimates.o) : set_new_model_estimates.f $(LIB)(opac_mod.o) $(LIB_MOD)(control_variable_mod.o)\
              $(LIB_MOD)(line_mod.o) $(LIB_MOD)(line_vec_mod.o) $(LIB_MOD)(mod_cmfgen.o) $(LIB_MOD)(ang_qw_mod.o)
	$(F90) -c $(FG) set_new_model_estimates.f
	ar ruv $(LIB) set_new_model_estimates.o
                                                                                                                      
$(LIB)(solve_for_pops.o) : solve_for_pops.f $(LIB_MOD)(control_variable_mod.o) $(LIB_MOD)(mod_cmfgen.o)\
                                  $(LIB_MOD)(ang_qw_mod.o)
	$(F90) -c $(FG) solve_for_pops.f
	ar ruv $(LIB) solve_for_pops.o

$(LIB)(sub_cmf_line.o) : sub_cmf_line.f LINEGEN.INC $(LIB)(opac_mod.o) $(LIB_MOD)(control_variable_mod.o)\
              $(LIB_MOD)(line_mod.o) $(LIB_MOD)(var_rad_mod.o) $(LIB_MOD)(radiation_mod.o) $(LIB_MOD)(mod_cmfgen.o)\
              $(LIB_MOD)(cmf_sob_mod.o) $(LIB_MOD)(steq_data_mod.o) $(LIB_MOD)(ang_qw_mod.o) $(LIB_MOD)(line_vec_mod.o)
	$(F90) -c $(FG) sub_cmf_line.f
	ar ruv $(LIB) sub_cmf_line.o
                                                                                                                      
$(LIB)(sub_sob_line_v3.o) : sub_sob_line_v3.f $(LIB)(opac_mod.o) $(LIB_MOD)(control_variable_mod.o)\
              $(LIB_MOD)(line_mod.o) $(LIB_MOD)(var_rad_mod.o) $(LIB_MOD)(radiation_mod.o) $(LIB_MOD)(mod_cmfgen.o)\
              $(LIB_MOD)(cmf_sob_mod.o) $(LIB_MOD)(steq_data_mod.o) $(LIB_MOD)(ang_qw_mod.o) $(LIB_MOD)(line_vec_mod.o)
	$(F90) -c $(FG) sub_sob_line_v3.f
	ar ruv $(LIB) sub_sob_line_v3.o
                                                                                                                      
$(LIB)(sup_to_full_v4.o) : sup_to_full_v4.f $(LIB_MOD)(mod_cmfgen.o)
	$(F90) -c $(FG) sup_to_full_v4.f
	ar ruv $(LIB) sup_to_full_v4.o

clean :
	rm -f *.o 
