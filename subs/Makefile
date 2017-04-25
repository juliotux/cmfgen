include ../Makefile_definitions

# Library where object files will/are stored

LIB=$(LIB_DIR)libsubs.a
LIB_MOD=$(LIB_DIR)libmod.a

OPTION=$(FH)

#
# Entries for library
#
OBJ = $(LIB)(add_par_to_full.o)\
      $(LIB)(adjust_line_freq_v2.o)\
      $(LIB)(ave_flips.o)\
      $(LIB)(undo_it.o)\
      $(LIB)(blkband.o)\
      $(LIB)(chk_eq_num.o)\
      $(LIB)(chk_col_name.o)\
      $(LIB)(cmf_form_sol.o)\
      $(LIB)(cmf_form_sol_v2.o)\
      $(LIB)(cnvt_fr_dc.o)\
      $(LIB)(cnvt_fr_dc_v2.o)\
      $(LIB)(colcool_sl_v3.o)\
      $(LIB)(comp_j_conv_v2.o)\
      $(LIB)(comp_obs_v2.o)\
      $(LIB)(comp_var_jrec.o)\
      $(LIB)(det_main_cont_freq.o)\
      $(LIB)(dop_pro.o)\
      $(LIB)(dp_zero.o)\
      $(LIB)(dspheredchi.o)\
      $(LIB)(edd_j_var_v4.o)\
      $(LIB)(edd_j_var_v6.o)\
      $(LIB)(eddline_var.o)\
      $(LIB)(elec_prep.o)\
      $(LIB)(equal.o)\
      $(LIB)(esopac.o)\
      $(LIB)(estau.o)\
      $(LIB)(estau_v2.o)\
      $(LIB)(eval_adiabatic.o)\
      $(LIB)(eval_ed.o)\
      $(LIB)(ew_format.o)\
      $(LIB)(ex_e1x_fun.o)\
      $(LIB)(exponx.o)\
      $(LIB)(extend3opac.o)\
      $(LIB)(extend_opac.o)\
      $(LIB)(extend_vtsigma.o)\
      $(LIB)(fdg_f_to_s_ns_v1.o)\
      $(LIB)(fg_comp.o)\
      $(LIB)(fg_j_cmf_v9.o)\
      $(LIB)(fg_j_cmf_v10.o)\
      $(LIB)(fix_dchi.o)\
      $(LIB)(fixpop.o)\
      $(LIB)(formsol.o)\
      $(LIB)(fqcomp_ibc.o)\
      $(LIB)(fqcomp_ibc_v2.o)\
      $(LIB)(fqnew.o)\
      $(LIB)(full_to_sup.o)\
      $(LIB)(gamma.o)\
      $(LIB)(gausel.o)\
      $(LIB)(gbf.o)\
      $(LIB)(gen_omega_rd_v2.o)\
      $(LIB)(genaccel.o)\
      $(LIB)(genangqw.o)\
      $(LIB)(genosc_v5.o)\
      $(LIB)(genosc_v6.o)\
      $(LIB)(get_col_ref.o)\
      $(LIB)(getelec_v2.o)\
      $(LIB)(gff.o)\
      $(LIB)(gnthydl.o)\
      $(LIB)(htrpwgt.o)\
      $(LIB)(hweight.o)\
      $(LIB)(hydanggaunt.o)\
      $(LIB)(hydcrossl.o)\
      $(LIB)(hydro_terms.o)\
      $(LIB)(hydro_terms_v5.o)\
      $(LIB)(hyper.o)\
      $(LIB)(ichrlen.o)\
      $(LIB)(impar.o)\
      $(LIB)(indexx.o)\
      $(LIB)(init_temp_v2.o)\
      $(LIB)(ins_line_obs_v2.o)\
      $(LIB)(ins_line_obs_v3.o)\
      $(LIB)(ins_line_obs_v5.o)\
      $(LIB)(ins_line_v6.o)\
      $(LIB)(integrate.o)\
      $(LIB)(interp_knowncoeef.o)\
      $(LIB)(jfeau_ibc.o)\
      $(LIB)(jfeau_ibc_v2.o)\
      $(LIB)(jgrey.o)\
      $(LIB)(jgrey_with_fvt.o)\
      $(LIB)(jtrpwgt.o)\
      $(LIB)(jweight.o)\
      $(LIB)(kivarnm.o)\
      $(LIB)(ktrpwgt.o)\
      $(LIB)(kweight.o)\
      $(LIB)(lam_vac.o)\
      $(LIB)(lamvacair.o)\
      $(LIB)(last2rays.o)\
      $(LIB)(last2vki.o)\
      $(LIB)(lin_interp.o)\
      $(LIB)(locsolut.o)\
      $(LIB)(lte_pop_sl.o)\
      $(LIB)(lte_pop_sl_v2.o)\
      $(LIB)(ltepop.o)\
      $(LIB)(ltepop_wld_v1.o)\
      $(LIB)(ltepop_wld_v2.o)\
      $(LIB)(lum_from_eta.o)\
      $(LIB)(lum_from_eta_v2.o)\
      $(LIB)(marupv.o)\
      $(LIB)(mat5pen.o)\
      $(LIB)(mess.o)\
      $(LIB)(mod_trap_derivatives.o)\
      $(LIB)(modifyba.o)\
      $(LIB)(mom_j_cmf_v5.o)\
      $(LIB)(mom_j_cmf_v6.o)\
      $(LIB)(mom_j_cmf_v7.o)\
      $(LIB)(momjbar.o)\
      $(LIB)(mon_int_funs_v2.o)\
      $(LIB)(mon_int_ins_v1.o)\
      $(LIB)(mon_interp.o)\
      $(LIB)(mon_interp_fast.o)\
      $(LIB)(mon_interp_fast_v2.o)\
      $(LIB)(monotonic_integrate.o)\
      $(LIB)(mult2d.o)\
      $(LIB)(multvec.o)\
      $(LIB)(nagminsol.o)\
      $(LIB)(new_r_scale_v2.o)\
      $(LIB)(newgsit.o)\
      $(LIB)(newjsold.o)\
      $(LIB)(newpertjd.o)\
      $(LIB)(newvkimd.o)\
      $(LIB)(newwmat.o)\
      $(LIB)(ngaccel.o)\
      $(LIB)(ngaccel_arb_ord.o)\
      $(LIB)(nordangqw.o)\
      $(LIB)(nordflux.o)\
      $(LIB)(nordtau.o)\
      $(LIB)(ntrpwgt.o)\
      $(LIB)(nweight.o)\
      $(LIB)(occupation_prob.o)\
      $(LIB)(omega_gen_v2.o)\
      $(LIB)(open_dir_acc_v1.o)\
      $(LIB)(par_fun_v4.o)\
      $(LIB)(pertjfeau_ibc.o)\
      $(LIB)(phys_con.o)\
      $(LIB)(poptoion.o)\
      $(LIB)(prrr_sl_v6.o)\
      $(LIB)(prrr_lowt_v1.o)\
      $(LIB)(qfromf.o)\
      $(LIB)(qkim.o)\
      $(LIB)(quad_x_gen_v5.o)\
      $(LIB)(quad_x_lit_v3.o)\
      $(LIB)(rd_f_to_s_ids.o)\
      $(LIB)(rd_f_to_s_ids_v2.o)\
      $(LIB)(rd_free_val.o)\
      $(LIB)(rd_log.o)\
      $(LIB)(rd_t_ed.o)\
      $(LIB)(rdgendie_v4.o)\
      $(LIB)(readba.o)\
      $(LIB)(regrid_b_on_ne.o)\
      $(LIB)(regrid_b_on_sph_tau.o)\
      $(LIB)(regrid_dchi.o)\
      $(LIB)(regrid_h.o)\
      $(LIB)(regrid_tx_r.o)\
      $(LIB)(regrid_t_ed.o)\
      $(LIB)(regrid_t_ed_v2.o)\
      $(LIB)(regrid_log_dc_v1.o)\
      $(LIB)(rext_coef_v2.o)\
      $(LIB)(rite_asc.o)\
      $(LIB)(rite_atmdes_v2.o)\
      $(LIB)(rite_gam_v2.o)\
      $(LIB)(rv_sn_model_v2.o)\
      $(LIB)(s15adf.o)\
      $(LIB)(same_n.o)\
      $(LIB)(scale_grey.o)\
      $(LIB)(scale_pops.o)\
      $(LIB)(scr_read_v2.o)\
      $(LIB)(set_case.o)\
      $(LIB)(set_cont_freq.o)\
      $(LIB)(set_cont_freq_v3.o)\
      $(LIB)(set_cont_freq_v4.o)\
      $(LIB)(set_dc_lte_v2.o)\
      $(LIB)(set_edge_freq_v2.o)\
      $(LIB)(set_t_ed_grid.o)\
      $(LIB)(set_x_freq.o)\
      $(LIB)(set_gm_all.o)\
      $(LIB)(set_z_pop.o)\
      $(LIB)(specify_it_cycle_v2.o)\
      $(LIB)(simq.o)\
      $(LIB)(sobew.o)\
      $(LIB)(sobew_grad.o)\
      $(LIB)(sobew_grad_v2.o)\
      $(LIB)(sobjbar_sim.o)\
      $(LIB)(solveba_v3.o)\
      $(LIB)(starfin.o)\
      $(LIB)(starnew.o)\
      $(LIB)(starpcyg_v3.o)\
      $(LIB)(starpcyg_v4.o)\
      $(LIB)(starrave.o)\
      $(LIB)(steqne.o)\
      $(LIB)(storeba.o)\
      $(LIB)(sup_to_full_v3.o)\
      $(LIB)(sup_to_full_v3b.o)\
      $(LIB)(tau.o)\
      $(LIB)(tcompd.o)\
      $(LIB)(term_io.o)\
      $(LIB)(tfeau.o)\
      $(LIB)(thomas.o)\
      $(LIB)(thomas_pone.o)\
      $(LIB)(thomas_rh.o)\
      $(LIB)(torscl.o)\
      $(LIB)(torscl_v2.o)\
      $(LIB)(torscl_v3.o)\
      $(LIB)(transposeba.o)\
      $(LIB)(tuvghd.o)\
      $(LIB)(tuvghd_rh.o)\
      $(LIB)(ungrid.o)\
      $(LIB)(up_tx_tvx.o)\
      $(LIB)(up_tx_tvx_noeps.o)\
      $(LIB)(up_tx_tvx_v2.o)\
      $(LIB)(up_tx_tvx_ext_v1.o)\
      $(LIB)(update_gam.o)\
      $(LIB)(updateu.o)\
      $(LIB)(uptvx_edd.o)\
      $(LIB)(uptx_edd.o)\
      $(LIB)(upvnot.o)\
      $(LIB)(var_formsol.o)\
      $(LIB)(var_mom_j_cmf_v6.o)\
      $(LIB)(var_mom_j_cmf_v7.o)\
      $(LIB)(var_mom_j_cmf_v8.o)\
      $(LIB)(var_mom_jext_cmf_v2.o)\
      $(LIB)(var_momjbar.o)\
      $(LIB)(var_op_v4.o)\
      $(LIB)(var_x_opa_eta_v2.o)\
      $(LIB)(varlamki.o)\
      $(LIB)(vkifeau_ibc.o)\
      $(LIB)(voigt_1p0.o)\
      $(LIB)(vsebyj_multi_v4.o)\
      $(LIB)(weight.o)\
      $(LIB)(wr2d_ma.o)\
      $(LIB)(wr2d.o)\
      $(LIB)(wr2d_v2.o)\
      $(LIB)(wr_ad_cool.o)\
      $(LIB)(wr_art_heat.o)\
      $(LIB)(wr_cooling_term.o)\
      $(LIB)(wr_x_info.o)\
      $(LIB)(wrcoolgen_v2.o)\
      $(LIB)(write_j_cmf_err.o)\
      $(LIB)(writedc.o)\
      $(LIB)(writedc_v2.o)\
      $(LIB)(writedc_v3.o)\
      $(LIB)(write_direct_info.o)\
      $(LIB)(write_direct_info_v3.o)\
      $(LIB)(write_vec.o)\
      $(LIB)(writv.o)\
      $(LIB)(writv_v2.o)\
      $(LIB)(wrrecomchk_v4.o)\
      $(LIB)(x02ajf.o)\
      $(LIB)(x_rrr_cool_v6.o)\
      $(LIB)(xcross.o)\
      $(LIB)(xvecd.o)\
      $(LIB)(zalongp.o)\
      $(LIB)(zero.o)\
      $(LIB)(zero_2d_mat.o)

all : $(LIB) UNIX
	cp -f $(INSTALL_DIR)perm_mod_file jnk.mod
	mv -f *.mod $(MOD_DIR)
	rm -f *.o
$(LIB) : $(OBJ)	

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f 
	$(F90) -c $(OPTION)  $<
	ar ruv $(LIB) $*.o

# 
# Specific dependency of source files which uses external 
# fortran moudles. (The order of compilation matters.)
#
DIS=$(LIB_MOD)(mod_lev_dis_blk.o)

$(LIB)(dspheredchi.o) : dspheredchi.f $(LIB)(mod_trap_derivatives.o) 
	$(F90) -c $(OPTION) dspheredchi.f
	ar ruv $(LIB) dspheredchi.o

$(LIB)(occupation_prob.o) : occupation_prob.f $(DIS)
	$(F90) -c $(OPTION) occupation_prob.f
	ar ruv $(LIB) occupation_prob.o

$(DIS):
	(cd ../new_main/mod_subs ; make)

UNIX: 
	(cd ../unix; make)

clean :
	rm -f *.o

clean_lib :
	rm -f $(LIB)
