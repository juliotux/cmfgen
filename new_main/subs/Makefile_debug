# Include local system definitions

include ../../Makefile_definitions

# Library where object files will/are stored

LIB=$(LIB_DIR)libdev.a

#
# Entries for library
#
OBJ = $(LIB)(add_par_to_full_v2.o)\
      $(LIB)(ba_update_v6.o)\
      $(LIB)(cmf_blkband_v2.o)\
      $(LIB)(create_iv_links_v2.o)\
      $(LIB)(determine_nse.o)\
      $(LIB)(do_ng_accel_band_v1.o)\
      $(LIB)(eval_adiabatic_v2.o)\
      $(LIB)(evalse_qwvj_v6.o)\
      $(LIB)(evalse_x_qwvj_v4.o)\
      $(LIB)(fixpop_in_ba_v2.o)\
      $(LIB)(fixpop_t.o)\
      $(LIB)(generate_full_matrix.o)\
      $(LIB)(rd_store_3int.o)\
      $(LIB)(read_ba_data_v2.o)\
      $(LIB)(set_ba_storage.o)\
      $(LIB)(set_imp_vec.o)\
      $(LIB)(solveba_v6.o)\
      $(LIB)(steq_multi_v7.o)\
      $(LIB)(steqne_v4.o)\
      $(LIB)(steq_advec_v1.o)\
      $(LIB)(store_ba_data_v2.o)\
      $(LIB)(var_op_v7.o)\
      $(LIB)(var_op_v8.o)\
      $(LIB)(var_x_opa_eta_v3.o)\
      $(LIB)(var_x_opa_eta_v4.o)\
      $(LIB)(vsebyj_multi_v6.o)\
      $(LIB)(vsebyj_x_v6.o)\
      $(LIB)(wr_asci_steq.o)\
      $(LIB)(fixpop_in_ba_v3.o)\
      $(LIB)(generate_full_matrix_v2.o)\
      $(LIB)(solveba_v7.o)\
      $(LIB)(read_bcd_mat.o)\
      $(LIB)(cmf_blkband_v3.o)

OBJ2 = $(LIB)(steq_data_mod.o)

all : $(LIB) AUTO CHG TWO
	rm -f *.o
$(LIB) : mod $(OBJ2) $(OBJ)	

mod: $(LIB)(steq_data_mod.o)
chg: $(LIB)(chg_exch_mod_v3.o)

# General rules to make obects in a library and object files

MGEN= ../mod_cmfgen.o
MCHG= $(LIB)(chg_exch_mod_v3.o)

$(LIB)(%.o) : %.f $(LIB)(steq_data_mod.o)  $(GEN)
	$(F90) -c $(FG) -I.. $< 
	ar ruv $(LIB) $*.o

# Specific dependency of source files which uses external 
# fortran moudles. (The order of compilation matters.)

$(LIB)(steq_data_mod.o) : steq_data_mod.f 
	$(F90) -c $(FG) steq_data_mod.f
	ar ruv $(LIB) steq_data_mod.o

$(LIB)(determine_nse.o) : determine_nse.f $(LIB)(steq_data_mod.o) $(MGEN) $(MCHG)
	$(F90) -c $(FG) -Ichg -I.. determine_nse.f
	ar ruv $(LIB) determine_nse.o

$(LIB)(create_iv_links_v2.o) : create_iv_links_v2.f $(LIB)(steq_data_mod.o) $(MGEN)
	$(F90) -c $(FG) -I.. create_iv_links_v2.f
	ar ruv $(LIB) create_iv_links_v2.o

$(LIB)(fixpop_in_ba_v3.o) : fixpop_in_ba_v3.f $(LIB)(steq_data_mod.o) $(MGEN)
	$(F90) -c $(FG) -I.. fixpop_in_ba_v3.f
	ar ruv $(LIB) fixpop_in_ba_v3.o

$(LIB)(generate_full_matrix_v2.o) : generate_full_matrix.f $(LIB)(steq_data_mod.o) $(MGEN)
	$(F90) -c $(FG) -I.. generate_full_matrix_v2.f
	ar ruv $(LIB) generate_full_matrix_v2.o

$(LIB)(set_ba_storage.o) : set_ba_storage.f $(LIB)(steq_data_mod.o) $(MGEN)
	$(F90) -c $(FG) -I.. set_ba_storage.f
	ar ruv $(LIB) set_ba_storage.o

$(MGEN):
	(cd .. ; make mod)

$(MCHG):
	(cd chg ; make mod)

AUTO:
	(cd auto; make)
CHG:
	(cd chg; make)
TWO:
	(cd two; make)

clean :
	rm -f *.o 

clean_lib :
	rm $(LIB)
