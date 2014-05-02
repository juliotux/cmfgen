# Include local system definitions

include ../../../Makefile_definitions

# Library where object files will/are stored
# FX =  -Mextend -Mbackslash    -Ktrap=fp

LIB=$(LIB_DIR)libdev.a
LIB_MOD=$(LIB_DIR)libmod.a

#
# Entries in a library
#
OBJ = $(LIB)(rd_chg_exch_v3.o)\
      $(LIB)(do_penning_ion.o)\
      $(LIB)(eval_chg_rates_v3.o)\
      $(LIB)(set_chg_exch_v4.o)\
      $(LIB)(set_chg_lev_id_v4.o)\
      $(LIB)(verify_chg_exch_v3.o)\
      $(LIB)(steq_ba_chg_exch_v3.o)\
      $(LIB)(wr_chg_cool_v3.o)

all : $(LIB)
	rm -f *.o

$(LIB) : $(OBJ)	

CG =  $(LIB_DIR)libmod.a(chg_exch_mod_v3.o)
ST =  $(LIB_DIR)libmod.a(steq_data_mod.o)
CMF = $(LIB_DIR)libmod.a(mod_cmfgen.o)

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f $(CG)
	$(F90) -c $(FG) $< 
	ar ruv $(LIB) $*.o

# Specific rules

$(LIB)(steq_ba_chg_exch_v3.o) : steq_ba_chg_exch_v3.f $(CG) $(ST) 
	$(F90) -c $(FG) steq_ba_chg_exch_v3.f
	ar ruv $(LIB) steq_ba_chg_exch_v3.o

$(LIB)(set_chg_lev_id_v4.o) : set_chg_lev_id_v4.f $(CG) $(CMF)
	$(F90) -c $(FG) set_chg_lev_id_v4.f
	ar ruv $(LIB) set_chg_lev_id_v4.o

$(CG):
	(cd ../mod_subs ; make)

$(ST):
	(cd ../mod_subs ; make)

$(MOD):
	(cd ../mod_subs ; make)

clean :
	rm -f *.o 

clean_lib :
	rm $(LIB)
