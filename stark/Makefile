# Include local system definitions

include ../Makefile_definitions

# Library where object files will/are stored

LIB=$(LIB_DIR)libsubs.a
LIB2=$(LIB_DIR)libsubs.a
                         
OBJ = prof_mod.o init_prof_module.o set_prof_v4.o set_prof_v5.o \
      get_prof_storage_limits.o griem_stark_mod.o griem_v2.o set_prof_limits_v3.o set_prof_limits_v4.o \
      conv_strk_v2.o rd_strk_list.o strk_bs_hhe.o \
      rd_bs_strk_tab.o voigtn.o

NEWF = strk_hei_ir.o strk_ds.o

all : $(LIB) $(LIB2)
	rm -f *.o
	cp $(INSTALL_DIR)perm_mod_file  jnk.mod
	mv -f *.mod $(MOD_DIR)

$(LIB) : $(LIB)($(OBJ))
$(LIB2) : $(LIB2)($(NEWF))

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f
	$(F90) -c $(FG)  $<
	ar ruv $(LIB) $*.o
                           
$(LIB2)(%.o) : %.f90
	$(F90) -c $(FFREE)  $<
	ar ruv $(LIB) $*.o
                           
# Dependencies

$(LIB)(prof_mod.o)             : prof_mod.f 
$(LIB)(init_prof_module.o)     : init_prof_module.f   $(LIB)(prof_mod.o) 
$(LIB)(set_prof_v4.o)          : set_prof_v4.f        $(LIB)(prof_mod.o)  $(LIB)(rd_strk_list.o)
$(LIB)(set_prof_limits_v2.o)   : set_prof_limits_v2.f $(LIB)(prof_mod.o) $(LIB)(rd_strk_list.o)
$(LIB)(strk_bs_hhe.o)          : strk_bs_hhe.f        $(LIB)(rd_bs_strk_tab.o)
                                         
.dummy : clean
clean :
	rm -f *.o
	rm -f $(TARGET) $(MASTER)
