#
# Include local system definitions

include ../Makefile_definitions

# Library where object files will/are stored

LIB=$(LIB_DIR)libtools.a
#
# Entries for library
#
OBJ = $(LIB)(usr_option.o)\
      $(LIB)(constants_mod.o)\
      $(LIB)(condit_string.o)\
      $(LIB)(count_occur.o)\
      $(LIB)(form_sp_num.o)\
      $(LIB)(gen_in.o)\
      $(LIB)(get_lu.o)\
      $(LIB)(mod_color_pen_def.o)\
      $(LIB)(new_gen_in.o)\
      $(LIB)(get_indx.o)\
      $(LIB)(read_keyword.o)\
      $(LIB)(rd_sing_vec_rvtj.o)\
      $(LIB)(set_case_fn.o)\
      $(LIB)(sve_file.o)\
      $(LIB)(update_keyword.o)\
      $(LIB)(upper_case.o)\
      $(LIB)(usr_hidden.o)\
      $(LIB)(wr_col_str.o)\
      $(LIB)(wr_fmt_string.o)\
      $(LIB)(wr_string.o)


all : $(LIB)
	cp -f $(INSTALL_DIR)perm_mod_file jnk.mod
	mv -f *.mod $(MOD_DIR)
	rm -f *.o
$(LIB) : $(OBJ)	

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f 
	$(F90) -c $(FG)  $<
	ar ruv $(LIB) $*.o

$(LIB)(usr_option.o) : usr_option.f $(LIB)(wr_string.o) 
$(LIB)(usr_hidden.o) : usr_hidden.f $(LIB)(wr_string.o) $(LIB)(sve_file.o)

#
#	$(F90) -c $(FG) usr_option.f 
#	ar ruv $(LIB) usr_option.o
#	$(F90) -c $(FG) usr_hidden.f 
#	ar ruv $(LIB) usr_hidden.o
# 
# basic dependency
#

# specific dependency of source files which uses external 
# fortran moudles. (The order of compilation matters.)

#char_functions_class.o) \
#             libplot.a(pixplot_module.o) 


clean :
	rm -f *.o 

clean_lib :
	rm -f *.a
