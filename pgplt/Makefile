# Include local system definitions

include ../Makefile_definitions

# Library where object files will/are stored

LIB=$(LIB_DIR)libpgplt.a
LIB_MOD=$(LIB_DIR)libmod.a

#
# Entries in a library
#
OBJ = $(LIB)(mod_curve_data.o)\
      $(LIB)(amoeba.o)\
      $(LIB)(amotry.o)\
      $(LIB)(change_color.o)\
      $(LIB)(change_pen.o)\
      $(LIB)(curve.o)\
      $(LIB)(cut_points_from_plot.o)\
      $(LIB)(define_more_pens.o)\
      $(LIB)(do_fill.o)\
      $(LIB)(do_many_ew.o)\
      $(LIB)(do_mult_gf.o)\
      $(LIB)(do_vec_op.o)\
      $(LIB)(draw_right_y_axis.o)\
      $(LIB)(ed_gaus_fit.o)\
      $(LIB)(ew_group.o)\
      $(LIB)(find_lines.o)\
      $(LIB)(fntick.o)\
      $(LIB)(gaus_romb.o)\
      $(LIB)(gaus_fit_er.o)\
      $(LIB)(get_gramon_min_max.o)\
      $(LIB)(get_pgi_file_info.o)\
      $(LIB)(get_titles.o)\
      $(LIB)(gramon_pgplot.o)\
      $(LIB)(hist_adj.o)\
      $(LIB)(justify_convert_v2.o)\
      $(LIB)(modify_pgi_ps.o)\
      $(LIB)(mon_num.o)\
      $(LIB)(monbord_v3.o)\
      $(LIB)(moninside.o)\
      $(LIB)(put_text.o)\
      $(LIB)(simp_ew.o)\
      $(LIB)(strip_slash_p.o)\
      $(LIB)(shrink_vectors.o)\
      $(LIB)(spacing.o)\
      $(LIB)(vectorpen.o)\
      $(LIB)(gaus_fit_data.o)\
      $(LIB)(do_gaus_fit.o)\
      $(LIB)(wr_gaus_fit.o)\
      $(LIB)(write_line_lambdas.o)\
      $(LIB)(gaus_fit_func.o)

all : $(LIB)
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

GEN = ../lib/libtools.a(gen_in.o)
COL=  ../lib/libtools.a(mod_color_pen_def.o)

$(LIB)(change_color.o)    : change_color.f  $(GEN)
$(LIB)(change_pen.o)      : change_pen.f    $(GEN)

$(LIB)(gramon_pgplot.o)   : gramon_pgplot.f $(LIB)(mod_curve_data.o) $(GEN) $(COL) 
$(LIB)(curve.o)           : curve.f         $(LIB)(mod_curve_data.o) $(GEN)
$(LIB)(vectorpen.o)       : vectorpen.f     $(LIB)(mod_curve_data.o) $(GEN)

$(LIB)(wr_gaus_fit.o)     : wr_gaus_fit.f   $(LIB)(gaus_fit_data.o)
$(LIB)(do_gaus_fit.o)     : do_gaus_fit.f   $(LIB)(gaus_fit_data.o)
$(LIB)(do_mult_gf.o)      : do_mult_gf.f    $(LIB)(gaus_fit_data.o)
$(LIB)(gaus_fit_func.o)   : gaus_fit_func.f $(LIB)(gaus_fit_data.o)
$(LIB)(gaus_fit_er.o)     : gaus_fit_er.f   $(LIB)(gaus_fit_data.o)
$(LIB)(gaus_romb.o)       : gaus_romb.f     $(LIB)(gaus_fit_data.o)
$(LIB)(ed_gaus_fit.o)     : ed_gaus_fit.f   $(LIB)(gaus_fit_data.o)
$(LIB)(find_lines.o)      : find_lines.f    $(LIB)(gaus_fit_data.o)

$(GEN):
	(cd ../tools; make )

$(COL):
	(cd ../tools; make )

clean :
	rm -f *.o

clean_lib :
	rm -f $(LIB)
