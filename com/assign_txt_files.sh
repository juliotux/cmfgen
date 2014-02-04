#
# Command procedure to perform soflinks in the current directory so that
# the help files for DISPGEN etc are available.
#
# For DISPGEN
#
ln -sf $cmfdist/txt_files/maingen_opt_desc.txt        MAINGEN_OPT_DESC
ln -sf $cmfdist/txt_files/maingen_options.txt         MAINGEN_OPTIONS
#
# For PLT_SPEC
#
ln -sf $cmfdist/txt_files/plt_spec_options.txt        PLT_SPEC_OPTIONS
ln -sf $cmfdist/txt_files/plt_spec_opt_desc.txt       PLT_SPEC_OPT_DESC
ln -sf $cmfdist/txt_files/h2_is_line_list.dat         H2_IS_LINE_LIST
ln -sf $cmfdist/txt_files/hi_is_line_list.dat         HI_IS_LINE_LIST
#
# For PLT_JH
#
ln -sf $cmfdist/txt_files/plt_jh_opt_desc.txt         PLT_JH_OPT_DESC
ln -sf $cmfdist/txt_files/plt_jh_options.txt      PLT_JH_OPTIONS
#
# For WR_F_TO_S
#
ln -sf $cmfdist/txt_files/wr_f_opt_desc.txt          WR_F_OPT_DESC
ln -sf $cmfdist/txt_files/wr_f_options.txt           WR_F_OPTIONS
#
# For MASS_SC
#
ln -sf $cmfdist/txt_files/sol_abund.dat              SOL_ABUND
