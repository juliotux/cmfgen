#
# Put these aliases in you .cshrc file. If in your directory,
# you will need to replace hillier by you username
##############################################################
##############################################################
#
# cmfdist will be defined for general use (e.g. cd cmfdist)
# Change to you local distribution directory

set cmfdist = ~jdh/cur_cmf
set CMFDIST = ~jdh/cur_cmf
set atomic  = ~jdh/atomic 
set ATOMIC  = ~jdh/atomic 
set sn_atomic  = ~jdh/sn_atomic 
set SN_ATOMIC  = ~jdh/sn_atomic 


# Used to assign help fles for DISPGEN, PLT_SPEC etc
# Creates soft links.
# Must be executed inplotting directiry.

alias astxt $cmfdist/com/assign_txt_files.sh

#
# Executable shorthands

alias cmfgen         $cmfdist/exe/cmfgen.exe
alias cmf_flux       $cmfdist/exe/cmf_flux.exe
alias dispgen        $cmfdist/exe/dispgen.exe
alias plt_spec       $cmfdist/exe/plt_spec.exe
alias plt_scr        $cmfdist/exe/plt_scr.exe
alias plt_jh         $cmfdist/exe/plt_jh.exe

alias append_dc      $cmfdist/exe/append_dc.exe
alias rewrite_dc     $cmfdist/exe/rewrite_dc.exe
alias n_col_merge    $cmfdist/exe/n_col_merge.exe
alias n_pair_merge   $cmfdist/exe/n_pair_merge.exe

#
# Performs an NG acceleration operating on the SCRTEMP file

alias do_ng        $cmfdist/exe/do_ng.exe


#
# Used to delete uneccesary files (e.g. Scratch files) after Model has 
# completed. 

alias clean        $cmfdist/com/clean.sh

#
# Used to rename *OUT files to *_IN

alias out2in       $cmfdist/com/out2in.sh

#
# Used to name *.for files to *.f

alias for2f        $cmfdist/com/for_to_f.sh

#
# Used to name xxx.inc file to XXX.INC

alias inc2inc      $cmfdist/com/inc2inc.sh

##################################################
#
# Used to compare *.f files in one directory to those in another

alias for_dif        $cmfdist/com/for_dif.sh

#
# Used to compare *.INC files in one directory to those in another

alias INC_dif        $cmfdist/com/INC_dif.sh

#
# Used to copy files required for the generation of a new model from a 
# a directory containg a completed model.

alias cpmod       $cmfdist/com/cpmod.sh

# General aliases to clean a model directory. Be careful --- these might
# delete needed files if used incorrectly.

alias rmrrr     'rm -vf *PRRR'
alias rmin      'rm -vf *_IN'
alias dfort     'rm -vf fort.*'
alias dsve      'rm -vf *.sve'
alias dlog      'rm -vf *.log'
alias dscratch  'rm -vf CSCRATCH* BCSCRATCH* DSCRATCH*'
#
# Removes soft links from current directory
#
alias rmlinks         "find * -type l -maxdepth 0 -exec rm -vf {} ';'"
#
# Removes soft links from curent directory, and directories below it.
#
alias rm_all_links   "find * -type l -maxdepth 10 -exec rm -vf {} ';'"

