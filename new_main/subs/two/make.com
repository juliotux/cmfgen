$
$set verify
$
$fort:=f90/extend
$!fort:=F90/DEBUG/LIS/NOOPT/CHECK=(BOUNDS,OVER)/G_FLOAT/SHOW=INCLUDE/EXTEND
$
$fort two_phot_mod.f
$
$fort steq_ba_two_phot_rate_v3.f /include=[-]
$
$fort two_phot_rate.f
$fort two_phot_opac.f
$fort set_two_phot_v2.f
$fort two_phot_var_opac.f
$
$lib [-]tsts *.obj
$
$set noverify
