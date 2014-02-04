$CD DISK$prog:[cmfmod.burst.a3_r100.obs]
$@[-]BATCH ASS
$copy cmf_flux_param_init.dat cmf_flux_param.dat
$
$RN disk$user:[code.work.obs]cmf_flux_v4
[-]RVTJ.DAT
60.0
$delete/noconfirm eddfactor.dat;*
$delete/noconfirm es_j_conv.dat;*
$rename obsframe.dat obs_fin.dat
$rename hydro.dat hydro_fin.dat
$modify cmf_flux_param.dat "BLANK        [GLOBAL_LINE]" "SOB          [GLOBAL_LINE]"
$modify cmf_flux_param.dat "F            [DO_SOB_LINES]" "T            [DO_SOB_LINES]"
$modify cmf_flux_param.dat "2            [NUM_ES]" "1            [NUM_ES]"
$
$RN disk$user:[code.work.obs]cmf_flux_v4
[-]RVTJ.DAT
60.0
$delete/noconfirm eddfactor.dat;*
$delete/noconfirm es_j_conv.dat;*
$delete/noconfirm j_comp.dat;*
$delete/noconfirm trans_info.dat;*
$delete/noconfirm cfdat_out.dat;*
$ddat
$rename obsframe.dat obs_cont.dat
