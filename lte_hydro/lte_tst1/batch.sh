#!/bin/tcsh

#**********************************************************************
#   Define directories (CAPITAL LETTERS)
#**********************************************************************

# MODEL is model directory
# ATOMIC is atomic data directory

setenv CMFGEN_PROG  $cmfdist/exe/cmfgen_dev.exe

#***********************************************************************
#    Atomic data soft links. All other shell commands should be
#    place after these links. These links are model dependent. Ther
#    are required by CMFGEN, CMF_FLUX, and DISPGEN
#
#    When N_F=N_S, no F_TO_S link is needed.
#***********************************************************************

#*****************************************************************************
# Generic
# -------
#*****************************************************************************
 ln -sf  $ATOMIC/misc/two_phot_data.dat                    TWO_PHOT_DATA
 ln -sf  $ATOMIC/misc/xray_phot_fits.dat                    XRAY_PHOT_FITS
 ln -sf $ATOMIC/misc/rs_xray_fluxes_sol.dat                 RS_XRAY_FLUXES
 ln -sf  $ATOMIC/HYD/I/5dec96/hyd_l_data.dat               HYD_L_DATA
 ln -sf  $ATOMIC/HYD/I/5dec96/gbf_n_data.dat               GBF_N_DATA

#*****************************************************************************
#    Hydrogen
#   ----------
#*****************************************************************************
 ln -sf  $ATOMIC/HYD/I/5dec96/hiphot.dat                   PHOTHI_A
 ln -sf  $ATOMIC/HYD/I/5dec96/hi_osc.dat                   HI_F_OSCDAT
 ln -sf  $ATOMIC/HYD/I/5dec96/hi_f_to_s_15.dat             HI_F_TO_S
 ln -sf  $ATOMIC/HYD/I/5dec96/hicol.dat                    HI_COL_DATA


#*****************************************************************************
#    Helium
#   --------
#*****************************************************************************
 ln -sf  $ATOMIC/HE/II/5dec96/he2phot.dat                  PHOTHe2_A
 ln -sf  $ATOMIC/HE/II/5dec96/he2_osc.dat                  He2_F_OSCDAT
 ln -sf  $ATOMIC/HE/II/5dec96/he2_f_to_s_22.dat            He2_F_TO_S
 ln -sf  $ATOMIC/HE/II/5dec96/he2col.dat                   He2_COL_DATA
#
 ln -sf  $ATOMIC/HE/I/5dec96/heiphot_a7.dat                PHOTHeI_A
 ln -sf  $ATOMIC/HE/I/5dec96/heioscdat_a7.dat              HeI_F_OSCDAT
 ln -sf  $ATOMIC/HE/I/5dec96/hei_f_to_s_a7_ext.dat            HeI_F_TO_S
 ln -sf  $ATOMIC/HE/I/5dec96/heicol.dat                    HeI_COL_DATA

# ln -sf  $ATOMIC/HE/I/5dec96/heiphot_a4.dat                PHOTHeI_A
# ln -sf  $ATOMIC/HE/I/5dec96/heioscdat_a4.dat              HeI_F_OSCDAT
# ln -sf  $ATOMIC/HE/I/5dec96/hei_f_to_s_a4.dat             HeI_F_TO_S
# ln -sf  $ATOMIC/HE/I/5dec96/heicol.dat                    HeI_COL_DATA


#*****************************************************************************
#    Carbon
#   --------
#*****************************************************************************
#
 ln -sf  $ATOMIC/CARB/II/tst/phot_sm_100.dat            PHOTC2_A
 ln -sf  $ATOMIC/CARB/II/tst/phot_data_B.dat            PHOTC2_B
 ln -sf  $ATOMIC/CARB/II/tst/c2osc_rev.dat              C2_F_OSCDAT
 ln -sf  $ATOMIC/CARB/II/tst/c2col.dat                  C2_COL_DATA
 ln -sf  $ATOMIC/CARB/II/tst/f_to_s_104.dat             C2_F_TO_S
 ln -sf  $ATOMIC/CARB/II/tst/c2_auto.dat                C2_AUTO_DATA
#
 ln -sf  $ATOMIC/CARB/III/17oct97/ciiiphot_sm_a.dat          PHOTCIII_A
 ln -sf  $ATOMIC/CARB/III/17oct97/ciiiphot_sm_b.dat          PHOTCIII_B
 ln -sf  $ATOMIC/CARB/III/17oct97/ciiiosc_st_split_big.dat   CIII_F_OSCDAT
 ln -sf  $ATOMIC/CARB/III/17oct97/ciii_f_to_s_split_big.dat  CIII_F_TO_S
 ln -sf  $ATOMIC/CARB/III/17oct97/ciiicol.dat                CIII_COL_DATA
 ln -sf  $ATOMIC/CARB/III/17oct97/dieciii_ic.dat              DIECIII
#
 ln -sf  $ATOMIC/CARB/IV/5dec96/civphot_a12.dat            PHOTCIV_A
 ln -sf  $ATOMIC/CARB/IV/5dec96/civosc_a12_split.dat       CIV_F_OSCDAT
 ln -sf  $ATOMIC/CARB/IV/5dec96/civ_f_to_s_split.dat       CIV_F_TO_S
 ln -sf  $ATOMIC/CARB/IV/5dec96/civcol.dat                 CIV_COL_DATA

 ln -sf  $ATOMIC/CARB/V/9may02/phot_smooth_3000                 PHOTCV_A
 ln -sf  $ATOMIC/CARB/V/9may02/fin_osc                          CV_F_OSCDAT
 ln -sf  $ATOMIC/CARB/V/9may02/f_to_s_nosplit_43                CV_F_TO_S
 ln -sf  $ATOMIC/CARB/V/9may02/col_guess.dat                    CV_COL_DATA

 ln -sf  $ATOMIC/CARB/VI/9may02/phot_smooth_3000              PHOTCSIX_A
 ln -sf  $ATOMIC/CARB/VI/9may02/fin_osc                       CSIX_F_OSCDAT
 ln -sf  $ATOMIC/CARB/VI/9may02/col_guess.dat                 CSIX_COL_DATA
 ln -sf  $ATOMIC/CARB/VI/9may02/f_to_s_30.dat                 CSIX_F_TO_S

#*****************************************************************************
#    Nitrogen
#   ----------
#*****************************************************************************
 ln -sf $ATOMIC/NIT/II/29jul98/n2osc_split.dat             N2_F_OSCDAT
 ln -sf $ATOMIC/NIT/II/29jul98/f_to_s_term.dat             N2_F_TO_S
 ln -sf $ATOMIC/NIT/II/29jul98/n2phot_a.dat                PHOTN2_A
 ln -sf $ATOMIC/NIT/II/29jul98/n2phot_b.dat                PHOTN2_B
 ln -sf $ATOMIC/NIT/II/29jul98/n2phot_c.dat                PHOTN2_C
 ln -sf $ATOMIC/NIT/II/29jul98/n2col.dat                   N2_COL_DATA
 ln -sf $ATOMIC/NIT/II/29jul98/n2die.dat                   DIEN2
#
 ln -sf  $ATOMIC/NIT/III/7jun99/niii_phot_op_3000.dat      PHOTNIII_A
 ln -sf  $ATOMIC/NIT/III/7jun99/niiiosc_op_split.dat       NIII_F_OSCDAT
 ln -sf  $ATOMIC/NIT/III/7jun99/niiicol.dat                NIII_COL_DATA
 ln -sf  $ATOMIC/NIT/III/7jun99/niii_f_to_s_split.dat      NIII_F_TO_S
#
 ln -sf  $ATOMIC/NIT/IV/5dec96/nivphot_a.dat               PHOTNIV_A
 ln -sf  $ATOMIC/NIT/IV/5dec96/nivphot_b.dat               PHOTNIV_B
 ln -sf  $ATOMIC/NIT/IV/5dec96/nivosc_ns_split.dat         NIV_F_OSCDAT
 ln -sf  $ATOMIC/NIT/IV/5dec96/niv_f_to_s_split.dat        NIV_F_TO_S
 ln -sf  $ATOMIC/NIT/IV/5dec96/nivcol.dat                  NIV_COL_DATA
 ln -sf  $ATOMIC/NIT/IV/5dec96/nivdie.dat                  DIENIV
#
 ln -sf  $ATOMIC/NIT/V/5dec96/nvphot_a12.dat               PHOTNV_A
 ln -sf  $ATOMIC/NIT/V/5dec96/nvosc_a12_split.dat          NV_F_OSCDAT
 ln -sf  $ATOMIC/NIT/V/5dec96/nv_f_to_s_split_sm.dat       NV_F_TO_S
 ln -sf  $ATOMIC/NIT/V/5dec96/nvcol.dat                    NV_COL_DATA

 ln -sf  $ATOMIC/NIT/VI/4feb05/phot_smooth_50              PHOTNSIX_A
 ln -sf  $ATOMIC/NIT/VI/4feb05/fin_osc_pack                NSIX_F_OSCDAT
 ln -sf  $ATOMIC/NIT/VI/4feb05/f_to_s_pack_43              NSIX_F_TO_S
 ln -sf  $ATOMIC/NIT/VI/4feb05/col_guess.dat               NSIX_COL_DATA

 ln -sf  $ATOMIC/NIT/VII/4feb05/phot_smooth_50             PHOTNSEV_A
 ln -sf  $ATOMIC/NIT/VII/4feb05/fin_osc                    NSEV_F_OSCDAT
 ln -sf  $ATOMIC/NIT/VII/4feb05/col_guess.dat              NSEV_COL_DATA

#*****************************************************************************
#    Oxygen
#   --------
#*****************************************************************************
 ln -sf $ATOMIC/OXY/II/3oct00/o2osc_fin.dat                O2_F_OSCDAT
 ln -sf $ATOMIC/OXY/II/3oct00/f_to_s_ls.dat                O2_F_TO_S
 ln -sf $ATOMIC/OXY/II/3oct00/phot_sm_3000.dat             PHOTO2_A
 ln -sf $ATOMIC/OXY/II/3oct00/o2col.dat                    O2_COL_DATA
#
 ln -sf  $ATOMIC/OXY/III/5dec96/oiiiphot_a.dat             PHOTOIII_A
 ln -sf  $ATOMIC/OXY/III/5dec96/oiiiphot_b.dat             PHOTOIII_B
 ln -sf  $ATOMIC/OXY/III/5dec96/oiiiosc_split.dat          OIII_F_OSCDAT
 ln -sf  $ATOMIC/OXY/III/5dec96/oiii_f_to_s_split.dat      OIII_F_TO_S
 ln -sf  $ATOMIC/OXY/III/5dec96/oiiicol.dat                OIII_COL_DATA
 ln -sf  $ATOMIC/OXY/III/5dec96/oiiidie.dat                DIEOIII
#
 ln -sf  $ATOMIC/OXY/IV/5dec96/oivphot_a.dat               PHOTOIV_A
 ln -sf  $ATOMIC/OXY/IV/5dec96/oivphot_b.dat               PHOTOIV_B
 ln -sf  $ATOMIC/OXY/IV/5dec96/oivosc_ns_split.dat         OIV_F_OSCDAT
 ln -sf  $ATOMIC/OXY/IV/5dec96/oiv_f_to_s_split.dat        OIV_F_TO_S
 ln -sf  $ATOMIC/OXY/IV/5dec96/oivcol.dat                  OIV_COL_DATA
 ln -sf  $ATOMIC/OXY/IV/5dec96/oivdoubdie_103.dat          DIEOIV
#
 ln -sf  $ATOMIC/OXY/V/5dec96/ovphot_a.dat                 PHOTOV_A
 ln -sf  $ATOMIC/OXY/V/5dec96/ovphot_b.dat                 PHOTOV_B
 ln -sf  $ATOMIC/OXY/V/5dec96/ovosc_ns_split.dat           OV_F_OSCDAT
 ln -sf  $ATOMIC/OXY/V/5dec96/ov_f_to_s_split_ext.dat      OV_F_TO_S
 ln -sf  $ATOMIC/OXY/V/5dec96/ovcol.dat                    OV_COL_DATA
 ln -sf  $ATOMIC/OXY/V/5dec96/ovdie.dat                    DIEOV
#

 ln -sf  $ATOMIC/OXY/VI/5dec96/osixphot_a12.dat            PHOTOSIX_A
 ln -sf  $ATOMIC/OXY/VI/5dec96/osixosc_a12_split.dat       OSIX_F_OSCDAT
 ln -sf  $ATOMIC/OXY/VI/5dec96/osix_f_to_s_split.dat       OSIX_F_TO_S
 ln -sf  $ATOMIC/OXY/VI/5dec96/osixcol.dat                 OSIX_COL_DATA

 ln -sf  $atomic/OXY/VII/20sep02/fin_osc_pack              OSEV_F_OSCDAT
 ln -sf  $atomic/OXY/VII/20sep02/phot_sm_3000              PHOTOSEV_A
 ln -sf  $atomic/OXY/VII/20sep02/f_to_s_pack_43            OSEV_F_TO_S
 ln -sf  $atomic/OXY/VII/20sep02/col_data                  OSEV_COL_DATA

 ln -sf  $atomic/OXY/VIII/20sep02/fin_osc                   OVIII_F_OSCDAT
 ln -sf  $atomic/OXY/VIII/20sep02/phot_sm_3000             PHOTOVIII_A
# ln -sf  $atomic/OXY/VIII/20sep02/                        OVIII_F_TO_S
 ln -sf  $atomic/OXY/VIII/20sep02/col_guess.dat            OVIII_COL_DATA

#*****************************************************************************
# Neon
# ----
#*****************************************************************************
 ln -sf  $ATOMIC/NEON/II/1dec99/phot_sm_3000.dat           PHOTNe2_A
 ln -sf  $ATOMIC/NEON/II/1dec99/fin_osc.dat                Ne2_F_OSCDAT
 ln -sf  $ATOMIC/NEON/II/1dec99/f_to_s_42.dat              Ne2_F_TO_S
 ln -sf  $ATOMIC/NEON/II/1dec99/col_data.dat               Ne2_COL_DATA
#
 ln -sf  $ATOMIC/NEON/III/1dec99/phot_sm_3000.dat          PHOTNeIII_A
 ln -sf  $ATOMIC/NEON/III/1dec99/fin_osc.dat               NeIII_F_OSCDAT
 ln -sf  $ATOMIC/NEON/III/1dec99/f_to_s_40.dat             NeIII_F_TO_S
 ln -sf  $ATOMIC/NEON/III/1dec99/col_data.dat              NeIII_COL_DATA
#
 ln -sf  $ATOMIC/NEON/IV/1dec99/phot_sm_3000.dat           PHOTNeIV_A
 ln -sf  $ATOMIC/NEON/IV/1dec99/fin_osc.dat                NeIV_F_OSCDAT
 ln -sf  $ATOMIC/NEON/IV/1dec99/f_to_s_45.dat              NeIV_F_TO_S
 ln -sf  $ATOMIC/NEON/IV/1dec99/col_data.dat               NeIV_COL_DATA
#
 ln -sf  $ATOMIC/NEON/V/20jun01/phot_sm_3000.dat             PHOTNeV_A
 ln -sf  $ATOMIC/NEON/V/20jun01/nevosc_rev.dat               NeV_F_OSCDAT
 ln -sf  $ATOMIC/NEON/V/20jun01/f_to_s_60.dat                NeV_F_TO_S
 ln -sf  $ATOMIC/NEON/V/20jun01/col_data.dat                 NeV_COL_DATA

 ln -sf  $ATOMIC/NEON/VI/20jun01/phot_sm_3000.dat           PHOTNeSIX_A
 ln -sf  $ATOMIC/NEON/VI/20jun01/neviosc_rev.dat            NeSIX_F_OSCDAT
 ln -sf  $ATOMIC/NEON/VI/20jun01/f_to_s_43.dat              NeSIX_F_TO_S
 ln -sf  $ATOMIC/NEON/VI/20jun01/col_data.dat               NeSIX_COL_DATA

 ln -sf  $ATOMIC/NEON/VII/20jun01/phot_sm_3000.dat          PHOTNeSEV_A
 ln -sf  $ATOMIC/NEON/VII/20jun01/neviiosc_rev.dat          NeSEV_F_OSCDAT
 ln -sf  $ATOMIC/NEON/VII/20jun01/f_to_s_46.dat             NeSEV_F_TO_S
 ln -sf  $ATOMIC/NEON/VII/20jun01/col_data.dat              NeSEV_COL_DATA

 ln -sf  $ATOMIC/NEON/VIII/20jun01/phot_sm_3000.dat         PHOTNeVIII_A
 ln -sf  $ATOMIC/NEON/VIII/20jun01/neviiiosc_rev.dat        NeVIII_F_OSCDAT
 ln -sf  $ATOMIC/NEON/VIII/20jun01/f_to_s_29.dat            NeVIII_F_TO_S
 ln -sf  $ATOMIC/NEON/VIII/20jun01/col_guess.dat            NeVIII_COL_DATA

#*****************************************************************************
#   Aluminium
#   ---------
#*****************************************************************************
 ln -sf  $ATOMIC/AL/III/5aug97/aliii_phot_a.dat            PHOTAlIII_A
 ln -sf  $ATOMIC/AL/III/5aug97/aliii_osc_split.dat         AlIII_F_OSCDAT
 ln -sf  $ATOMIC/AL/III/5aug97/aliii_f_to_s_sm.dat         AlIII_F_TO_S
 ln -sf  $ATOMIC/AL/III/5aug97/col_guess.dat               AlIII_COL_DATA

#*****************************************************************************
#    Silicon
#   --------
#*****************************************************************************
 ln -sf  $ATOMIC/SIL/III/5dec96b/phot_op.dat                PHOTSkIII_A
 ln -sf  $ATOMIC/SIL/III/5dec96b/osc_op_split_rev.dat       SkIII_F_OSCDAT
 ln -sf  $ATOMIC/SIL/III/5dec96b/f_to_s_split.dat           SkIII_F_TO_S
 ln -sf  $ATOMIC/SIL/III/5dec96b/col_guess.dat              SkIII_COL_DATA
#
 ln -sf  $ATOMIC/SIL/IV/5dec96/phot_op.dat                 PHOTSkIV_A
 ln -sf  $ATOMIC/SIL/IV/5dec96/osc_op_split.dat            SkIV_F_OSCDAT
 ln -sf  $ATOMIC/SIL/IV/5dec96/f_to_s_split.dat            SkIV_F_TO_S
 ln -sf  $ATOMIC/SIL/IV/5dec96/col_guess.dat               SkIV_COL_DATA

#*****************************************************************************
# Phosphorous
# -------
#*****************************************************************************
# ln -sf  $ATOMIC/PHOS/III/15feb01/phot_data.dat             PHOTPIII_A
# ln -sf  $ATOMIC/PHOS/III/15feb01/osc_op.dat                PIII_F_OSCDAT
# ln -sf  $ATOMIC/PHOS/III/15feb01/f_to_s_36.dat             PIII_F_TO_S
# ln -sf  $ATOMIC/PHOS/III/15feb01/col_guess.dat             PIII_COL_DATA
#
 ln -sf  $ATOMIC/PHOS/IV/15feb01/phot_data_a.dat            PHOTPIV_A
 ln -sf  $ATOMIC/PHOS/IV/15feb01/phot_data_b.dat            PHOTPIV_B
 ln -sf  $ATOMIC/PHOS/IV/15feb01/pivosc_rev.dat             PIV_F_OSCDAT
 ln -sf  $ATOMIC/PHOS/IV/15feb01/f_to_s_36.dat              PIV_F_TO_S
 ln -sf  $ATOMIC/PHOS/IV/15feb01/col_guess.dat              PIV_COL_DATA
#
 ln -sf  $ATOMIC/PHOS/V/15feb01/phot_data.dat               PHOTPV_A
 ln -sf  $ATOMIC/PHOS/V/15feb01/pvosc_rev.dat               PV_F_OSCDAT
 ln -sf  $ATOMIC/PHOS/V/15feb01/f_to_s_16.dat               PV_F_TO_S
 ln -sf  $ATOMIC/PHOS/V/15feb01/col_guess.dat               PV_COL_DATA
#
# ln -sf  $ATOMIC/PHOS/VI/15feb01/phot_op.dat                PHOTPSIX_A
# ln -sf  $ATOMIC/PHOS/VI/15feb01/osc_op.dat                 PSIX_F_OSCDAT
# ln -sf  $ATOMIC/PHOS/VI/15feb01/f_to_s_25.dat              PSIX_F_TO_S
# ln -sf  $ATOMIC/PHOS/VI/15feb01/col_guess.dat              PSIX_COL_DATA

#*****************************************************************************
#   Sulpher
#   -------
#*****************************************************************************
 ln -sf  $ATOMIC/SUL/III/3oct00/phot_sm_3000.dat          PHOTSIII_A
 ln -sf  $ATOMIC/SUL/III/3oct00/siiiosc_fin.dat           SIII_F_OSCDAT
 ln -sf  $ATOMIC/SUL/III/3oct00/f_to_s_127.dat            SIII_F_TO_S
 ln -sf  $ATOMIC/SUL/III/3oct00/col_siii.dat              SIII_COL_DATA
#
 ln -sf  $ATOMIC/SUL/IV/3oct00/phot_sm_3000.dat           PHOTSIV_A
 ln -sf  $ATOMIC/SUL/IV/3oct00/sivosc_fin.dat             SIV_F_OSCDAT
 ln -sf  $ATOMIC/SUL/IV/3oct00/f_to_s_69.dat              SIV_F_TO_S
 ln -sf  $ATOMIC/SUL/IV/3oct00/col_siv.dat                SIV_COL_DATA
#
 ln -sf  $ATOMIC/SUL/V/3oct00/phot_sm_3000.dat            PHOTSV_A
 ln -sf  $ATOMIC/SUL/V/3oct00/svosc_fin.dat               SV_F_OSCDAT
 ln -sf  $ATOMIC/SUL/V/3oct00/f_to_s_50.dat               SV_F_TO_S
 ln -sf  $ATOMIC/SUL/V/3oct00/col_sv.dat                  SV_COL_DATA
#
 ln -sf  $ATOMIC/SUL/VI/3oct00/phot_sm_3000.dat           PHOTSSIX_A
 ln -sf  $ATOMIC/SUL/VI/3oct00/sviosc_fin.dat             SSIX_F_OSCDAT
 ln -sf  $ATOMIC/SUL/VI/3oct00/f_to_s_33.dat              SSIX_F_TO_S
 ln -sf  $ATOMIC/SUL/VI/3oct00/col_guess.dat              SSIX_COL_DATA

#*****************************************************************************
# Chlorine
# --------
#*****************************************************************************

 ln -sf  $ATOMIC/CHL/IV/15feb01/phot_data.dat              PHOTClIV_A
 ln -sf  $ATOMIC/CHL/IV/15feb01/f_to_s_58.dat              ClIV_F_TO_S
 ln -sf  $ATOMIC/CHL/IV/15feb01/clivosc_fin.dat            ClIV_F_OSCDAT
 ln -sf  $ATOMIC/CHL/IV/15feb01/col_data.dat               ClIV_COL_DATA

 ln -sf  $ATOMIC/CHL/V/15feb01/phot_data.dat               PHOTClV_A
 ln -sf  $ATOMIC/CHL/V/15feb01/f_to_s_41.dat               ClV_F_TO_S
 ln -sf  $ATOMIC/CHL/V/15feb01/clvosc_fin.dat              ClV_F_OSCDAT
 ln -sf  $ATOMIC/CHL/V/15feb01/col_data.dat                ClV_COL_DATA

 ln -sf  $ATOMIC/CHL/VI/15feb01/phot_data.dat              PHOTClSIX_A
 ln -sf  $ATOMIC/CHL/VI/15feb01/f_to_s_32.dat              ClSIX_F_TO_S
 ln -sf  $ATOMIC/CHL/VI/15feb01/clviosc_rev.dat            ClSIX_F_OSCDAT
 ln -sf  $ATOMIC/CHL/VI/15feb01/col_guess.dat              ClSIX_COL_DATA

 ln -sf  $ATOMIC/CHL/VII/15feb01/phot_data.dat             PHOTClSEV_A
 ln -sf  $ATOMIC/CHL/VII/15feb01/f_to_s_42.dat             ClSEV_F_TO_S
 ln -sf  $ATOMIC/CHL/VII/15feb01/clviiosc_rev.dat          ClSEV_F_OSCDAT
 ln -sf  $ATOMIC/CHL/VII/15feb01/col_guess.dat             ClSEV_COL_DATA

#*****************************************************************************
# Argon
# -----
#*****************************************************************************
 ln -sf  $ATOMIC/ARG/III/1dec99/phot_sm_3000.dat           PHOTArIII_A
 ln -sf  $ATOMIC/ARG/III/1dec99/fin_osc.dat                ArIII_F_OSCDAT
 ln -sf  $ATOMIC/ARG/III/1dec99/f_to_s_32.dat              ArIII_F_TO_S
 ln -sf  $ATOMIC/ARG/III/1dec99/col_data.dat               ArIII_COL_DATA
#
 ln -sf  $ATOMIC/ARG/IV/1dec99/phot_sm_3000.dat            PHOTArIV_A
 ln -sf  $ATOMIC/ARG/IV/1dec99/fin_osc.dat                 ArIV_F_OSCDAT
 ln -sf  $ATOMIC/ARG/IV/1dec99/f_to_s_50.dat               ArIV_F_TO_S
 ln -sf  $ATOMIC/ARG/IV/1dec99/col_data.dat                ArIV_COL_DATA
#
 ln -sf  $ATOMIC/ARG/V/1dec99/phot_sm_3000.dat             PHOTArV_A
 ln -sf  $ATOMIC/ARG/V/1dec99/fin_osc.dat                  ArV_F_OSCDAT
 ln -sf  $ATOMIC/ARG/V/1dec99/f_to_s_64.dat                ArV_F_TO_S
 ln -sf  $ATOMIC/ARG/V/1dec99/col_data.dat                 ArV_COL_DATA

 ln -sf  $ATOMIC/ARG/VI/15feb01/phot_sm_3000.dat           PHOTArSIX_A
 ln -sf  $ATOMIC/ARG/VI/15feb01/arviosc_rev.dat            ArSIX_F_OSCDAT
 ln -sf  $ATOMIC/ARG/VI/15feb01/f_to_s_30.dat              ArSIX_F_TO_S
 ln -sf  $ATOMIC/ARG/VI/15feb01/col_data.dat               ArSIX_COL_DATA

 ln -sf  $ATOMIC/ARG/VII/15feb01/phot_sm_3000.dat          PHOTArSEV_A
 ln -sf  $ATOMIC/ARG/VII/15feb01/arviiosc_rev.dat          ArSEV_F_OSCDAT
 ln -sf  $ATOMIC/ARG/VII/15feb01/f_to_s_46.dat             ArSEV_F_TO_S
 ln -sf  $ATOMIC/ARG/VII/15feb01/col_data.dat              ArSEV_COL_DATA

 ln -sf  $ATOMIC/ARG/VIII/15feb01/phot_sm_3000.dat         PHOTArVIII_A
 ln -sf  $ATOMIC/ARG/VIII/15feb01/arviiiosc_rev.dat        ArVIII_F_OSCDAT
 ln -sf  $ATOMIC/ARG/VIII/15feb01/f_to_s_33.dat            ArVIII_F_TO_S
 ln -sf  $ATOMIC/ARG/VIII/15feb01/col_guess.dat            ArVIII_COL_DATA

#*****************************************************************************
#Calcium
#-------
#*****************************************************************************

 ln -sf  $ATOMIC/CA/III/10apr99/phot_smooth.dat            PHOTCaIII_A
 ln -sf  $ATOMIC/CA/III/10apr99/osc_op_sp.dat              CaIII_F_OSCDAT
 ln -sf  $ATOMIC/CA/III/10apr99/f_to_s.dat                 CaIII_F_TO_S
 ln -sf  $ATOMIC/CA/III/10apr99/col_guess.dat              CaIII_COL_DATA
#
 ln -sf  $ATOMIC/CA/IV/10apr99/phot_smooth.dat             PHOTCaIV_A
 ln -sf  $ATOMIC/CA/IV/10apr99/osc_op_sp.dat               CaIV_F_OSCDAT
 ln -sf  $ATOMIC/CA/IV/10apr99/f_to_s.dat                  CaIV_F_TO_S
 ln -sf  $ATOMIC/CA/IV/10apr99/col_guess.dat               CaIV_COL_DATA
#
 ln -sf  $ATOMIC/CA/V/10apr99/phot_smooth.dat              PHOTCaV_A
 ln -sf  $ATOMIC/CA/V/10apr99/osc_op_sp.dat                CaV_F_OSCDAT
 ln -sf  $ATOMIC/CA/V/10apr99/f_to_s.dat                   CaV_F_TO_S
 ln -sf  $ATOMIC/CA/V/10apr99/col_guess.dat                CaV_COL_DATA

 ln -sf  $ATOMIC/CA/VI/10apr99/phot_smooth.dat             PHOTCaSIX_A
 ln -sf  $ATOMIC/CA/VI/10apr99/osc_op_sp.dat               CaSIX_F_OSCDAT
 ln -sf  $ATOMIC/CA/VI/10apr99/f_to_s.dat                  CaSIX_F_TO_S
 ln -sf  $ATOMIC/CA/VI/10apr99/col_guess.dat               CaSIX_COL_DATA

 ln -sf  $ATOMIC/CA/VII/10apr99/phot_smooth.dat            PHOTCaSEV_A
 ln -sf  $ATOMIC/CA/VII/10apr99/osc_op_sp.dat              CaSEV_F_OSCDAT
 ln -sf  $ATOMIC/CA/VII/10apr99/f_to_s.dat                 CaSEV_F_TO_S
 ln -sf  $ATOMIC/CA/VII/10apr99/col_guess.dat              CaSEV_COL_DATA

#*****************************************************************************
#    Iron
#   ------
#*****************************************************************************
# 
 ln -sf  $ATOMIC/FE/IV/18oct00/phot_sm_3000.dat            PHOTFeIV_A
# ln -sf $ATOMIC/FE/IV/18oct00/f_to_s_100.dat              FeIV_F_TO_S
 ln -sf  /home/lyrebird/jdh/cmfmod/10Lac/M44/FeIV_F_TO_S_new             FeIV_F_TO_S
 ln -sf  /home/lyrebird/jdh/cmfmod/10Lac/M44/FeIV_F_OSCDAT_new           FeIV_F_OSCDAT
# ln -sf  $ATOMIC/FE/IV/18oct00/feiv_osc.dat               FeIV_F_OSCDAT
 ln -sf  $ATOMIC/FE/IV/18oct00/col_data.dat                FeIV_COL_DATA
#
#ln -sf  $ATOMIC/FE/IV/23oct95/feiv_phot_sm_op.dat         PHOTFeIV_A
#ln -sf  $ATOMIC/FE/IV/23oct95/feiv_gf_kb.dat              FeIV_F_OSCDAT
#ln -sf  $ATOMIC/FE/IV/23oct95/feiv_f_to_s_67.dat          FeIV_F_TO_S
#ln -sf  $ATOMIC/FE/IV/23oct95/col_guess.dat               FeIV_COL_DATA
#
 ln -sf  $ATOMIC/FE/V/18oct00/phot_sm_3000.dat              PHOTFeV_A
 ln -sf  $ATOMIC/FE/V/18oct00/f_to_s_139.dat                FeV_F_TO_S
 ln -sf  $ATOMIC/FE/V/18oct00/fev_osc.dat                   FeV_F_OSCDAT
 ln -sf  $ATOMIC/FE/V/18oct00/col_guess.dat                 FeV_COL_DATA
#
#ln -sf  $ATOMIC/FE/V/23oct95/fev_phot_sm_op.dat           PHOTFeV_A
#ln -sf  $ATOMIC/FE/V/23oct95/fev_f_to_s_46.dat            FeV_F_TO_S
#ln -sf  $ATOMIC/FE/V/23oct95/fev_gf_kb.dat                FeV_F_OSCDAT
#ln -sf  $ATOMIC/FE/V/23oct95/col_guess.dat                FeV_COL_DATA
#
#
 ln -sf  $ATOMIC/FE/VI/18oct00/phot_sm_3000.dat            PHOTFeSIX_A
 ln -sf  $ATOMIC/FE/VI/18oct00/f_to_s_67.dat               FeSIX_F_TO_S
 ln -sf  $ATOMIC/FE/VI/18oct00/fevi_osc.dat                FeSIX_F_OSCDAT
 ln -sf  $ATOMIC/FE/VI/18oct00/col_data.dat                FeSIX_COL_DATA
#
 ln -sf  $ATOMIC/FE/VII/18oct00/phot_sm_3000.dat          PHOTFeSEV_A
 ln -sf  $ATOMIC/FE/VII/18oct00/f_to_s_69.dat             FeSEV_F_TO_S
 ln -sf  $ATOMIC/FE/VII/18oct00/fevii_osc.dat             FeSEV_F_OSCDAT
 ln -sf  $ATOMIC/FE/VII/18oct00/col_guess.dat             FeSEV_COL_DATA

 ln -sf  $ATOMIC/FE/VIII/8may97/feviii_phot_op.dat        PHOTFeVIII_A
 ln -sf  $ATOMIC/FE/VIII/8may97/feviii_f_to_s_53.dat      FeVIII_F_TO_S
 ln -sf  $ATOMIC/FE/VIII/8may97/feviii_osc_kb_rk.dat      FeVIII_F_OSCDAT
 ln -sf  $ATOMIC/FE/VIII/8may97/col_guess.dat             FeVIII_COL_DATA

#*****************************************************************************
#  Nickel
#  ------
#*****************************************************************************

 ln -sf  $ATOMIC/NICK/IV/18oct00/phot_data.dat                PHOTNkIV_A
# ln -sf  $ATOMIC/NICK/IV/18oct00/f_to_s_115.dat              NkIV_F_TO_S
 ln -sf  /home/lyrebird/jdh/cmfmod/10Lac/M44/NkIV_F_TO_S_new                NkIV_F_TO_S
 ln -sf  $ATOMIC/NICK/IV/18oct00/nkiv_osc.dat                 NkIV_F_OSCDAT
 ln -sf  $ATOMIC/NICK/IV/18oct00/col_guess.dat                NkIV_COL_DATA

 ln -sf  $ATOMIC/NICK/V/18oct00/phot_data.dat                 PHOTNkV_A
 ln -sf  $ATOMIC/NICK/V/18oct00/f_to_s_152.dat                NkV_F_TO_S
 ln -sf  $ATOMIC/NICK/V/18oct00/nkv_osc.dat                   NkV_F_OSCDAT
 ln -sf  $ATOMIC/NICK/V/18oct00/col_guess.dat                 NkV_COL_DATA

 ln -sf  $ATOMIC/NICK/VI/18oct00/phot_data.dat                PHOTNkSIX_A
 ln -sf  $ATOMIC/NICK/VI/18oct00/f_to_s_62.dat                NkSIX_F_TO_S
 ln -sf  $ATOMIC/NICK/VI/18oct00/nkvi_osc.dat                 NkSIX_F_OSCDAT
 ln -sf  $ATOMIC/NICK/VI/18oct00/col_guess.dat                NkSIX_COL_DATA

 ln -sf  $ATOMIC/NICK/VII/18oct00/phot_data.dat               PHOTNkSEV_A
 ln -sf  $ATOMIC/NICK/VII/18oct00/f_to_s_61.dat               NkSEV_F_TO_S
 ln -sf  $ATOMIC/NICK/VII/18oct00/nkvii_osc.dat               NkSEV_F_OSCDAT
 ln -sf  $ATOMIC/NICK/VII/18oct00/col_guess.dat               NkSEV_COL_DATA

 ln -sf  $ATOMIC/NICK/VIII/11jun01/phot_data.dat              PHOTNkVIII_A
 ln -sf  $ATOMIC/NICK/VIII/11jun01/f_to_s_48.dat              NkVIII_F_TO_S
 ln -sf  $ATOMIC/NICK/VIII/11jun01/nkviii_osc.dat             NkVIII_F_OSCDAT
 ln -sf  $ATOMIC/NICK/VIII/11jun01/col_guess.dat              NkVIII_COL_DATA

 ln -sf  $ATOMIC/NICK/IX/11jun01/phot_data.dat                PHOTNkIX_A
 ln -sf  $ATOMIC/NICK/IX/11jun01/f_to_s_48.dat                NkIX_F_TO_S
 ln -sf  $ATOMIC/NICK/IX/11jun01/nkix_osc.dat                 NkIX_F_OSCDAT
 ln -sf  $ATOMIC/NICK/IX/11jun01/col_guess.dat                NkIX_COL_DATA
#
#*****************************************************************************
#
# END OF ATOMIC DATA SOFT LINKS
#
#*****************************************************************************

#If we pass a paremeter to the BATCH file, we assume that we just want to
#ln -sf  data files. In this case we exit. 

if ( $1 == '') then

else
    echo "Data files assigned"
    exit
endif

#
#***********************************************************************
#    Input files: Only needed is different from defaults
#***********************************************************************

# ln -sf  $MODEL/vadat.dat                                  VADAT
# ln -sf  $MODEL/input.dat                                  IN_ITS
# ln -sf  $MODEL/cfdat_in.dat                               CFDAT

#***********************************************************************
#    Departure coeficient input files. T_IN is required for a new
#    model with GRID=FALSE. Other links not needed unless diferent
#    from default.
#***********************************************************************


 ln -sf  He2_IN                 T_IN

#
#***********************************************************************
#    Output files: If disk space is at a preimu, these two files might
#     be better stored on an alternative disk.
#***********************************************************************

#ln -sf /usr/limey/jdh/BAMAT            BAMAT
#ln -sf /usr/limey/jdh/BAION            BAION

#***********************************************************************
#***********************************************************************
#    Run program
#***********************************************************************
#***********************************************************************

#Putting time stamp etc on program

rm -f batch.log

echo "Program started on:" > 'batch.log'
date >> 'batch.log'
echo "Machine name is :" >> 'batch.log'
uname -n >> 'batch.log'
(exec nice $CMFGEN_PROG  >>& 'batch.log')&
wait
#
echo "Program finished on:" >> 'batch.log'
date >> 'batch.log'


#***********************************************************************
# Execute the command in next_job.sh. This will allow another job to
# be started.
#***********************************************************************

# $NEXT_JOB &
