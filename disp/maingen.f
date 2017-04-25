!
! Main subroutine to examine model output from CMFGEN. Routine is called by
! DISPGEN. All Model and Atomic data is primarily read in by DISPGEN, and
! is contained in the module MOD_DISP
!
	SUBROUTINE MAINGEN(RMDOT,LUM,
	1                    ND,NP,NC,
	1                    N_MAX,ND_MAX,NC_MAX,NP_MAX,
	1                    N_LINE_MAX,N_PLT_MAX)
	USE MOD_DISP
	USE MOD_USR_OPTION
	USE MOD_USR_HIDDEN
	USE MOD_WR_STRING
	USE MOD_LEV_DIS_BLK
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
! Altered  01-MAr-2016 : Adjusted PLNID -- only compute continuum opacity every 1009 km/s.
!                          (Results in considerable speed up [1-Feb-2016]).
! Altered  09-Jun-2015 : Added PLT_PRFS to call.
!                          Altered DTDP option.
! Altered  07-Jul-2011 : Included commands to look at f, and Mdot in "shell" models.
! Altered  15-Mar-2011 :  Section for plotting photoionization cross-sections
!                           removed to subroutine. Can now plot all ground-state
!                           cross section for a species.
!                         TCMF option installed to use GREY T from CMFGEN calculation.
! Altered  13-Jun-2010 :  Rayleigh scattering included in OPACITIES.INC.
!                         ESEC replaced by CHI_SCAT in may locations.
!                         Still some issues to be decided (i.e., Raleigh scattering to opacity file?)
!  	                  Changed to SOBEW_GRAD_V2 -- changed call. Adapted to work for very low opacities.
!                         Now include variable INNER_BND_METH.
!
! Altered  15-Aug-2003 :  QF option installed. Allows ion column densities to be compared.
!                         Option added to PLTPHOT to allow photoioization cross-section
!                           to be compted at a particular electron desnity. This allows the
!                           effects of level dissolution to be seen.
! Altered  28-Mar-2003 :  MAX_ION replaced bu NUM_IONS everywhere.
! Altered  01-Jub-1998 :  MOD_LEV_DIS_BLK replaces  include file.
!                         REXT_COEF changed to REXT_COEF_V2
! altered  12/16/96  DLM  Changed USRINP call to USR_OPTION.  USR_OPTION
!                         is a generic subroutine call which is defined in
!                         MOD_USR_OPTION.  The interface in MOD_USR_OPTION
!                         determines which usr_option subroutine to call by 
!                         examining the types of variables passed to USR_OPTION.
!
!                         .sve files created by SVE_FILE routines
!
!                         Also cleaned up a few places:
!                           1) Changed all comparisons to input value X to compare 
!                              to X(1:?) => needed by new .sve routines.
!                           2) Changed LOGR and LINR options to XLOGR and XLINR
!                              to be consistant => needed by new .sve routines.
!
! altered   2/27/97  DLM  Added box= option to create .box files
!                         Added hidden options for reading of .sve file
!
! 
!
	INTEGER NC,NP,ND
	INTEGER N_MAX,ND_MAX,NC_MAX,NP_MAX
	INTEGER N_LINE_MAX,N_PLT_MAX
!
	REAL*8 LUM,RMDOT
!
	INTEGER, PARAMETER :: NLF_MAX=1001
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	DOUBLE PRECISION CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
!
! 
!
	REAL*8 EINA
	REAL*8 OSCIL
	REAL*8 GUPDIE
	REAL*8 GLDGU
	REAL*8 EDGEDIE
	REAL*8 GLOW
	REAL*8 NUST(ND)
! 
!
! Arrays containing LINE frequencies in numerical order.
!
	REAL*8 VEC_FREQ(N_LINE_MAX)
	REAL*8 VEC_OSCIL(N_LINE_MAX)
	REAL*8 VEC_EINA(N_LINE_MAX)
	REAL*8 VEC_DP_WRK(N_LINE_MAX)
	INTEGER VEC_INDX(N_LINE_MAX)
	INTEGER VEC_ION_INDX(N_LINE_MAX)
	INTEGER VEC_MNL_F(N_LINE_MAX)
	INTEGER VEC_MNUP_F(N_LINE_MAX)
	INTEGER VEC_INT_WRK(N_LINE_MAX)
	CHARACTER*6 VEC_SPEC(N_LINE_MAX)
	CHARACTER*6 VEC_CHAR_WRK(N_LINE_MAX)
	CHARACTER*6 VEC_TRANS_TYPE(N_LINE_MAX)
	CHARACTER*80, ALLOCATABLE :: VEC_TRANS_NAME(:)
	CHARACTER*80, ALLOCATABLE :: TMP_TRANS_NAME(:)
	INTEGER N_LINE_FREQ
	INTEGER MNL_F
	INTEGER MNUP_F
!
! Angular quadrature weights.
!
	REAL*8 P(NP)				!Impact parameters
	REAL*8 JQW(ND,NP),HMIDQW(ND,NP)
	REAL*8 KQW(ND,NP),NMIDQW(ND,NP)
	REAL*8 HQW(ND,NP)
!
! Note that PFDOP contains the lineprofile in doppler units, whilst
! PF is use for the line profile in frequency units.
!
	REAL*8 LINE_PRO(NLF_MAX),LFQW(NLF_MAX),PF(NLF_MAX)
	REAL*8 ERF(NLF_MAX),PFDOP(NLF_MAX)
	REAL*8 JNU(ND,NLF_MAX+1),HNU(ND,NLF_MAX+1)
	REAL*8 FEDD(ND,NLF_MAX+1),GEDD(ND,NLF_MAX+1)
	REAL*8 HBC_VEC(3,NLF_MAX+1),NBC_VEC(3,NLF_MAX+1)
	REAL*8 INBC_VEC(NLF_MAX+1)
	REAL*8 TDOP,VTURB
!
	REAL*8 JBAR(ND),ZNET(ND)
	REAL*8 CHIROSS(ND),TAUROSS(ND)
!
! Variables required to compute TGREY and for interpolations.
!
	REAL*8 JQWEXT(NP_MAX,NP_MAX),KQWEXT(NP_MAX,NP_MAX),PEXT(NP_MAX)
	REAL*8 TA(NP_MAX),TB(NP_MAX),TC(NP_MAX)
	REAL*8 Z(NP_MAX),DTAU(NP_MAX),XM(NP_MAX),RJ(NP_MAX)
	REAL*8 CHI(NP_MAX),REXT(NP_MAX),dCHIdr(NP_MAX)
	REAL*8 INBC,HBC,HBCNEW,NBC,FA(NP_MAX),GAM(NP_MAX),GAMH(NP_MAX)

	REAL*8 RJEXT(NP_MAX),FEXT(NP_MAX),Q(NP_MAX),FOLD(NP_MAX)
	REAL*8 CHIEXT(NP_MAX),ETAEXT(NP_MAX),ESECEXT(NP_MAX)
	REAL*8 SOURCEEXT(NP_MAX)
	REAL*8 ZETAEXT(NP_MAX),THETAEXT(NP_MAX)
	REAL*8 COEF(0:3,NP_MAX)
	REAL*8 TGREY(NP_MAX)
	INTEGER GRID(NP_MAX),INDX(NP_MAX)
	LOGICAL INACCURATE,REXT_COMPUTED,GREY_COMP
	LOGICAL GREY_WITH_V_TERMS
	LOGICAL GREY_REL
	LOGICAL PLANE_PARALLEL_NOV
	LOGICAL FLUSH_FILE
!
	REAL*8 XV(N_PLT_MAX),XNU(N_PLT_MAX)
	REAL*8 YV(N_PLT_MAX),ZV(N_PLT_MAX),WV(N_PLT_MAX)
!
! 
!
! Collisional matrices.
!
	REAL*8 OMEGA_F(N_MAX,N_MAX)
	REAL*8 dln_OMEGA_F_dlnT(N_MAX,N_MAX)
	REAL*8 OMEGA_S(N_MAX,N_MAX)
	REAL*8 dln_OMEGA_S_dlnT(N_MAX,N_MAX)
	EXTERNAL OMEGA_GEN_V3
!
	REAL*8 UNIT_VEC(N_MAX)
	REAL*8 ZERO_VEC(N_MAX)
!
	REAL*8 DOP_PRO
	REAL*8 S15ADF,XCROSS_V2
	EXTERNAL JTRPWGT,HTRPWGT,KTRPWGT,NTRPWGT
	EXTERNAL JWEIGHT,HWEIGHT,KWEIGHT,NWEIGHT,XCROSS
!
! Photoionization cross-section routines.
!
	EXTERNAL SUB_PHOT_GEN
!
	CHARACTER*30 UC
	EXTERNAL UC
	REAL*8 KEV_TO_HZ,ANG_TO_HZ
	REAL*8 PI
	REAL*8 C_CMS
	REAL*8 C_KMS
!
	CHARACTER(LEN=6) METHOD
	CHARACTER(LEN=10) TYPE_ATM
	CHARACTER(LEN=10) INNER_BND_METH
!
	REAL*8, ALLOCATABLE :: CHI_PAR(:,:)
	REAL*8, ALLOCATABLE :: ETA_PAR(:,:)
	REAL*8 YMAPV(MAX_ION)			!Defined so use instead of NUM_IONS
	REAL*8 XMAPV(ND)
	LOGICAL LAST_NON_ZERO
	CHARACTER*10 LOC_ION_ID
!
	REAL*8, ALLOCATABLE :: CHI_LAM(:,:)
	REAL*8, ALLOCATABLE :: ETA_LAM(:,:)
	REAL*8, ALLOCATABLE :: CHI_TOT_LAM(:)
	REAL*8, ALLOCATABLE :: ETA_TOT_LAM(:)
	REAL*8, ALLOCATABLE :: LAM_VEC(:)
	INTEGER NLAM
	INTEGER NION
!
! Local variables
!
	REAL*8 FB(ND,ND),WM(ND,ND)
	REAL*8 GB(ND),U(ND),VB(ND),VC(ND),CHIL(ND),ETAL(ND)
	REAL*8 SOURCE(ND),TCHI(ND),ZETA(ND),THETA(ND)
	REAL*8 ETA(ND),ETA_WITH_ES(ND)
	REAL*8 ESEC(ND),EMHNUKT(ND),VT(ND)
	REAL*8 CHI_SCAT(ND),CHI_RAY(ND)
	REAL*8 AV(ND),CV(ND)
	REAL*8 FORCE_MULT(ND_MAX)
	LOGICAL DO_DPTH(ND)
!
! Used for computing distrbution of lines with respect to the Sobolev Line
! optical depth.
!
	INTEGER NBINS
	REAL*8 DELTA_TAU
	REAL*8 TAU_BEG
	REAL*8 TAU_MIN
	REAL*8 TAU_CONSTANT
	LOGICAL WEIGHT_NV
	LOGICAL MEAN_TAU
	LOGICAL RADIAL_TAU
!
	INTEGER I,J,K,L
	INTEGER IMIN,IMAX
	INTEGER ISPEC,ID
	INTEGER NLF,ML
	INTEGER IOS,NFREQ,R_INDX
	INTEGER NL,NUP
	INTEGER NEXT_LOC
	INTEGER PHOT_ID
	INTEGER DPTH_INDX
	INTEGER LINE_INDX
!
	INTEGER, PARAMETER :: ITAU_GRT_LIM=6
	INTEGER TAU_GRT_LOGX(-ITAU_GRT_LIM:ITAU_GRT_LIM)
!
	REAL*8 DTDR,DBB,S1,IC
	REAL*8 HFLUX
	REAL*8 EXC_EN,EDGE_FREQ
	REAL*8 RVAL,TAU_VAL,ED_VAL
	REAL*8 XDIS,YDIS,DIS_CONST
	REAL*8 LAM_ST,LAM_EN,DEL_NU
	REAL*8 NU_ST,NU_EN
	REAL*8 FREQ_RES,FREQ_MAX
	REAL*8 T1,T2,T3,T4
	REAL*8 EK_EJECTA
	REAL*8 TMP_ED
	REAL*8 TAU_LIM
	REAL*8 DEPTH_LIM
	REAL*8 TEMP,TSTAR,NEW_RSTAR,NEW_VSTAR
	REAL*8 VSM_DIE_KMS
	REAL*8 DIST_KPC
!
	INTEGER NPINS,NCX,NDX,NPX
	INTEGER ND_TMP
	REAL*8 FREQ,FL,FL_SAVE
	REAL*8 TAU_SOB
	REAL*8  TMP_GION
	EQUIVALENCE (FREQ,FL)
!
!
	REAL*8 SCLHT,VCORE,VPHOT,VINF1,V_BETA1,V_EPS1
	REAL*8 VINF2,V_BETA2,V_EPS2
!
! Storage for constants used for evaluating the level disolution.
! Required by SUP_TO_FULL and LTE_POP_WLD.
!
	LOGICAL LEVEL_DISSOLUTION
!
	CHARACTER*80 NAME,XAXIS,YAXIS,XAXSAV
!
	INTEGER, PARAMETER :: NSC=31
	COMMON/TOPBORD/ SCED(NSC),XED(NSC),NXED,TOPLABEL
	REAL*8 SCED,XED
	INTEGER NXED
	CHARACTER(LEN=30) TOPLABEL
	DATA SCED/2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,
	1         10.5,11.0,11.5,12.0,12.5,13.0,13.5,14.0,14.5,15,15.5,
	1         16,16.5,17.0/
!
	INTEGER LEV(10)
	INTEGER IDEPTH(10)
	LOGICAL FLAG,LINV,TRAPFORJ,JONS,JONLY,IN_R_SUN,SONLY
	LOGICAL ELEC,DIF,SCALE,THICK,NORAD,ROSS,INC_RAY_SCAT
	LOGICAL SPEC_FRAC,RADIAL
	LOGICAL LST_DEPTH_ONLY
	LOGICAL DIE_REG,DIE_WI
	LOGICAL VALID_VALUE
	LOGICAL FOUND
	LOGICAL NEW_FORMAT,NEW_FILE
	LOGICAL NEW_FREQ
	LOGICAL DO_TAU
	LOGICAL DO_KAP
	LOGICAL DONE_LINE
	LOGICAL LINE_STRENGTH
	LOGICAL PLT_J,PLT_H,PLT_LF,PLT_FM
!
	REAL*8 RED_EXT
	REAL*8 VDOP_FG_FRAC
	REAL*8 VDOP_MOM_FRAC
	CHARACTER*30 TWO_PHOT_OPTION
	CHARACTER*30 FREQ_INPUT
	CHARACTER*10 SOL_OPT
	CHARACTER*10 FG_SOL_OPT
	CHARACTER*10 N_TYPE
	LOGICAL KEV_INPUT,HZ_INPUT,ANG_INPUT
	LOGICAL LINE_BL,FULL_ES,EDDC,HAM,SKIPEW
	LOGICAL THK_CONT,THK_LINE,LIN_DC,LINX,LINY
	LOGICAL FIRST_RATE
	LOGICAL L_TRUE,L_FALSE
	DATA L_TRUE/.TRUE./
	DATA L_FALSE/.FALSE./
!
	LOGICAL XRAYS
	REAL*8 FILL_FAC_XRAYS,T_SHOCK,V_SHOCK
	REAL*8 FILL_VEC_SQ(ND)
!
! Equivalent width and flux variables.
!
	INTEGER CNT
	REAL*8 ERR(ND),EWACC
	REAL*8 EW,EWOLD,MOMEW,CONT_INT,MOMCONT_INT
!
	REAL*8 AMASS
!
	CHARACTER MAIN_OPT_STR*80
	CHARACTER X*20
	CHARACTER XOPT*10
	CHARACTER XSPEC*10
!
	CHARACTER TYPE*10
	CHARACTER(LEN=80) STRING
	CHARACTER(LEN=10) COL_OPT
	CHARACTER(LEN=80) FILENAME
!
	CHARACTER FMT*80
	CHARACTER*120 DEFAULT
	CHARACTER*120 DESCRIPTION
	CHARACTER(LEN=1) FIRST_COLR
!
! External functions
!
	EXTERNAL WR_PWD
	CHARACTER*120 WR_PWD
!
	REAL*8 GFF,GBF,LAMVACAIR,SPEED_OF_LIGHT,BOLTZMANN_CONSTANT
	REAL*8 FUN_PI,SECS_IN_YEAR,MASS_SUN,LUM_SUN,ATOMIC_MASS_UNIT
	REAL*8 ASTRONOMICAL_UNIT,RAD_SUN,TEFF_SUN,STEFAN_BOLTZ
	INTEGER GET_INDX_DP
	EXTERNAL GFF,GBF,LAMVACAIR,SPEED_OF_LIGHT,GET_INDX_DP
	EXTERNAL FUN_PI,SECS_IN_YEAR,MASS_SUN,LUM_SUN,ATOMIC_MASS_UNIT
	EXTERNAL ASTRONOMICAL_UNIT,RAD_SUN,TEFF_SUN,BOLTZMANN_CONSTANT,STEFAN_BOLTZ
!
	LOGICAL PRESENT			!indicates whether species is linked
	LOGICAL EQUAL
!
	INTEGER, PARAMETER :: LU_IN=10		!File input (open/close)
	INTEGER, PARAMETER :: LU_OUT=11		!File input (open/close)
	INTEGER, PARAMETER :: LU_LOG=7
	INTEGER, PARAMETER :: LU_CROSS=15
	INTEGER, PARAMETER :: LU_NET=17
	INTEGER, PARAMETER :: LU_REC=18
	INTEGER, PARAMETER :: LU_COL=20		!Collisonal
	INTEGER, PARAMETER :: T_IN=5                 !Terminal input
	INTEGER, PARAMETER :: T_OUT=6                !Terminal output
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
	INTEGER, PARAMETER :: ITEN=10
	INTEGER, PARAMETER :: ITHREE=3
!
	REAL*8, PARAMETER :: RZERO=0.0D0
	REAL*8, PARAMETER :: RONE=1.0D0
	REAL*8, PARAMETER :: RTWO=2.0D0
	REAL*8, PARAMETER :: RTHREE=3.0D0
!
! 
!
	XAXIS=' '
	XAXSAV=' '
	YAXIS=' '
	METHOD='LOGLOG'
	FIRST_COLR='T'
	DPTH_INDX=ND/2
	COL_OPT='CL'
!
! NB: ZERO_FLUX condition is equivalent to DIF=.TRUE. and DBB=0.
!     Thus, for routine compatibility, we set DBB=0 when DIF is false.
!
	IF(V(1) .GT. 10000.0D0)THEN
	  INNER_BND_METH='ZERO_FLUX'		!Default for SN model
	  DIF=.FALSE.
	ELSE
	  INNER_BND_METH='DIFFUSION'		!Default for normal stars
	  DIF=.TRUE.
	END IF
	TRAPFORJ=.TRUE.
	TSTAR=T(ND)
!
	TYPE_ATM='P'
	T1=LOG(MASS_DENSITY(5)/MASS_DENSITY(1))/LOG(R(1)/R(5))
	WRITE(6,'(A)')RED_PEN
	IF(T1 .LT. 2)THEN
	  WRITE(6,'(A)')' Power for density law variation at outer boudary is shallower than r^{-n} where n=2'
	  WRITE(6,'(A)')' Shown below is the power based on depths 1 and I'
	  WRITE(6,*)' '
	  DO I=2,MIN(15,ND)
	    WRITE(6,'(X,I4,ES12.3)')I,LOG(MASS_DENSITY(I)/MASS_DENSITY(1))/LOG(R(1)/R(I))
	  END DO
	  CALL USR_OPTION(T1,'POW','2','Power for rhj propto r^{-n}')
	  WRITE(6,*)' '
	END IF
!
	IF(T1 .LT. 100)THEN
	  WRITE(TYPE_ATM(2:6),'(F5.2)')T1
	  WRITE(6,'(A)')' For optical depth calculations we will assume a power density distribution'
	  WRITE(6,'(A,F5.2)')' at the outer boundary. Density exponent for atmosphere is ',T1
	ELSE
	  WRITE(TYPE_ATM(2:9),'(ES8.2)')T1
	  WRITE(6,'(A)')' For optical depth calculations we will assume a power density distribution'
	  WRITE(6,'(A,ES8.2)')' at the outer boundary. Density exponent for atmosphere is ',T1
	END IF
	WRITE(6,'(A)')DEF_PEN
!
! If the grey temperature is available in CMFGEN, we don't need to (although we can with the
! GREY option) compute it.
!
	GREY_COMP=.FALSE.
	IF(CMFGEN_TGREY(ND) .EQ. 0.0D0)THEN
	  GREY_COMP=.FALSE.
	ELSE
	  GREY_COMP=.TRUE.
	  TGREY(1:ND)=CMFGEN_TGREY(1:ND)
	END IF
	FIRST_RATE=.TRUE.
!
	XRAYS=.TRUE.                    !Changed def .FALSE.
	FILL_FAC_XRAYS=1.0D-100
	T_SHOCK=300.0D0
	V_SHOCK=200.0D0
	FILL_VEC_SQ(:)=0.0D0
!
	LEVEL_DISSOLUTION=.TRUE.
	PI=FUN_PI()
	DIST_KPC=1.0D0
	VSM_DIE_KMS=3000.0D0
	REXT_COMPUTED=.FALSE.
!
	LST_DEPTH_ONLY=.FALSE.
!
	LAM_ST=0.0D0
	LAM_EN=1.0D+06
!
	ZERO_VEC(:)=0.0D0
	UNIT_VEC(:)=1.0D0
!
! Conversion factor from Kev to units of 10^15 Hz.
! Conversion factor from Angstroms to units of 10^15 Hz.
!
	KEV_TO_HZ=0.241838E+03
	ANG_TO_HZ=SPEED_OF_LIGHT()*1.0D-07  	!10^8/10^15
	C_CMS=SPEED_OF_LIGHT()
	C_KMS=1.0D-05*SPEED_OF_LIGHT()
!
	ANG_INPUT=.TRUE.
	HZ_INPUT=.FALSE.
	KEV_INPUT=.FALSE.
	FREQ_INPUT='Units are Angstroms'
!
! DTDR is used in the diffusion approximation inner boundary condition.
! In the MAIN code it is fixxd by the ROSSELAND mean opacity, and the
! LUMINOSITY.
!
	DTDR=(T(ND)-T(ND-1))/(R(ND-1)-R(ND))
!
	THK_CONT=.TRUE.
	THK_LINE=.TRUE.
	LINE_BL=.TRUE.
!
! Take log of R (default X axis).
!
	DO I=1,ND
	  XV(I)=DLOG10(R(I)/R(ND))
	END DO
	XAXIS='Log(r/R\d*\u)'
	XAXSAV=XAXIS
	DEFAULT=WR_PWD()
	CALL RD_STRING(NAME,'Title',DEFAULT,
	1    'Title to be placed on all graphs')
!
! Allocate arrays:
!
	IOS=0
	ALLOCATE (CHI_PAR(ND,NUM_IONS),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (ETA_PAR(ND,NUM_IONS),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(T_OUT,*)'Error in MAINGEN'
	  WRITE(T_OUT,*)'Unable to allocate CHI_PAR or  ETA_PAR'
	  STOP
	END IF
!
! 
!
! Determine the population of the ground state, summing up over J levels that
! are within 2000 cm^-1 of each other. These are used to evaluate the
! free-free rates, and also provide a check on effective recombination
! rates.
!
	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES)THEN
	    ATM(ID)%DXzV(1:ND)=ATM(ID)%DXzV_F(1:ND)
	    IF(ATM(ID+1)%XzV_PRES)THEN
	      DO J=2,ATM(ID+1)%NXzV_F
	        IF((ATM(ID+1)%EDGEXzV_F(J)-ATM(ID+1)%EDGEXzV_F(1))
	1          *1.0D+015/C_CMS .LT. 2000)THEN
	          DO I=1,ND
	            ATM(ID)%DXzV(I)=ATM(ID)%DXzV(I)+ATM(ID+1)%XzV_F(J,I)
	          END DO
	        END IF
	      END DO
	    END IF
	  END IF
	END DO
!
! 
!
! Determine impact parameters. Calculate all required Angular quadrature
! weights.
!
	CALL IMPAR(P,R,R(ND),NC,ND,NP)
!
! The HMIDQW and NMIDQW qw are defined at the MIDPOINTS of the mesh, hence
! the option TRUE.
!
	IF(TRAPFORJ)THEN
	  CALL GENANGQW(JQW,R,P,WM(1,1),WM(1,3),WM(1,5)
	1         ,NC,ND,NP,JTRPWGT,.FALSE.)
	  CALL GENANGQW(HQW,R,P,WM(1,1),WM(1,3),WM(1,5)
	1        ,NC,ND,NP,HTRPWGT,.FALSE.)
	  CALL GENANGQW(KQW,R,P,WM(1,1),WM(1,3),WM(1,5)
	1        ,NC,ND,NP,KTRPWGT,.FALSE.)
	  CALL GENANGQW(HMIDQW,R,P,WM(1,1),WM(1,3),WM(1,5)
	1        ,NC,ND,NP,HTRPWGT,.TRUE.)
	  CALL GENANGQW(NMIDQW,R,P,WM(1,1),WM(1,3),WM(1,5)
	1        ,NC,ND,NP,NTRPWGT,.TRUE.)
	ELSE
	  CALL GENANGQW(JQW,R,P,WM(1,1),WM(1,3),WM(1,5)
	1         ,NC,ND,NP,JWEIGHT,.FALSE.)
	  CALL GENANGQW(HQW,R,P,WM(1,1),WM(1,3),WM(1,5)
	1        ,NC,ND,NP,HWEIGHT,.FALSE.)
	  CALL GENANGQW(KQW,R,P,WM(1,1),WM(1,3),WM(1,5)
	1        ,NC,ND,NP,KWEIGHT,.FALSE.)
	  CALL GENANGQW(HMIDQW,R,P,WM(1,1),WM(1,3),WM(1,5)
	1        ,NC,ND,NP,HWEIGHT,.TRUE.)
	  CALL GENANGQW(NMIDQW,R,P,WM(1,1),WM(1,3),WM(1,5)
	1        ,NC,ND,NP,NWEIGHT,.TRUE.)
	END IF
!
! 
! Compute vector constants for evaluating the level disolution. These
! constants are the same for all species. These are stored in a common block,
! and are required by SUP_TO_FULL and LTE_POP_WLD. NB: POPION was read in from
! the RVTJ file.
!
	CALL COMP_LEV_DIS_BLK(ED,POPION,T,LEVEL_DISSOLUTION,ND)
!
! Set 2-photon data with current atomic models and populations.
!
        DO ID=1,NUM_IONS-1
          I=ID
          CALL SET_TWO_PHOT_DISP_V3(ION_ID(ID), ID,
	1       ATM(ID)%XzVLTE_F,        ATM(ID)%NXzV_F,
	1       RONE,                    ATM(ID)%XzVLEVNAME_F,
	1       ATM(ID)%EDGEXzV_F,       ATM(ID)%GXzV_F,
	1       RONE,                    ATM(ID)%NXzV_F,  ND,
	1       ATM(ID)%ZXzV,            IONE,   ATM(ID)%XzV_PRES)
        END DO
!
!
! Compute LTE populations
!
	INCLUDE 'EVAL_LTE_FULL.INC'
!
! 
!
! Set up lines that will be treated with the continuum calculation.
! This section of code is also used by the code treating purely lines
! (either single transition Sobolev or CMF, or overlapping Sobolev).
!
! To define the line transitions we need to operate on the FULL atom models.
! We thus perform separate loops for each species.
!
	ML=0			!Initialize line counter.
	ALLOCATE (VEC_TRANS_NAME(N_LINE_MAX))
!
	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES)THEN
	    DO MNL_F=1,ATM(ID)%NXzV_F-1
	      DO MNUP_F=MNL_F+1,ATM(ID)%NXzV_F
	        IF(ATM(ID)%AXzV_F(MNL_F,MNUP_F) .NE. 0)THEN
	          ML=ML+1
	          IF(ML .GT. N_LINE_MAX)THEN
	            WRITE(T_OUT,*)'Error in MAINGEN'
	            WRITE(T_OUT,*)'N_LINE_MAX is too small'
	            STOP
	          END IF
	          VEC_FREQ(ML)=ATM(ID)%EDGEXzV_F(MNL_F)-ATM(ID)%EDGEXzV_F(MNUP_F)
	          VEC_SPEC(ML)=ION_ID(ID)
	          VEC_MNL_F(ML)=MNL_F
	          VEC_MNUP_F(ML)=MNUP_F
	          VEC_ION_INDX(ML)=ID
	          VEC_OSCIL(ML)=ATM(ID)%AXzV_F(MNL_F,MNUP_F)
	          VEC_EINA(ML)=ATM(ID)%AXzV_F(MNUP_F,MNL_F)
	          VEC_TRANS_NAME(ML)=TRIM(VEC_SPEC(ML))//
	1           '('//TRIM(ATM(ID)%XzVLEVNAME_F(MNUP_F))//'-'//
	1           TRIM(ATM(ID)%XzVLEVNAME_F(MNL_F))//')'
	        END IF
	      END DO
	    END DO
	  END IF
	END DO
	WRITE(6,*)'Total number of line transition is',ML
!
!
! 
!
!
! Get lines and arange in numerically decreasing frequency. This will
! allow us to consider line overlap, and to include lines with continuum
! frequencies so that the can be handled automatically.
!
	N_LINE_FREQ=ML
!
	CALL INDEXX(N_LINE_FREQ,VEC_FREQ,VEC_INDX,L_FALSE)
	CALL SORTDP(N_LINE_FREQ,VEC_FREQ,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_OSCIL,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_EINA,VEC_INDX,VEC_DP_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_MNL_F,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_ION_INDX,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_MNUP_F,VEC_INDX,VEC_INT_WRK)
	CALL SORTCHAR(N_LINE_FREQ,VEC_SPEC,VEC_INDX,VEC_CHAR_WRK)
	CALL SORTCHAR(N_LINE_FREQ,VEC_TRANS_TYPE,VEC_INDX,VEC_CHAR_WRK)
	ALLOCATE (TMP_TRANS_NAME(N_LINE_FREQ))
	CALL SORTCHAR(N_LINE_FREQ,VEC_TRANS_NAME,VEC_INDX,TMP_TRANS_NAME)
	DEALLOCATE (TMP_TRANS_NAME)
!
! 
!
! Entering main secion of the subroutine which controls which OPTIONS
! are executed. Above statements are initialization only.
!
! This message will only be printed once
!
	WRITE(T_OUT,*)
	WRITE(T_OUT,"(8X,A)")'(default is to write file '//
	1    'MAIN_OPT_STR.sve)'
	WRITE(T_OUT,"(8X,A)")'(append sve=filename to '//
	1    'write a new .sve file)'
	WRITE(T_OUT,"(8X,A)")'(box=filename to '//
	1    'write a .box file containing several .sve files)'
	WRITE(T_OUT,"(8X,A)")'(.filename to read .sve file)'
	WRITE(T_OUT,"(8X,A)")'(#filename to read .box file)'
	WRITE(T_OUT,*)
!
!
! This call resets the .sve algorithm.  Specifically it sets the next
! input answer to be a main option, and all subsequent inputs to be
! sub-options.  
!
 3	CALL SVE_FILE('RESET')
!
	DESCRIPTION=' '				!As obvious.
	MAIN_OPT_STR='  '
	DEFAULT='GR'
	CALL USR_OPTION(MAIN_OPT_STR,'OPTION',DEFAULT,DESCRIPTION)
!
!   If the main option begins with a '.', a previously
!   written .sve file is read.
!
!   If the main option begins with a '#', a previously
!   written .box file is read.
!
!   If sve= is apended to the end of this main option, a new .sve file
!   is opened with the given name and the main option and all subsequent
!   sub-options are written to this file.
!
!   If box= is input then a .box file is created, which contains the name
!   of several .sve files to process.
!
!   If only a main option is given, the option and subsequent sub-options
!   are saved in a file called 'main option.sve'.  All following main
!   options are saved in separate files.
!
! Remove variable changes from main option.
!
	I=INDEX(MAIN_OPT_STR,'(')
	IF(I .EQ. 0)THEN
	  X=UC(TRIM(MAIN_OPT_STR))
	ELSE
	  X=UC(MAIN_OPT_STR(1:I-1))	!Remove line variables.
	END IF
!
! Remove possile file names etc.
!          
	I=INDEX(X,' ')
	IF(I .EQ. 0)THEN
	  X=UC(TRIM(X))
	ELSE
	  X=UC(X(1:I-1))	!Remove file names
	END IF
!
! NB: X contains the full option
!     XOPT contains the main part of the option with the species identifier
!     XSPEC is the species identifier (eg. HE2)
!
! _ is presently used to separate the OPTION (e.g. DC) from the species
! (He2). This is the only location it needs ti specified.
! 
	I=INDEX(X,'_')
	IF(I .EQ. 0)THEN
	  XSPEC=' '
	  XOPT=X
	ELSE
	  J=INDEX(X,' ')
	  XOPT=X(1:I-1)
	  XSPEC=X(I+1:J)
	END IF
!
! Check whether the species is 
!                             (i)  In code
!                             (ii) Available
!
	IF(XSPEC .NE. ' ')THEN
	  PRESENT=.FALSE.
	  IF(XSPEC .EQ. 'ALL')THEN	!Required for SMOOTH option
	    PRESENT=.TRUE.
	  END IF
!
	  DO ISPEC=1,NSPEC
	    IF(XSPEC .EQ. UC(SPECIES(ISPEC)))THEN
	      PRESENT=.TRUE.
	      IF(POPDUM(ND,ISPEC) .EQ. 0)THEN
	        WRITE(T_OUT,*)'Error --- this  species is unavailable'
	        GOTO 1
	      END IF
	    END IF
	   END DO
!
	  DO ID=1,NUM_IONS
	    IF(XSPEC .EQ. UC(ION_ID(ID)))THEN
	      PRESENT=.TRUE.
	      IF(.NOT. ATM(ID)%XzV_PRES)THEN
	        WRITE(T_OUT,*)'Error --- this  ionization stage is unavailable'
	        GOTO 1
	      END IF
	    END IF
	  END DO
!
	  IF(.NOT. PRESENT)THEN
	    WRITE(T_OUT,*)'Error --- unknown species or ionization stage'
	    WRITE(T_OUT,*)'      --- ',XSPEC
	    WRITE(T_OUT,*)'Species may not be linked into DISPLAY package'
	    GOTO 1
	  END IF
	END IF
!
! 
!
! This IF statment is a preliminary to later work. It is
! used when line opacities are required. Only of use
! for options operating on lines. (First call is input NL and NUP).
!
! Line options.
!
	IF(XOPT .EQ. 'NETR'
	1          .OR. XOPT .EQ. 'MOMR'
	1          .OR. XOPT .EQ. 'SOBR'
	1          .OR. XOPT .EQ. 'EW'
	1          .OR. XOPT .EQ. 'LAM'
	1          .OR. XOPT .EQ. 'SRCE'
	1          .OR. XOPT .EQ. 'SRCEBB'
	1          .OR. XOPT .EQ. 'SRCEJC'
	1          .OR. XOPT .EQ. 'EP'
	1          .OR. XOPT .EQ. 'CHIL'
	1          .OR. XOPT .EQ. 'TAUL'
	1          .OR. XOPT .EQ. 'WRL'
	1          .OR. XOPT .EQ. 'BETA')THEN
!
	   CALL GET_LINE_INDICES(VEC_FREQ,VEC_MNL_F,VEC_MNUP_F,VEC_SPEC,VEC_TRANS_NAME,
	1                           N_LINE_FREQ,XSPEC,NL,NUP,FLAG)
	   LEV(1)=NL; LEV(2)=NUP
	   IF(.NOT. FLAG)GOTO 1
!
! Compute line opacity and emissivity.
!
	  FLAG=.FALSE.
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES .AND. XSPEC .EQ. UC(ION_ID(ID)))THEN
	      IF(LEV(1) .LE. ATM(ID)%NXzV_F .AND. LEV(2) .LE. ATM(ID)%NXzV_F)THEN
	        FREQ=ATM(ID)%EDGEXzV_F(LEV(1))-ATM(ID)%EDGEXzV_F(LEV(2))
	        OSCIL=ATM(ID)%AXzV_F(LEV(1),LEV(2))
	        DO I=1,ND
	          T1=ATM(ID)%W_XzV_F(LEV(2),I)/ATM(ID)%W_XzV_F(LEV(1),I)
	          GLDGU=ATM(ID)%GXzV_F(LEV(1))/ATM(ID)%GXzV_F(LEV(2))
	          CHIL(I)=OPLIN*ATM(ID)%AXzV_F(LEV(1),LEV(2))*(
	1           T1*ATM(ID)%XzV_F(LEV(1),I)-GLDGU*ATM(ID)%XzV_F(LEV(2),I) )
	          ETAL(I)=EMLIN*FREQ*ATM(ID)%AXzV_F(LEV(2),LEV(1))*ATM(ID)%XzV_F(LEV(2),I)
	        END DO
	        CALL WRITE_LINE(LEV(1),LEV(2),FREQ,ION_ID(ID))
	        AMASS=AT_MASS(SPECIES_LNK(ID))
	        FLAG=.TRUE.
	      ELSE
	        WRITE(T_OUT,*)'Levels outside valid range'
	        WRITE(T_OUT,*)'Maximum level is',ATM(ID)%NXzV_F
	        GOTO 1
	      END IF
	    END IF
	  END DO
!
	  IF(.NOT. FLAG)THEN
	    WRITE(T_OUT,*)' Invalid population type or species unavailable.'
	    GOTO 1
	  END IF
!
	  IF(OSCIL .EQ. 0)THEN
	    WRITE(T_OUT,*)'Oscillator zero for this transition'
	    WRITE(T_OUT,*)'Rerturning to main loop for new option'
	    GOTO 1
	  END IF
	  IF(.NOT. FLAG)THEN
	    WRITE(T_OUT,*)'Levels outside permitted range'
	    WRITE(T_OUT,*)'Rerturning to main loop for new option'
	    GOTO 1
	  END IF
!
! Adjust the line opacities and emissivities for the effect of clumping.
!
	  DO I=1,ND
	    CHIL(I)=CHIL(I)*CLUMP_FAC(I)
	    ETAL(I)=ETAL(I)*CLUMP_FAC(I)
	  END DO
!
	END IF
!
! 
!
	IF(XOPT .EQ. 'DIE')THEN
	  CALL USR_OPTION(LEV,1,1,'LEVS',' ','Lower level ')
	  CALL USR_OPTION(EINA,'EINA',' ','Einstein A coefficient')
	  CALL USR_OPTION(GUPDIE,'G_UP',' ','G for Autoionizing level')
	  CALL USR_OPTION(FL,'FL',' ','Transition frequency')
!
	  FLAG=.FALSE.
	  DO ID=1,NUM_IONS
	    IF(XSPEC .EQ. UC(ION_ID(ID)))THEN
	      EDGEDIE=ATM(ID)%EDGEXzV_F(LEV(1))-FL
	      CALL LTEPOP(NUST,ED,ATM(ID)%DXzV_F,GUPDIE,EDGEDIE,T,
	1                 ATM(ID)%GIONXzV_F,1,ND)
	      GLOW=ATM(ID)%GXzV_F(LEV(1))
	      DO I=1,ND
	        TA(I)=ATM(ID)%XzV_F(LEV(1),I)*ATM(ID)%W_XzV_F(LEV(1),I)
	      END DO
	      FLAG=.TRUE.
	    END IF
	  END DO
!
	  IF(.NOT. FLAG)THEN
	    WRITE(T_OUT,*)' Invalid Population type'
	    GOTO 1
	  END IF
!
	  GLDGU=GLOW/GUPDIE
	  OSCIL=EINA*EMLIN/( GLDGU*OPLIN*TWOHCSQ*(FL**2) )
	  T1=OSCIL*OPLIN
	  T2=FL*EINA*EMLIN
!
! The effect of clumping is allowed for.
!
	  DO I=1,ND
	    CHIL(I)=T1*( TA(I)-GLDGU*NUST(I) )*CLUMP_FAC(I)
	    ETAL(I)=T2*NUST(I)*CLUMP_FAC(I)
	  END DO
	END IF
!
! 
!
! **************************************************************************
! **************************************************************************
!
! Continuum options requiring an INPUT frequency, and an indication of whether
! electron scattering is to be included in the opacity. Electron scattering 
! opacity should be included for all radiation field options.
!
! **************************************************************************
! **************************************************************************
!
	IF( XOPT .EQ. 'OP'
	1          .OR. XOPT .EQ. 'KAPPA'
	1          .OR. XOPT .EQ. 'ETA'
	1          .OR. XOPT .EQ. 'XTAUC'
	1          .OR. XOPT .EQ. 'TAUC'
	1          .OR. XOPT .EQ. 'NEWRG'
	1          .OR. XOPT .EQ. 'DTAUC'
	1          .OR. XOPT .EQ. 'WRC'
	1          .OR. XOPT .EQ. 'JEXT'
	1          .OR. XOPT .EQ. 'CJBB')THEN
	  CALL USR_OPTION(FREQ,'LAM','0.0',FREQ_INPUT)
	  IF(FREQ .EQ. 0)THEN
	    CALL ESOPAC(CHI,ED,ND)		!Electron scattering
	    DO I=1,ND
	       CHI(I)=CHI(I)*CLUMP_FAC(I)
	    END DO
	    ELEC=.TRUE.
	  ELSE
	    IF(KEV_INPUT)THEN
	      FREQ=FREQ*KEV_TO_HZ
	    ELSE IF(ANG_INPUT)THEN
	      FREQ=ANG_TO_HZ/FREQ
	    END IF
	    CALL USR_OPTION(ELEC,'ELEC','T','Include electron scattering?')
	    CALL USR_OPTION(INC_RAY_SCAT,'RAY','F','Include Rayeigh scattering')
	  END IF
	END IF
!
! 
!
! **************************************************************************
! **************************************************************************
!
! Compute continuum opacity for options which require continuum
! opacity.
!
! **************************************************************************
! **************************************************************************
!
	IF(  (XOPT .EQ. 'NETR'
	1          .OR. XOPT .EQ. 'MOMR'
	1          .OR. XOPT .EQ. 'SOBR'
	1          .OR. XOPT .EQ. 'EW'
	1          .OR. XOPT .EQ. 'SRCE'
	1          .OR. XOPT .EQ. 'EP'
	1          .OR. XOPT .EQ. 'BETA'
	1          .OR. XOPT .EQ. 'ETA'
	1          .OR. XOPT .EQ. 'OP'
	1          .OR. XOPT .EQ. 'KAPPA'
	1          .OR. XOPT .EQ. 'TAUC'
	1          .OR. XOPT .EQ. 'DTAUC'
	1          .OR. XOPT .EQ. 'XTAUC'
	1          .OR. XOPT .EQ. 'NEWRG'
	1          .OR. XOPT .EQ. 'WRC'
	1          .OR. XOPT .EQ. 'WRL'
	1          .OR. XOPT .EQ. 'JEXT'
	1          .OR. XOPT .EQ. 'CJBB'
	1          .OR. XOPT .EQ. 'DIE')
	1          .AND. FL .NE. 0  )THEN
!
! Compute continuum opacity and emissivity at the line frequency.
!
	  INCLUDE 'OPACITIES.INC'
!	  INCLUDE 'HOME:[jdh.cmf.carb.test]OPACITIES.INC'
!
! Compute DBB and DDBBDT for diffusion approximation. DBB=dB/dR
! and DDBBDT= dB/dTR .
!
	  T1=HDKT*FL/T(ND)
	  T2=1.0D0-EMHNUKT(ND)
	  DBB=TWOHCSQ*( FL**3 )*T1*DTDR/T(ND)*EMHNUKT(ND)/(T2**2)
	  IF(.NOT. DIF)DBB=0.0D0
!
! Solve for the continuum radiation field at the line frequency.
!
	  THICK=.TRUE.
	  S1=ETA(1)/CHI(1)
!
! Adjust the opacities and emissivities for the influence of clumping.
!
	  DO I=1,ND
	    ETA(I)=ETA(I)*CLUMP_FAC(I)
	    CHI(I)=CHI(I)*CLUMP_FAC(I)
	    ESEC(I)=ESEC(I)*CLUMP_FAC(I)
	    CHI_RAY(I)=CHI_RAY(I)*CLUMP_FAC(I)
	    CHI_SCAT(I)=CHI_SCAT(I)*CLUMP_FAC(I)
	  END DO
!
	END IF
! 
!
! **************************************************************************
! **************************************************************************
!
! Define reinded radius grid, and evaluate quadrature weights for angular
! integrations.
!
! **************************************************************************
! **************************************************************************
!
	IF(   (XOPT .EQ. 'GREY' .OR. XOPT .EQ. 'JEXT'
	1                         .OR. XOPT .EQ. 'EWEXT'
	1                         .OR. XOPT .EQ. 'INTERP')
	1    .AND. .NOT. REXT_COMPUTED)THEN
	  REXT_COMPUTED=.TRUE.
	  CALL USR_HIDDEN(NPINS,'NPINS','2','0, 1 or 2')
	  NDX=(ND-1)*NPINS+ND
	  I=ND-10		!Parabolic interp for > I
!
! NTERP has been set to ND
!
	  CALL REXT_COEF_V2(REXT,COEF,INDX,NDX,R,GRID,ND,
	1                     NPINS,.TRUE.,I,IONE,ND)
	  NCX=NC+10
	  NPX=NDX+NCX-2
	  CALL IMPAR(PEXT,REXT,R(ND),NCX,NDX,NPX)
!
! Note that the T? vectors (here used as dummy variables) must be at least
! NPX long.
!
	  IF(TRAPFORJ)THEN
	    CALL GENANGQW(JQWEXT,REXT,PEXT,TA,TB,TC,NCX,NDX,NPX,
	1                 JTRPWGT,.FALSE.)
	    CALL GENANGQW(KQWEXT,REXT,PEXT,TA,TB,TC,NCX,NDX,NPX,
	1                 KTRPWGT,.FALSE.)
	  ELSE
	    CALL GENANGQW(JQWEXT,REXT,PEXT,TA,TB,TC,NCX,NDX,NPX,
	1                 JWEIGHT,.FALSE.)
	    CALL GENANGQW(KQWEXT,REXT,PEXT,TA,TB,TC,NCX,NDX,NPX,
	1                 KWEIGHT,.FALSE.)
	  END IF
	END IF
!
! 
!
	IF(  (XOPT .EQ. 'NETR'
	1          .OR. XOPT .EQ. 'MOMR'
	1          .OR. XOPT .EQ. 'SOBR') .AND. FIRST_RATE)THEN
	  FIRST_RATE=.FALSE.
	  CALL GEN_ASCI_OPEN(LU_NET,'RATES','UNKNOWN',' ',' ',IZERO,IOS)
	END IF
!
	IF(XOPT .EQ. 'NETR'
	1          .OR. XOPT .EQ. 'MOMR')THEN
!
	  CALL USR_OPTION(NLF,'NLF','25','Number of requency points in line profile?')
	  IF(NLF .GT. NLF_MAX)THEN
	    WRITE(T_OUT,*)'Error - NLF is too large - Maximum values is',
	1	 NLF_MAX
	    GOTO 1
	  END IF
!
! Set frequencies equally spaced in frequency.
!
	  T1=12.0D0/(NLF-1)
	  DO I=1,NLF
	    PFDOP(I)=6.0D0-(I-1)*T1
	  END DO
!
! ERF is used in computin the Sobolev incident intensity at the
! outer boudary. ERF = int from -inf to x of e(-x^2)/sqrt(pi).
! Note that ERF is not the error function. ERF is related to the
! complementary error function by ERF =-0.5D0 . erfc(X).
! S15ADF is a NAG routine which returns erfc(x).
!
! The incident Sobolev intensity is S[ 1.0-exp(-tau(sob)*ERF) ]
!
	  J=0
	  DO I=1,NLF
	    ERF(I)=-0.5D0*S15ADF(PFDOP(I),J)
	  END DO
!
	  CALL USR_OPTION(TDOP,'TDOP','0.0',
	1    'Doppler temperature for Line profile (units 10^4K)?')
	  CALL USR_OPTION(VTURB,'VTURB','10.0',
	1      'Turbulent velcity for Line profile (units km/s)?')
!
	  T1=4.286299D-05*SQRT( TDOP/AMASS + (VTURB/12.85)**2 )
	  DO I=1,NLF
	    PF(I)=PFDOP(I)*T1
	  END DO
	  CALL TRAPUNEQ(PF,LFQW,NLF)
	  T1=0.0
	  DO ML=1,NLF
	    LFQW(ML)=LFQW(ML)*FL*1.0D+15
	    LINE_PRO(ML)=DOP_PRO(PF(ML),FL,TDOP,VTURB,AMASS)
	    T1=T1+LFQW(ML)*LINE_PRO(ML)
	  END DO
!
! Normalize frequency weights.
!
	  DO ML=1,NLF
	    LFQW(ML)=LFQW(ML)/T1
	  END DO
!
! Ensure total opacity is positive.
!
	  J=(NLF+1)/2.0D0
	  DO I=1,ND
	    T2=CHI(I)+CHIL(I)*LINE_PRO(J)
	    IF( T2 .LT. 0.0D0)THEN
	      T1=3.0D-10*CHIL(I)*R(I)/V(I)/FL
	      WRITE(T_OUT,
	1      "(' Negative total opacity at d= ',I3,' :  Tausob=',
	1      1P,E12.4,' :   X/Xc=',1P,E12.4)")I,T1,T2/CHI(I)
	      CHIL(I)=-0.8D0*CHI(I)/LINE_PRO(J)
	    END IF
	  END DO
	END IF
!
! 
!
! This section compute the CONTINUUM mean intensity using EDDINGTON factors, or
! a full ray-by ray solution.
!
	IF(XOPT .EQ. 'NETR'
	1          .OR. XOPT .EQ. 'MOMR'
	1          .OR. XOPT .EQ. 'SOBR'
	1          .OR. XOPT .EQ. 'SRCEJC'
	1          .OR. XOPT .EQ. 'CJBB'
	1          .OR. XOPT .EQ. 'DIE'
	1          .OR. ( XOPT .EQ. 'EW' .AND. XOPT .NE. 'EWEXT' ))THEN
!
	  EDDC=.TRUE.
	  DEFAULT='T'
	  IF(XOPT .EQ. 'MOMR')THEN
	    EDDC=.TRUE.
	    DEFAULT='T'
	  ENDIF
	  CALL USR_OPTION(EDDC,'EDDC',DEFAULT,'Use Eddington factors to compute Jc ?')
!
	  IF(EDDC)THEN
	    S1=ZETA(1)
	    DO I=1,ND
	      SOURCE(I)=ZETA(I)
	    END DO
	    INACCURATE=.TRUE.
	    HBC=0.99D0
	    THK_CONT=.FALSE.
	    DO WHILE(INACCURATE)
	      CALL FQCOMP_V2(TA,TB,TC,XM,DTAU,R,Z,P,Q,FEXT,
	1            SOURCE,CHI,DCHIDR,JQW,KQW,DBB,HBC,
	1            INBC,IC,S1,THK_CONT,INNER_BND_METH,NC,ND,NP,METHOD)
	      T1=DBB
	      IF(INNER_BND_METH .EQ. 'ZERO_FLUX' .OR. INNER_BND_METH .EQ. 'HOLLOW')T1=0.0D0
	      CALL JFEAUNEW(TA,TB,TC,DTAU,R,RJ,Q,FEXT,
	1            ZETA,THETA,CHI,T1,IC,HBC,
	1            INBC,THK_CONT,DIF,ND,METHOD)
	      S1=ZETA(1)+THETA(1)*RJ(1)
	      INACCURATE=.FALSE.
	      T1=0.0D0
	      DO I=1,ND
	       T1=MAX(ABS(FOLD(I)-FEXT(I)),T1) 
	       FOLD(I)=FEXT(I)
	       SOURCE(I)=ZETA(I)+THETA(I)*RJ(I)
	      END DO
	      IF(T1 .GT. 1.0E-05)INACCURATE=.TRUE.
	      WRITE(T_OUT,'('' Maximum fractional change is'',1PE11.4)')T1
	    END DO
	    WRITE(6,*)GAM(1),RJ(1)
	    DO I=1,ND
	      WRITE(26,'(ES18.8,6ES14.4)')R(I),ETA(I),CHI(I),THETA(I)*CHI(I),ZETA(I),THETA(I),RJ(I)
	    END DO
!
	  ELSE
!
! Compute J. DBB and S1 have been previously evaluated. VT is used
! for dCHI_dR.
!
	    CALL NEWJSOLD(TA,TB,TC,XM,WM,FB,RJ,DTAU,R,Z,P,
	1               ZETA,THETA,CHI,VT,JQW,
	1               THK_CONT,DIF,DBB,IC,NC,ND,NP,METHOD)
	  END IF
!
! Increment emissivity due to electron scatering from the continuum.
! it is assumed that electron scatering from the line does not
! contribute signifcantly to the emissivity.
!
	  DO I=1,ND
	    ETA_WITH_ES(I)=ETA(I)+RJ(I)*CHI_SCAT(I)
	  END DO
!
	END IF
! 
!
! ****************************************************************************
! ****************************************************************************
!
! Begin the individual option section.
!
! *****************************************************************************
! *****************************************************************************
!
! ****************************************************************************
! ****************************************************************************
!
! Set up X axis options -  Log(R/RP),R/RP,LOG(T(rossland)) or Ne .
!
! ****************************************************************************
! ****************************************************************************
!
	IF(XOPT .EQ. 'XLOGR')THEN
	  DO I=1,ND
	    XV(I)=LOG10(R(I)/R(ND))
	  END DO
	  XAXIS='Log(r/R\d*\u)'
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XLINR')THEN
	  DO I=1,ND
	    XV(I)=R(I)/R(ND)
	  END DO
	  XAXIS='r/R\d*\u'
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XAU')THEN
	  FLAG=.FALSE.
	  CALL USR_HIDDEN(FLAG,'LIN','F','Linear axis (def=F)')
	  T1=1.0D+10/ASTRONOMICAL_UNIT()
	  XV(1:ND)=T1*R(1:ND)
	  IF(FLAG)THEN
	    XAXIS='r(AU)'
	  ELSE
	    XV(1:ND)=LOG10(XV(1:ND))
	    XAXIS='Log r(AU)'
	  END IF
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XRSUN')THEN
	  FLAG=.FALSE.
	  CALL USR_HIDDEN(FLAG,'LIN','F','Linear axis (def=F)')
	  T1=1.0D+10/RAD_SUN()
	  XV(1:ND)=T1*R(1:ND)
	  IF(FLAG)THEN
	    XAXIS='r(R\dsun\u)'
	  ELSE
	    XV(1:ND)=LOG10(XV(1:ND))
	    XAXIS='Log r(R\d\sun\u)'
	  END IF
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XARC')THEN
!
! NB: 1 AU corresponds to exactly 1" at 1kpc.
!
	  FLAG=.FALSE.
	  DEFAULT=WR_STRING(DIST_KPC)
	  CALL USR_OPTION(DIST_KPC,'DIST',' ','Distance to star in kpc')
	  CALL USR_HIDDEN(FLAG,'LIN','F','Linear axis (def=F)')
	  T1=1.0D-03*(1.0D+10/ASTRONOMICAL_UNIT())/DIST_KPC
	  XV(1:ND)=T1*R(1:ND)
	  IF(FLAG)THEN
	    XAXIS='r(")'
	  ELSE
	    XV(1:ND)=LOG10(XV(1:ND))
	    XAXIS='Log r(")'
	  END IF
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XED')THEN
	  DO I=1,ND
	    XV(I)=LOG10(ED(I))
	  END DO
	  XAXIS='Log(N\de\u)'
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XTEMP')THEN
	  DO I=1,ND
	    XV(I)=T(I)
	  END DO
	  XAXIS='T(10\u4\dK)'
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XN')THEN
	  DO I=1,ND
	    XV(I)=I
	  END DO
	  XAXIS='I'
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XATOM')THEN
	  DO I=1,ND
	    XV(I)=LOG10(POP_ATOM(I))
	  END DO
	  XAXIS='Log(N\di\u)'
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XVEL')THEN
	  DO I=1,ND
	    XV(I)=V(I)
	  END DO
	  XAXIS='V(km s\u-1\d)'
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XLOGV')THEN
	  CALL DLOGVEC(V,XV,ND)
	  XAXIS='Log V(kms\u-1\d)'
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XCOLD')THEN
!
	  DO I=1,ND
	    ZETA(I)=1.0D+10*MASS_DENSITY(I)*CLUMP_FAC(I)
	  END DO
	  CALL TORSCL(TA,ZETA,R,TB,TC,ND,METHOD,TYPE_ATM)
	  DO I=1,ND
	    XV(I)=DLOG10(TA(I))
	  END DO
	  XAXIS='Log m(gm cm\u-2\d)'
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XMASS')THEN
!
	  CALL USR_OPTION(FLAG,'MDIR','T','Integrate inwards from outer boundary')
	  DO I=1,ND
	    ZETA(I)=4.0D0*PI*1.0D+30*MASS_DENSITY(I)*R(I)*R(I)*CLUMP_FAC(I)
	  END DO
	  CALL TORSCL(TA,ZETA,R,TB,TC,ND,METHOD,TYPE_ATM)
	  IF(FLAG)THEN
	    WRITE(6,*)'Mass in integrated from outer boudary.'
	    WRITE(6,*)'Normalizing mass (in Msun) is',TA(ND)/1.989D+33
	    DO I=1,ND
	      XV(I)=DLOG10(TA(I)/TA(ND))
	    END DO
	    XAXIS='Log Mass Fraction'
	  ELSE
	    WRITE(6,*)'Mass in integrated from inner boudary.'
	    WRITE(6,*)'Total mass (in Msun) is',TA(ND)/1.989D+33
	    DO I=1,ND
	      XV(I)=(TA(ND)-TA(I))/1.989D+33
	    END DO
	    XAXIS='Mass (Msun)'
	  END IF
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XFLUX')THEN
	    DO I=1,ND
	      TA(I)=ABS(CLUMP_FAC(I)*FLUX_MEAN(I))
	    END DO
	    CALL TORSCL(TB,TA,R,XM,TC,ND,METHOD,TYPE_ATM)
	    WRITE(T_OUT,*)'Flux mean opacity set positive'
	    WRITE(T_OUT,*)'Flux optical depth is : ',TB(ND)
	    XV(1:ND)=LOG10(TB(1:ND))
	    XAXIS='Log(\gt\dFlux\u)'
	    XAXSAV=XAXIS
!
	ELSE  IF(XOPT .EQ. 'XTAUC')THEN
	  IF(.NOT. ELEC)THEN
	    DO I=1,ND
	      CHI(I)=CHI(I)-ESEC(I)
	    END DO
	  END IF
	  CALL TORSCL(TA,CHI,R,TB,TC,ND,METHOD,TYPE_ATM)
	  DO I=1,ND
	    XV(I)=DLOG10(TA(I))
	  END DO
	  XAXIS='Log(\gt\dc\u)'
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XCAKT')THEN
	  CALL USR_HIDDEN(T1,'VTH','10.0','Thermal doppler velocity (km/s)')
	  DO I=1,ND
	    T2=(SIGMA(I)+1.0D0)*V(I)/R(I)
	    XV(I)=DLOG10(6.65D-15*ED(I)*CLUMP_FAC(I)*T1/T2)
	  END DO
	  XAXIS='Log(t)'
!
	ELSE IF(XOPT .EQ. 'SET-ATM')THEN
	  CALL USR_OPTION(TYPE_ATM,'ATM',' ','Type of atmosphere: EXP, WIND or POW')
	  TYPE_ATM=UC(TYPE_ATM)
	  IF(TYPE_ATM .EQ. 'EXP')THEN
	  ELSE IF(TYPE_ATM .EQ. 'POW')THEN
	    TYPE_ATM='P'
	    T1=LOG(MASS_DENSITY(5)/MASS_DENSITY(1))/LOG(R(1)/R(5))
	    IF(T1 .LT. 100)THEN
	      WRITE(TYPE_ATM(2:6),'(F5.2)')T1
	      WRITE(6,'(A,F5.2)')'Density exponent for atmospheres is ',T1
	    ELSE
	      WRITE(TYPE_ATM(2:9),'(ES8.2)')T1
	      WRITE(6,'(A,ES8.2)')'Density exponent for atmospheres is ',T1
	    END IF
	  ELSE
	    TYPE_ATM=' '
	    WRITE(T_OUT,*)'Atmosphere is assumed to have a wind'
	  ENDIF
!
	ELSE IF(XOPT .EQ. 'SET-METH')THEN
	  CALL USR_OPTION(METHOD,'METH',' ','Method for derivative evaluation: LOGLOG, LOGMON, ZERO')
	  METHOD=UC(METHOD)
	  IF(METHOD .NE. 'ZERO' .AND.
	1       METHOD .NE. 'LOGLOG' .AND.
	1       METHOD .NE. 'LOGMON')THEN
	    METHOD='LOGLOG'
	    WRITE(T_OUT,*)'Did not recognize option'
	    WRITE(T_OUT,*)'Setting METHOD=LOGLOG'
	  ENDIF
!
	ELSE IF(XOPT .EQ. 'SET-IBC')THEN
	  DEFAULT=INNER_BND_METH
	  CALL USR_OPTION(INNER_BND_METH,'IBC',DEFAULT,
	1      'Inner boundary condition: DIFFUSION, SCHUSTER, ZERO_FLUX, HOLLOW')
	  IF(INNER_BND_METH .NE. 'DIFFUSION' .AND. INNER_BND_METH .NE. 'SCHUSTER' .AND.
	1       INNER_BND_METH .NE. 'ZERO_FLUX' .AND. INNER_BND_METH .NE. 'HOLLOW')THEN
	    INNER_BND_METH='DIFFSION'
	    WRITE(T_OUT,*)'Did not recognize option'
	    WRITE(T_OUT,*)'Setting INNER_BND_METH to DIFFUSION'
	  ENDIF
	  DIF=.TRUE.
	  IF(INNER_BND_METH .NE. 'DIFFUSION')DIF=.FALSE.
! 
! **************************************************************************
! **************************************************************************
!
! Information options. Provide a means of obtaining information about the
! model without using the debug option.
!
! **************************************************************************
! **************************************************************************
!

        ELSE IF(XOPT .EQ. 'WRFREQ')THEN
	  OPEN(UNIT=10,FILE='EDGE_FREQ',STATUS='NEW')
	    DO ID=1,NUM_IONS
	     IF(ATM(ID)%XzV_PRES)CALL EDGEWR(ATM(ID)%EDGEXzV_F,
	1       ATM(ID)%NXzV_F,UC(ION_ID(ID)),10)
	    END DO
	  CLOSE(UNIT=10)
	  WRITE(T_OUT,*)'Edge frequencie written to unit 10'
! 
        ELSE IF(XOPT .EQ. 'WRN')THEN
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES)WRITE(6,'(1X,A,T43,A,I4)')
	1        'Number of levels in FULL'//TRIM(ION_ID(ID)),
	1        ' model atom is:',ATM(ID)%NXzV_F
	  END DO
!
        ELSE IF(XOPT .EQ. 'WRID')THEN
!
	  CALL USR_HIDDEN(LEV,2,2,'LIMS','1,0','Imin, Imax')
!
	  FLAG=.FALSE.
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES .AND. XSPEC .EQ. UC(ION_ID(ID)))THEN
	      IF(LEV(2) .EQ. 0)LEV(2)=ATM(ID)%NXzV_F
	      WRITE(T_OUT,'(5(2X,A16,1X,I3))')(TRIM(ATM(ID)%XzVLEVNAME_F(I)),I,I=LEV(1),LEV(2))
	      FLAG=.TRUE.
	    END IF
	  END DO
	  IF(.NOT. FLAG)THEN
	    WRITE(T_OUT,*)'Error --- invalid species identification'
	  END IF
!
	ELSE IF(XOPT .EQ. 'WRLST')THEN
	  OPEN(UNIT=LU_OUT,FILE='LST_LEVL_INFO',STATUS='UNKNOWN')
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES)THEN
	      K=0
	      DO J=2,ATM(ID)%NXzV_F-1
	        DO L=1,J-1
	          IF(ATM(ID)%AXzV_F(J,L) .GE. 1.0D-20)K=K+1
	        END DO
	      END DO
	      I=ATM(ID)%NXzV_F
	      WRITE(LU_OUT,'(A7,3X,I4,I8,3X,A30)')ION_ID(ID),ATM(ID)%NXzV_F,K,ATM(ID)%XzVLEVNAME_F(I)
	    END IF
	  END DO
	  CLOSE(LU_OUT)
!
	ELSE IF(XOPT .EQ. 'WRTAB')THEN
	  OPEN(UNIT=10,FILE='NEW_TABLE',STATUS='UNKNOWN')
	  DO I=1,ND
	    TA(I)=CLUMP_FAC(I)*ROSS_MEAN(I)*(R(ND)/R(I))**2
	  END DO
	  CALL TORSCL(TAUROSS,TA,R,TB,TC,ND,METHOD,TYPE_ATM)
	  DO I=1,ND
	    WRITE(10,'(ES12.6,ES12.4,4X,F8.4,2X,ES12.4)')1.0E+10*R(I),POP_ATOM(I),T(I),TAUROSS(I)
	  END DO
!
	ELSE IF(XOPT .EQ. 'CHKA')THEN
!
	  WRITE(6,*)' '
          WRITE(6,*)BLUE_PEN//'Writing out levels with no decay routes'//DEF_PEN
	  WRITE(6,*)' '
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES)THEN
	      DO J=2,ATM(ID)%NXZV_F
	        T1=SUM(ATM(ID)%AXzV_F(J,1:J-1))
	        IF(T1 .EQ. 0.0D0)THEN
	          WRITE(6,'(A6,3X,I4,4X,A)')ION_ID(ID),J,TRIM(ATM(ID)%XzVLEVNAME_F(J))
	          WRITE(30,'(A6,3X,I4,4X,A)')ION_ID(ID),J,TRIM(ATM(ID)%XzVLEVNAME_F(J))
	        END IF
	      END DO
	    END IF
	  END DO
!
! 
!
! Change units for input of frequencies/waevengths.
!
	ELSE IF(XOPT .EQ. 'XKEV')THEN
	  KEV_INPUT=.TRUE.
	  HZ_INPUT=.FALSE.
	  ANG_INPUT=.FALSE.
	  WRITE(T_OUT,*)'Frequencies are now input in keV'
	  FREQ_INPUT='Units are keV'
	ELSE IF(XOPT .EQ. 'XHZ')THEN
	  KEV_INPUT=.FALSE.
	  HZ_INPUT=.TRUE.
	  ANG_INPUT=.FALSE.
	  WRITE(T_OUT,*)'Frequencies are now input in Hz'
	  FREQ_INPUT='Units are 10^15 Hz'
	ELSE IF(XOPT .EQ. 'XANG')THEN
	  KEV_INPUT=.FALSE.
	  HZ_INPUT=.FALSE.
	  ANG_INPUT=.TRUE.
	  WRITE(T_OUT,*)'Frequencies are now input in Angstroms.'
	  FREQ_INPUT='Units are Angstroms'
!
! 
!
	ELSE IF(XOPT .EQ. 'DIE')THEN
	  DO I=1,ND
	    ZNET(I)=1.0D0-RJ(I)*CHIL(I)/ETAL(I)
	  END DO
	  WRITE(LU_NET,40001)'Znet for DIELECTRONIC transition'
	  WRITE(LU_NET,40002)LEV(1)
	  WRITE(LU_NET,40003)(ZNET(I),I=1,ND)
	  WRITE(T_OUT,40003)(ZNET(I),I=1,ND)
!
! 
	ELSE IF(XOPT .EQ. 'NETR')THEN
!
	  CALL USR_OPTION(SOL_OPT,'FORM','FG','Method: FORM, FG, or HAM?')
!
	  IF(SOL_OPT(1:2) .EQ. 'FG')THEN
	    N_TYPE='N_ON_J'
	    CALL USR_OPTION(N_TYPE,'NTYPE','N_ON_J','N_ON_J .OR. G_ONLY?')
	    IF(UC(N_TYPE(1:1)) .EQ. 'N')N_TYPE='N_ON_J'
	    IF(UC(N_TYPE(1:1)) .EQ. 'G')N_TYPE='G_ONLY'
!
	    CALL USR_OPTION(IDEPTH,10,1,'DEPTHS','0','Depths to be plotted [0 plots JBAR] ')
	    CALL USR_OPTION(RED_EXT,'RED_EXT','14.0D0','Redward extension in Doppler widths' )
	    CALL USR_OPTION(VDOP_FG_FRAC,'FG_FRAC','0','Doppler fraction for FG[V10] ')
	    CALL USR_OPTION(VDOP_MOM_FRAC,'MOM_FRAC','0','Doppler fraction for MOM[V6] ')
	    CALL USR_HIDDEN(ELEC,'THK','T','Thick outer boundary in continuum & line')
	    CALL USR_OPTION(FG_SOL_OPT,'SOL_METH','INT/INS','FG solution method: INT/INS or DIFF/INS')
!
	    PLT_J=.FALSE.; PLT_H=.FALSE.; PLT_LF=.FALSE.; PLT_LF=.FALSE.; PLT_FM=.FALSE.
	    IF(IDEPTH(1) .EQ. 0)THEN
	    ELSE
	      CALL USR_OPTION(PLT_J,'PLT_J','F','Plot J')
	      CALL USR_OPTION(PLT_H,'PLT_H','F','Plot H')
	      CALL USR_OPTION(PLT_LF,'PLT_LF','F','Plot chi(L).H')
	      IF(PLT_J .OR. PLT_H .OR. PLT_FM)THEN
	      ELSE
	        CALL USR_OPTION(PLT_FM,'PLT_FM','T','Plot force multiplier')
	      END IF
	    END IF
!
! Note that CLUMP_FAC has already been included in ETA, CHI, ESEC, ETAL, and CHIL.
!
	    CALL COMP_JBAR(ETA,CHI,CHI_SCAT,ETAL,CHIL,
	1                T,V,SIGMA,R,P,
	1                JQW,HMIDQW,KQW,NMIDQW,
	1                JBAR,DIF,DTDR,IC,ELEC,
	1                VTURB,VDOP_FG_FRAC,VDOP_MOM_FRAC,
	1                RED_EXT,AMASS,FL,
	1                METHOD,FG_SOL_OPT,N_TYPE,
	1                LUM,PLT_J,PLT_H,PLT_LF,PLT_FM,
	1                LEV,10,NLF,NC,NP,ND)
	        DO I=1,ND
	          ZNET(I)=1.0D0-JBAR(I)*CHIL(I)/ETAL(I)
	        END DO
	      YAXIS=' '
	      IF(PLT_J)YAXIS='R\u2\dJ'
	      IF(PLT_H .AND. YAXIS .NE. ' ')YAXIS=TRIM(YAXIS)//'; R\u2\ddH'
	      IF(PLT_H .AND. YAXIS .EQ. ' ')YAXIS='R\u2\dH'
	      IF(PLT_LF .AND. YAXIS .NE. ' ')YAXIS=TRIM(YAXIS)//'; R\u2\dLF'
	      IF(PLT_LF .AND. YAXIS .EQ. ' ')YAXIS='R\u2\d^LF'
	      IF(PLT_LF .AND. YAXIS .EQ. ' ')YAXIS='MT'
	   ELSE
	    IDEPTH(1)=0
	    CALL USR_OPTION(FULL_ES,'ALLES','T',
	1        'Include line photons scatterd in resonace zone ?')
	    CALL USR_HIDDEN(SKIPEW,'SKIPEW','F','Skip accurate EW computation ?')
	    CALL USR_HIDDEN(EWACC,'EWACC','0.05','required % accuracy in EW')
!
	    INACCURATE=.TRUE.
	    CNT=0
	    EWOLD=1.0D+37
	    CALL DP_ZERO(GAM,ND)
	    DO WHILE(INACCURATE)
!
! The .FALSE. option is use to indicate that we do not require the
! approximate ANR operator to be computed. Use TA for APPROX_LAM
! which is not computed anyway.
!
	      IF(SOL_OPT(1:3) .EQ. 'HAM')THEN
	        CALL HAM_FORMSOL(ETA_WITH_ES,CHI,CHI_SCAT,CHIL,ETAL,V,SIGMA,R,P,
	1                  JBAR,ZNET,TA,TB,TC,.FALSE.,
	1                  RJ,AV,CV,GAM,
	1                  EW,CONT_INT,LINE_BL,FULL_ES,
	1                  JQW,HMIDQW,
	1                  PF,LINE_PRO,LFQW,ERF,FL,DIF,DBB,IC,
	1                  1.0D0,THK_LINE,THK_CONT,NLF,NC,NP,ND,METHOD)
	      ELSE
	        CALL FORMSOL(ETA_WITH_ES,CHI,CHI_SCAT,CHIL,ETAL,V,SIGMA,R,P,
	1                  JBAR,ZNET,TA,TB,TC,.FALSE.,
	1                  RJ,AV,CV,GAM,
	1                  EW,CONT_INT,LINE_BL,FULL_ES,
	1                  JQW,HMIDQW,
	1                  PF,LINE_PRO,LFQW,ERF,FL,DIF,DBB,IC,
	1                  1.0D0,THK_LINE,THK_CONT,NLF,NC,NP,ND,METHOD)
	      END IF
	      WRITE(T_OUT,180)EW,CONT_INT,LAMVACAIR(FL)
180	      FORMAT(1X,'EW=',1PE14.6,3X,'Ic=',E14.6,3X,'Lam=',E14.6)
!
! Use a NG acceleration every four iterations to speed up convergence.
! We use FEDD for JPREV since it is not used with FORMSOL.
!
	      CNT=CNT+1
	      L=5-MOD(CNT,4)
	      IF(L .EQ. 5)L=1
	      CALL SETFORNG(GAM,FEDD,L,ND)
	      IF(L .EQ. 1)CALL NGACCEL(GAM,FEDD,ND,.TRUE.)
!
	      ERR(CNT)=200.0D0*(EW-EWOLD)/(EW+EWOLD)
	      EWOLD=EW
	      IF( ABS(ERR(CNT)) .LT. EWACC )INACCURATE=.FALSE.
	      IF(CNT .GT. 20)THEN
	        WRITE(T_OUT,'(1X,'' Max changes for looop '',I3)')
	        WRITE(T_OUT,'(1P,(2X,5E12.3))')(ERR(I),I=1,CNT)
	        INACCURATE=.FALSE.
	      END IF
	      IF(SKIPEW)INACCURATE=.FALSE.
!
	    END DO
	  END IF
!
	  DO I=1,ND
	    YV(I)=LOG10(JBAR(I))
	  END DO
	  IF(IDEPTH(1) .EQ. 0)THEN
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='Log(Jbar)'
	  END IF
!
	  WRITE(LU_NET,40001)'FORMSOL Transfer Solution'
	  WRITE(LU_NET,40002)LEV(1),LEV(2)
	  WRITE(LU_NET,40003)(ZNET(I),I=1,ND)
! 
	ELSE IF(XOPT .EQ. 'MOMR')THEN
!
	  CALL USR_HIDDEN(HAM,'HAM','F','Use second order differencing ?')
	  CALL USR_OPTION(FULL_ES,'ALLES','T',
	1      'Include line photons scatterd in resonace zone ?')
	  CALL USR_HIDDEN(SKIPEW,'SKIPEW','F','Skip accurate EW computation ?')
!
	  INACCURATE=.TRUE.
	  CNT=0
	  EWOLD=1.0D+37
	  CALL DP_ZERO(JNU,ND*(NLF+1))
	  DO WHILE(INACCURATE)
!
! As ETAL and CHIL are known, the FORMAL solution gives FEDD and GEDD
! directly. We only need to iterate to obtian a reliable EW value in
! the presnece of electron scattering. We use ETAEDD for consistency
! purposes - ETAEDD contains the continuum electron scattering term
! computed using Eddington Factors.
!
	    IF(HAM)THEN
	      CALL FG_HAM(ETA_WITH_ES,CHI,CHI_SCAT,RJ,
	1                  CHIL,ETAL,V,SIGMA,R,P,
	1                  JNU,HNU,FEDD,GEDD,
	1                  JQW,HMIDQW,KQW,NMIDQW,
	1                  INBC_VEC,HBC_VEC,NBC_VEC,
	1                  PF,LINE_PRO,LFQW,ERF,FL,
	1                  EW,CONT_INT,LINE_BL,
	1                  DIF,DBB,IC,METHOD,
	1                  THK_CONT,THK_LINE,NLF,NC,NP,ND)
	      WRITE(T_OUT,*)'FGHAM_EW=',EW,'Ic=',CONT_INT
!
! We use TA for RADEQ, and TB for the FLUX vectors returned by the
! EW computation.
!
	      CALL MOMHAM(ETA_WITH_ES,CHI,CHI_SCAT,THETA,RJ,CHIL,ETAL,
	1                  V,SIGMA,R,JBAR,ZNET,
	1                  JNU,HNU,FEDD,GEDD,
	1                  HBC_VEC,INBC_VEC,NBC_VEC,TA,TB,
	1                  PF,LINE_PRO,LFQW,FL,DIF,DBB,IC,METHOD,
	1                  MOMEW,MOMCONT_INT,LINE_BL,FULL_ES,
	1                  NLF,NC,NP,ND)
	      WRITE(T_OUT,180)MOMEW,MOMCONT_INT,LAMVACAIR(FL)
	    ELSE
	      CALL FG_COMP(ETA_WITH_ES,CHI,CHI_SCAT,RJ,
	1                  CHIL,ETAL,V,SIGMA,R,P,
	1                  JNU,HNU,FEDD,GEDD,
	1                  JQW,HMIDQW,KQW,NMIDQW,
	1                  INBC_VEC,HBC_VEC,NBC_VEC,
	1                  PF,LINE_PRO,LFQW,ERF,FL,
	1                  EW,CONT_INT,LINE_BL,
	1                  DIF,DBB,IC,METHOD,
	1                  THK_CONT,THK_LINE,NLF,NC,NP,ND)
	      WRITE(T_OUT,*)'FGCOMP_EW=',EW,'Ic=',CONT_INT
!
! We use TA for RADEQ, and TB for the FLUX vectors returned by the
! EW computation.
!
	      CALL MOMJBAR(ETA_WITH_ES,CHI,CHI_SCAT,THETA,RJ,CHIL,ETAL,
	1                  V,SIGMA,R,JBAR,ZNET,
	1                  JNU,HNU,FEDD,GEDD,
	1                  HBC_VEC,INBC_VEC,NBC_VEC,TA,TB,
	1                  PF,LINE_PRO,LFQW,FL,DIF,DBB,IC,METHOD,
	1                  MOMEW,MOMCONT_INT,LINE_BL,FULL_ES,
	1                  NLF,NC,NP,ND)
	      WRITE(T_OUT,180)MOMEW,MOMCONT_INT,LAMVACAIR(FL)
	    END IF
!
	    CNT=CNT+1
	    ERR(CNT)=200.0D0*(MOMEW-EWOLD)/(MOMEW+EWOLD)
	    EWOLD=MOMEW
	    IF( ABS(ERR(CNT)) .LT. 0.001 )INACCURATE=.FALSE.
	    IF(CNT .GT. 20)THEN
	      WRITE(T_OUT,'(1X,'' Max changes for looop '',I3)')
	      WRITE(T_OUT,'(1P,(2X,5E12.3))')(ERR(I),I=1,CNT)
	      INACCURATE=.FALSE.
	    END IF
	    IF(SKIPEW)INACCURATE=.FALSE.
!
	  END DO
!
	  WRITE(LU_NET,40001)'MOMHAM Transfer Solution'
	  WRITE(LU_NET,40002)LEV(1),LEV(2)
	  WRITE(LU_NET,40003)(ZNET(I),I=1,ND)
!
! 
!
	ELSE IF( XOPT .EQ. 'SOBR' )THEN
!
	  CALL USR_HIDDEN(NORAD,'NORAD','F','Ignore radiation filed ?')
!
	  DO I=1,ND
	    SOURCE(I)=ETA_WITH_ES(I)/CHI(I)
	  END DO
!
! We use U and GB for BETA and BETAC respectively.
!
	  CALL SOBJBAR(SOURCE,CHI,CHIL,ETAL,
	1                V,SIGMA,R,P,JQW,
	1                JBAR,ZNET,VB,VC,U,GB,
	1                FL,DIF,DBB,IC,THK_CONT,NLF,NC,NP,ND,METHOD)
!
	  DO I=1,ND
	    YV(I)=LOG10(JBAR(I))
	  END DO
	  YAXIS='Log(Jbar)'
	  CALL DP_CURVE(ND,XV,YV)
!
	  IF(NORAD)THEN
	    WRITE(LU_NET,40001)
	1           'Sobolev Approx but radiation field not included'
	    WRITE(LU_NET,40002)LEV(1),LEV(2)
	    WRITE(LU_NET,40003)(U(I),I=1,ND)
	  ELSE
	    WRITE(LU_NET,40001)'Sobolev Approximation'
	    WRITE(LU_NET,40002)LEV(1),LEV(2)
	    WRITE(LU_NET,40003)(ZNET(I),I=1,ND)
	  END IF
!
! 
!
! This section recomputes the LTE populations. Useful for adjusting
! the Temperature and hence the population levels to assist in obtaining
! a converged model.
!
	ELSE IF(XOPT .EQ. 'HMIN')THEN
	  DO I=1,ND
	    YV(I)=2.07D-22*ED(I)*ATM(1)%XzV_F(1,I)*EXP(HDKT*0.754D0/4.135667D0/T(I))
	    YV(I)=LOG(YV(I)/2.0D0/(T(I)**1.5D0))
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
!
	ELSE IF(XOPT .EQ. 'E')THEN
	  INCLUDE 'EVAL_LTE_FULL.INC'
!
! 
!
	ELSE IF(XOPT .EQ. 'EXTRAP')THEN
!
	  CALL USR_OPTION(LEV,ITWO,ITWO,'EWIN','1,1',
	1        'Depth band to be extrapolated over')
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES)THEN
	      DO I=LEV(1),LEV(2)
	        T1=POP_ATOM(I)/POP_ATOM(LEV(2)+1)
	        DO J=1,ATM(ID)%NXzV_F
	          ATM(ID)%XzV_F(J,I)=T1*ATM(ID)%XzV_F(J,LEV(2)+1)
	        END DO
	        ATM(ID)%DXzV_F(I)=T1*ATM(ID)%DXzV_F(LEV(2)+1)
	      END DO
	    END IF
	  END DO
	  T(LEV(1):LEV(2))=T(LEV(2)+1)
!
	  INCLUDE 'EVAL_LTE_FULL.INC'
!
	ELSE IF(XOPT .EQ. 'WRDC' .OR. XOPT .EQ. 'WRPOP' .OR.
	1       XOPT .EQ. 'WRLTE' .OR. XOPT .EQ. 'WRTX')THEN
	  IF(XOPT .EQ. 'WRDC')THEN
	    TYPE='DC'
	  ELSE IF(XOPT .EQ. 'WRTX')THEN
	    TYPE='TX'
	  ELSE
	    TYPE='POP'
	  END IF
	  CALL USR_HIDDEN(STRING,'EXT',TYPE,'File appendage')
	  WRITE(T_OUT,*)STRING
!
! Set depth vector. The default allows will write all depths. We
! can omit ceartain depths using the OWIN option.
!
	  DO_DPTH(1:ND)=.TRUE.
	  CALL USR_HIDDEN(LEV,ITEN,ITWO,'OWIN','0,0,0,0,0,0,0,0,0,0',
	1        'Depth bands (pairs)NOT to be output [Max of 5]')
	  DO I=1,9,2
	    IF(LEV(I) .GT. 0 .AND. LEV(I+1) .GE. LEV(I) .AND.
	1        LEV(I+1) .LE. ND)THEN
	      DO_DPTH(LEV(I):LEV(I+1))=.FALSE.
	    END IF
	  END DO      
!
! Write out departure coefficients to ASCI file.
! NB - 1 refers to dimension of DHYD (i.e. DHYD(1,nd)
!      1 refers to format for output.
!      1,NHY - For use with HEI.
!
	  IF(XOPT .EQ. 'WRLTE')THEN
	    DO ID=1,NUM_IONS
	      IF(ATM(ID)%XzV_PRES)THEN
	        FILENAME=TRIM(ION_ID(ID))//TRIM(STRING)
	        CALL NEW_WRITEDC_V5(ATM(ID)%XzVLTE_F,ATM(ID)%LOG_XzVLTE_F,ATM(ID)%W_XzV_F,
	1             ATM(ID)%EDGEXzV_F,ATM(ID)%GXzV_F,ATM(ID)%NXzV_F,
	1             ATM(ID)%DXzV_F,ATM(ID)%GIONXzV_F,IONE,R,T,ED,V,CLUMP_FAC,
	1             DO_DPTH,LUM,ND,FILENAME,TYPE,IONE)
	      END IF
	    END DO
	  ELSE
	    DO ID=1,NUM_IONS
	      IF(ATM(ID)%XzV_PRES)THEN
	        FILENAME=TRIM(ION_ID(ID))//TRIM(STRING)
	        CALL NEW_WRITEDC_V5(ATM(ID)%XzV_F,ATM(ID)%LOG_XzVLTE_F,ATM(ID)%W_XzV_F,
	1             ATM(ID)%EDGEXzV_F,ATM(ID)%GXzV_F,ATM(ID)%NXzV_F,
	1             ATM(ID)%DXzV_F,ATM(ID)%GIONXzV_F,IONE,R,T,ED,V,CLUMP_FAC,
	1             DO_DPTH,LUM,ND,FILENAME,TYPE,IONE)
	      END IF
	    END DO
	  END IF
!
!
! 
!
! This page computes the Rosseland mean opacity from the temperature
! distribution and the population levels. An option allows the electron
! scattering opacity to be excluded from the calculations. The Rossland
! optical depth scale is given in TAUROSS, and TA is a working vector. The
! Rossland opacity is given in CHIROSS. Not written as a subroutine to allow
! easy inclusion of additional opacity sources.
!
	ELSE IF(XOPT .EQ. 'SXROSS' .OR. XOPT .EQ. 'SYROSS')THEN
	  DO I=1,ND
	    TA(I)=CLUMP_FAC(I)*ROSS_MEAN(I)*(R(ND)/R(I))**2
	  END DO
	  CALL TORSCL(TAUROSS,TA,R,TB,TC,ND,METHOD,TYPE_ATM)
	  WRITE(T_OUT,*)'Spherical rossland optical depth is : ',TAUROSS(ND)
	  IF(X .EQ. 'SXROSS')THEN
	    DO I=1,ND
	      XV(I)=LOG10(TAUROSS(I))
	    END DO
	    XAXIS='Log(\gt\dSRoss\u)'
	    XAXSAV=XAXIS
	  ELSE IF(XOPT .EQ. 'SYROSS')THEN
	    DO I=1,ND
	      YV(I)=LOG10(TAUROSS(I))
	    END DO
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='Log(\gt\dSRoss\u)'
	  END IF
!
	ELSE IF(XOPT .EQ. 'XROSS' .OR. XOPT .EQ. 'YROSS' .OR.
	1                                 XOPT .EQ. 'GREY')THEN
	  IF(ROSS_MEAN(1) .EQ. 0.0D0)THEN
	    WRITE(T_OUT,*)'Rosseland mean opacity unavailabe'
	    GOTO 1
	  END IF
	  DO I=1,ND
	    CHIROSS(I)=CLUMP_FAC(I)*ROSS_MEAN(I)
	  END DO
	  CALL TORSCL(TAUROSS,CHIROSS,R,TB,TC,ND,METHOD,TYPE_ATM)
	  ROSS=.TRUE.
	  WRITE(T_OUT,*)'Rossland optical depth is : ',TAUROSS(ND)
!
! 
	  IF(XOPT .EQ. 'XROSS')THEN
	    DO I=1,ND
	      XV(I)=LOG10(TAUROSS(I))
	    END DO
	    XAXIS='Log(\gt\dRoss\u)'
	    XAXSAV=XAXIS
	  ELSE IF(XOPT .EQ. 'YROSS')THEN
	    DO I=1,ND
	      YV(I)=LOG10(TAUROSS(I))
	    END DO
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='Log(\gt\dRoss\u)'
	  ELSE
!
! We have two different methods avaialable to compute the Grey Temperature
! distribution. One takes into account the velocity terms, and is needed for
! supernovae calculations.
!
	    IF(V(1) .GT. 5.0D+03)THEN
	      CALL USR_OPTION(GREY_REL,'GREY_REL','T','Use for relativistic solution')
	    ELSE
	      CALL USR_OPTION(GREY_REL,'GREY_REL','F','Use for relativistic solution')
	    END IF
	    CALL USR_HIDDEN(GREY_WITH_V_TERMS,'VT','F','File appendage')
	    CALL USR_HIDDEN(ELEC,'LOGT','F','Log of T?')
	    IF(R(1)/R(ND) .LT. 2.0 .AND. V(1) .LT. 1.0D0)THEN
	      CALL USR_OPTION(PLANE_PARALLEL_NOV,'PPNOV','T','Plane parallel model no V?')
	    ELSE
	      CALL USR_HIDDEN(PLANE_PARALLEL_NOV,'PPNOV','F','Plane parallel model no V?')
	    END IF
!
! Compute Grey temperature distribution. 
!
! Will use FA for F, GAM for NEWRJ, and GAMH for NEWRK
!
	    YAXIS='T(10\u4 \dK)'
	    IF(GREY_REL)THEN
!
	      TA(1:ND)=0.3D0              !FEDD: Initial guess
	      TB(1:ND)=0.0D0              !H_ON_J: Initial guess
	      TC(1:ND)=0.0D0              !dlnJdlnR=0.0D0
	      GAMH(1:ND)=0.0D0            !Old FEDD
	      XM(1:ND)=0.0D0              !As grey solution, not needed (ETA)
	      NEW_FREQ=.TRUE.
	      HBC=0.7D0; INBC=0.1D0
	      WRITE(T_OUT,*)'Using MOM_JREL_GREY_V1 for grey solution'
!
! Note
!   HFLUX=LUM*Lsun/16/(PI*PI)/10**2/R**2 (10**2 for 1/R**2).
!   DBB = dBdR = 3.Chi.L/16(piR)**2 
!   DBB is used for the lower boundary diffusion approximation. 
!
	      HFLUX=3.826D+13*LUM/(4.0D0*PI*R(ND))**2
	      DBB=3.0D0*HFLUX*CHIROSS(ND)
	      IF(.NOT. DIF)DBB=0.0D0
	      T1=1.0D0
	      DO WHILE(T1 .GT. 1.0D-05)
	        CALL MOM_JREL_GREY_V1(XM,CHIROSS,CHIROSS,V,SIGMA,R,
	1              TB,TA,TC,RJ,HNU,HBC,INBC,
	1              DIF,DBB,IC,METHOD,
	1              L_TRUE,L_TRUE,NEW_FREQ,ND)
!
	        CALL FGREY_NOREL_V1(TA,TB,RJ,CHIROSS,R,V,SIGMA,
	1              P,JQW,HMIDQW,KQW,LUM,IC,METHOD,
	1              HBC,INBC,DIF,ND,NC,NP)
	        T1=0.0D0
	        DO I=1,ND
	          T1=MAX(ABS(TA(I)-GAMH(I)),T1)
	          GAMH(I)=TA(I)
	        END DO
	        NEW_FREQ=.FALSE.
	        WRITE(T_OUT,*)'Current grey iteration accuracy is',T1
	      END DO
	      CALL TORSCL(TA,CHIROSS,R,TB,TC,ND,METHOD,TYPE_ATM)
	      DO I=1,ND
	        ZV(I)=DLOG10(TA(I))
	        TGREY(I)=((3.14159265D0/5.67D-05*RJ(I))**0.25D0)*1.0D-04
	        YV(I)=TGREY(I)
	      END DO
	      GREY_COMP=.TRUE.
	      IF(ELEC)THEN
	        YV(1:ND)=LOG10(YV(1:ND))
	        YAXIS='Log T(10\u4 \dK)'
	      END IF
	      CALL DP_CURVE(ND,ZV,YV)
	      XAXIS='Log(\gt\dRoss\u)'
	      XAXSAV=XAXIS
!

	  ELSE IF(GREY_WITH_V_TERMS)THEN
	      T2=1.0D-05          !Accuracy to converge f
              CALL JGREY_WITH_FVT(RJ,TB,CHIROSS,R,V,SIGMA,
	1                  P,JQW,HMIDQW,KQW,NMIDQW,
	1                  LUM,METHOD,DIF,IC,
	1                  T2,ND,NC,NP)
	      CALL TORSCL(TA,CHIROSS,R,TB,TC,ND,METHOD,TYPE_ATM)
	      DO I=1,ND
	        ZV(I)=DLOG10(TA(I))
	        TGREY(I)=((3.14159265D0/5.67D-05*RJ(I))**0.25D0)*1.0D-04
	        YV(I)=TGREY(I)
	      END DO
	      GREY_COMP=.TRUE.
!
	      IF(ELEC)THEN
	        YV(1:ND)=LOG10(YV(1:ND))
	        YAXIS='Log T(10\u4 \dK)'
	      END IF
	      CALL DP_CURVE(ND,ZV,YV)
	      XAXIS='Log(\gt\dRoss\u)'
	      XAXSAV=XAXIS
	    ELSE IF(PLANE_PARALLEL_NOV)THEN
	      FOLD=0.3333D0
	      HBC=0.7D0; NBC=0.0D0; FL=1.0D0; INBC=0.1D0
!
! Note HFLUX=LUM*Lsun/16/(PI*PI)/10**2 (10**2 for 1/R**2).
! DBB - dBdR = 3. Chi. L/16(piR)**2 and is used for the lower boundary
! diffusion approximation. Since we are dealing with a plane-parallel
! atmopshere, we divide HFLUX by R*^2.
!
              HFLUX=3.826D+13*LUM/16.0D0/PI**2/R(ND)/R(ND)
              DBB=3.0D0*CHIROSS(ND)*HFLUX
              IF(.NOT. DIF)DBB=0.0D0
	      T1=1000.0D0
!
! Compute radial (vertical) optical depth increments.
!
	      CALL DERIVCHI(TB,CHIROSS,R,ND,METHOD)
	      CALL NORDTAU(DTAU,CHIROSS,R,R,TB,ND)
!
! Compute the solution vector. Note that the units need to be
! eventually included. The following follows direcly from d2K/d2Tau=0.
!
	      DO WHILE(T1 .GT. 1.0D-08)
	        T2=FOLD(1)/HBC
	        RJ(1)=HFLUX/HBC
	        DO I=2,ND
	          T2=T2+DTAU(I-1)
	          RJ(I)=T2*HFLUX/FOLD(I)
	        END DO
!
! Will use TA for IPLUS.
!
	        SOURCE(1:ND)=RJ
	        CALL FCOMP_PP_V2(R,TC,GAMH,SOURCE,CHIROSS,TA,HBC,
	1               NBC,INBC,DBB,IC,THK_CONT,DIF,ND,NC,METHOD)
	        T1=0.0D0
                DO I=1,ND
                  T1=MAX(ABS(FOLD(I)-GAMH(I)),T1)
                  FOLD(I)=GAMH(I)
                END DO
              END DO
!
! Compute the temperature distribution, and the Rossland optical depth scale.
! Assumes LTE. NB sigma=5.67E-05 and the factor of 1.0E-04 is to convert
! T from units of K to units of 10^4 K. 
!
	      GREY_COMP=.TRUE.
	      DO I=1,ND
	        TGREY(I)=((3.14159265D0/5.67D-05*RJ(I))**0.25D0)*1.0D-04
	      END DO
	      WRITE(T_OUT,*)'Use TGREY to plot T' 
!
	    ELSE
!
! We use a finer grid here. The finer grid has been  previously defined.
!
! Interpolate CHIROSS onto a finer grid.
!
	      DO I=1,NDX
	        CHI(I)=0.0D0
	        DO J=0,3
	          CHI(I)=CHI(I)+COEF(J,I)*DLOG( CHIROSS(J+INDX(I)) )
	        END DO
	        CHI(I)=DEXP(CHI(I))
	      END DO
!
	      DO I=1,NDX
	        FA(I)=1.0D0/3.0D0
	      END DO
	      HBC=1.0D0
	      T1=1000.0
	      DO WHILE(T1 .GT. 1.0D-05)
	        CALL JGREY(TA,TB,TC,XM,DTAU,REXT,Z,PEXT,RJ,
	1          GAM,GAMH,Q,FA,CHI,dCHIdR,
	1          JQWEXT,KQWEXT,LUM,HBC,HBCNEW,NCX,NDX,NPX,METHOD)
	        T1=0.0D0
	        DO I=1,NDX
	          T1=MAX(ABS(FA(I)-GAMH(I)),T1)
	          FA(I)=GAMH(I)
	        END DO
	        T1=MAX(ABS(HBC-HBCNEW),T1)
	        HBC=HBCNEW
	        WRITE(T_OUT,'('' Maximum change is'',1PE11.4)')T1
	      END DO
!
! Compute the temperature distribution, and the Rossland optical depth scale.
! Assumes LTE. NB sigma=5.67E-05 and the factor of 1.0E-04 is to convert
! T from units of K to units of 10^4 K. We use ZV for XV so as not to corrupt
! XV.
!
	      CALL TORSCL(TA,CHI,REXT,TB,TC,NDX,METHOD,TYPE_ATM)
	      DO I=1,NDX
	        ZV(I)=DLOG10(TA(I))
	        YV(I)=((3.14159265D0/5.67D-05*RJ(I))**0.25D0)*1.0D-04
	      END DO
!
! Store grey temperature on Model Grid.
!
	      GREY_COMP=.TRUE.
	      DO I=1,ND
	        TGREY(I)=YV(GRID(I))
	      END DO
!
	      CALL USR_HIDDEN(ELEC,'LOGT','F','Log of T?')
	      IF(ELEC)THEN
	        DO I=1,NDX
	          YV(I)=LOG10(YV(I))
	        END DO
	        YAXIS='Log T(10\u4 \dK)'
	      END IF
	      CALL DP_CURVE(NDX,ZV,YV)
	      XAXIS='Log(\gt\dRoss\u)'
	      XAXSAV=XAXIS
	    END IF
!
! Set REXT(1)=0 to insure that REXT is recomputed for the interpolation
! section.
!
!	    REXT(1)=0.0
	  END IF
!
! 
!
! ****************************************************************************
! ****************************************************************************
!
! Y axis options.
!
! ****************************************************************************
! ****************************************************************************
!
	ELSE IF(XOPT .EQ. 'YFLUX')THEN
	  CALL USR_OPTION(ELEC,'dFLUX','F','Plot dTAU (flux mean opacity)')
	  DO I=1,ND
	    TA(I)=ABS(CLUMP_FAC(I)*FLUX_MEAN(I))
	  END DO
	  CALL TORSCL(TB,TA,R,XM,TC,ND,METHOD,TYPE_ATM)
	  WRITE(T_OUT,*)'Flux optical depth is : ',TB(ND)
!
	  IF(ELEC)THEN
	    DO I=1,ND-1
	      IF(TB(I) .GT. 0)THEN
	        YV(I)=LOG10(TB(I+1)/TB(I))
	      ELSE
	        YV(I)=10.0D0
	      END IF
	    END DO
	    CALL DP_CURVE(ND-1,XV,YV)
	    YAXIS='dLog \gt(Flux)'
	  ELSE
	    DO I=1,ND
	      IF(TB(I) .GT. 0)THEN
	        YV(I)=LOG10(TB(I))
	      ELSE
	        YV(I)=0.0D0
	      END IF
	    END DO
	    YAXIS='Log \gt(Flux)'
	    CALL DP_CURVE(ND,XV,YV)
	  END IF
!
	ELSE IF(XOPT .EQ. 'ROSS')THEN
	  WRITE(T_OUT,*)'Volume filling factor not allowed for.'
	  CALL USR_OPTION(ELEC,'KAPPA','T','Mass absorption coefficient?')
	  IF(ROSS_MEAN(1) .NE. 0.0D0)THEN
	    IF(ELEC)THEN
	      DO I=1,ND
	        YV(I)=1.0D-10*ROSS_MEAN(I)/MASS_DENSITY(I)
	      END DO
	      YAXIS='Rosseland Mean Opacity (cm\u2\d/g)'           ! (cm\u-1\d)'
	    ELSE
	      CALL USR_OPTION(ELEC,'ON_NE','T','Normalize by the electron scattering opacity?')
	      IF(ELEC)THEN
	        DO I=1,ND
	          YV(I)=ROSS_MEAN(I)/(6.65D-15*ED(I))
	        END DO
	      YAXIS='Rosseland Mean Opacity/ \gsNe'           ! (cm\u-1\d)'
	      ELSE
	        DO I=1,ND
	          YV(I)=DLOG10(ROSS_MEAN(I))-10
	        END DO
	        YAXIS='Rosseland Mean Opacity (cm\u-1\d)'
	      END IF
	    END IF
	    CALL DP_CURVE(ND,XV,YV)
	  ELSE
	    WRITE(T_OUT,*)'Rosseland opacity not available.'
	  END IF
!
	ELSE IF(XOPT .EQ. 'YLOGR')THEN
	  DO I=1,ND
	    YV(I)=DLOG10(R(I)/R(ND))
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='Log(r/R\d*\u)'
!
	ELSE IF(XOPT .EQ. 'ED')THEN
	  DO I=1,ND
	    YV(I)=DLOG10(ED(I))
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='Log(N\de\u)'
!
	ELSE IF(XOPT .EQ. 'YDEN')THEN
	  DO I=1,ND
	    YV(I)=DLOG10(MASS_DENSITY(I))
	  END DO
	  YAXIS='\gr(gm cm\u-3\d)'
	  CALL DP_CURVE(ND,XV,YV)
!
	ELSE IF(XOPT .EQ. 'YCOLD')THEN
	  DO I=1,ND
	    ZETA(I)=1.0D+10*MASS_DENSITY(I)*CLUMP_FAC(I)
	  END DO
	  CALL TORSCL(TA,ZETA,R,TB,TC,ND,METHOD,TYPE_ATM)
	  DO I=1,ND
	    YV(I)=DLOG10(TA(I))
	  END DO
	  YAXIS='m(gm cm\u-2\d)'
	  CALL DP_CURVE(ND,XV,YV)
!
! Used to verify that model has a constant mass-loss rate. It can also be used to
! see Mdot variation in time dependent model.
!
	ELSE IF(XOPT .EQ. 'YMDOT')THEN
	  T1=4.0D+25*PI*365.25D0*24.0D0*3600.0D0/MASS_SUN()
!	  WRITE(6,*)T1
!	  WRITE(6,*)MASS_DENSITY(1),CLUMP_FAC(1),R(1),V(1)
	  DO I=1,ND
	    WRITE(26,*)T1*MASS_DENSITY(I)*CLUMP_FAC(I)*R(I)*R(I)*V(I)
	    YV(I)=LOG10( T1*MASS_DENSITY(I)*CLUMP_FAC(I)*R(I)*R(I)*V(I) )
	  END DO
	  WRITE(6,'(A,ES14.4,A)')'Mass loss rate is ',10**(YV(1)),' Msun/yr'
	  WRITE(6,'(A,ES14.4,A)')'     Log(Mdot) is ',YV(1)
	  YAXIS='Mass Loss rate(Msun/yr)'
	  CALL DP_CURVE(ND,XV,YV)
!
! Designed to look at Mdot in shell models. Since these models are not smooth, we avearge
! (sum) data over the smoothing step size.
!
	ELSE IF(XOPT .EQ. 'YAVMDOT')THEN
	  CALL USR_OPTION(T4,'dLOGR','0.1','Smoothing step in log R space')
	  T1=4.0D+30*PI/MASS_SUN()
	  DO I=1,ND
	    ZETA(I)=T1*MASS_DENSITY(I)*CLUMP_FAC(I)*R(I)*R(I)
	  END DO
	  CALL TORSCL(TA,ZETA,R,TB,TC,ND,METHOD,TYPE_ATM)
	  WRITE(6,*)'Done mass determination'
	  DO I=1,ND
	    ZETA(I)=1.0D+05/V(I)
	  END DO
	  CALL TORSCL(ZV,ZETA,R,TB,TC,ND,METHOD,TYPE_ATM)
	  WRITE(6,*)'Done time determination'
	  DO I=1,ND-1
	    T1=365.25D0*24.0D0*3600.0D0
	    J=MIN(I+1,ND)
	    DO K=I+1,ND
	      IF( LOG10(R(I)/R(J)) .GT. T4)EXIT
	      J=K
	    END DO
	    YV(I)=LOG10( T1*(TA(J)-TA(I))/(ZV(J)-ZV(I)) )
	  END DO
	  T1=4.0D+25*PI*365.25D0*24.0D0*3600.0D0/MASS_SUN()
	  I=ND
	  YV(I)=LOG10( T1*MASS_DENSITY(I)*CLUMP_FAC(I)*R(I)*R(I)*V(I) )
!
	  YAXIS='Mass Loss rate(Msun/yr)'
	  CALL DP_CURVE(ND,XV,YV)
!
	ELSE IF(XOPT .EQ. 'YF')THEN
	  YAXIS='Clumping factor'
	  CALL DP_CURVE(ND,XV,CLUMP_FAC)
!
! Designed to look at f, the clumping factor, in shell models. Since these models are not 
! smooth, we avearge (sum) data over the smoothing step size.
!
	ELSE IF(XOPT .EQ. 'YAVF')THEN
	  CALL USR_OPTION(T4,'dLOGR','0.1','Smoothing step in log R space')
!
	  T1=4.0D+30*PI
	  DO I=1,ND
	    ZETA(I)=T1*MASS_DENSITY(I)*R(I)*R(I)*CLUMP_FAC(I)
	  END DO
	  CALL TORSCL(TA,ZETA,R,TB,TC,ND,METHOD,TYPE_ATM)
!
	  DO I=1,ND
	    ZETA(I)=ZETA(I)*MASS_DENSITY(I)
	  END DO
	  CALL TORSCL(ZV,ZETA,R,TB,TC,ND,METHOD,TYPE_ATM)
!
! Compute dV
!
	  T1=4.0D+30*PI
	  DO I=1,ND
	    ZETA(I)=T1*R(I)*R(I)
	  END DO
	  CALL TORSCL(XM,ZETA,R,TB,TC,ND,METHOD,TYPE_ATM)
!
	  CALL DETERM_CLUM_POS(MASS_DENSITY,R,ND,TC,K)
	  J=1
	  IF(K .GT. 2)THEN
	    L=(TC(2)-TC(1))/2
	  ELSE
	    L=1
	    WRITE(6,*)'L=',L
	  END IF
	  DO I=1,ND-5
	    IMIN=MAX(1,I-L)
	    IMAX=MIN(ND,I+L)
	    YV(I)=(TA(IMAX)-TA(IMIN))*(TA(IMAX)-TA(IMIN))/(ZV(IMAX)-ZV(IMIN))/(XM(IMAX)-XM(IMIN))
	    IF(I .GT. TC(K)+L/2)THEN
	      DO J=I,ND-1
	        IMAX=J+1; IMIN=J
	        YV(J)=(TA(IMAX)-TA(IMIN))*(TA(IMAX)-TA(IMIN))/(ZV(IMAX)-ZV(IMIN))/(XM(IMAX)-XM(IMIN))
	      END DO
	      YV(ND)=1.0D0
	      EXIT
	    END IF
	  END DO
!
!	  DO I=1,ND-1
!	    J=MIN(I+1,ND)
!	    DO K=I+1,ND
!	      IF(LOG10(R(I)/R(J)) .GT. T4)EXIT
!	      J=K
!	    END DO
!	    YV(I)=(TA(J)-TA(I))*(TA(J)-TA(I))/(ZV(J)-ZV(I))/(XM(J)-XM(I))
!	  END DO
!	  YV(ND)=1.0D0
!
	  YAXIS='Clumping factor'
	  CALL DP_CURVE(ND,XV,YV)
!
	ELSE IF(XOPT .EQ. 'YMASS')THEN
	  T1=4.0D+30*PI/MASS_SUN()
	  DO I=1,ND
	    ZETA(I)=T1*MASS_DENSITY(I)*CLUMP_FAC(I)*R(I)*R(I)
	  END DO
	  CALL TORSCL(TA,ZETA,R,TB,TC,ND,METHOD,TYPE_ATM)
	  YAXIS='M(\dsun\u)'
	  WRITE(6,'(A,ES10.3,A)')'Mass of envelope (ejecta) is',TA(ND),' Msun'
	  CALL DP_CURVE(ND,XV,TA)
!
	ELSE IF(XOPT .EQ. 'IMASS')THEN
	  K=0
	  WRITE(6,'(A)')' '
	  DO ISPEC=1,NSPEC
	    IF( (XSPEC .EQ. SPECIES(ISPEC) .OR. XSPEC .EQ. 'ALL') .AND.  POPDUM(ND,ISPEC) .GT. 0.0D0)THEN
	      T1=4.0D+30*PI*AT_MASS(ISPEC)*ATOMIC_MASS_UNIT()/MASS_SUN()
	      DO I=1,ND
	        ZETA(I)=T1*POPDUM(I,ISPEC)*CLUMP_FAC(I)*R(I)*R(I)
	      END DO
	      CALL TORSCL(TA,ZETA,R,TB,TC,ND,METHOD,TYPE_ATM)
	      CALL DP_CURVE(ND,XV,TA)
	      K=K+1
	      WRITE(6,'( A,A4,A,ES9.2,A)',ADVANCE='NO')' Mass of ',TRIM(SPECIES(ISPEC)),' is',TA(ND),' Msun'
	      IF(MOD(K,2) .NE. 0)WRITE(6,'(10X)',ADVANCE='NO')
	      IF(MOD(K,2) .EQ. 0)WRITE(6,'(A)')' '
	    END IF
	    IF(MOD(K,2) .NE. 0. .AND. XSPEC .EQ. 'ALL')WRITE(6,'(A)')' '
	  END DO
	  YAXIS='Mass(M\d\(9)\u)'
!
	ELSE IF(XOPT .EQ. 'LOGT')THEN
	  CALL DLOGVEC(T,YV,ND)
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='Log(T[10\u4\dK])'
!
	ELSE IF(XOPT .EQ. 'TGREY')THEN
	  IF(GREY_COMP)THEN
	    CALL DP_CURVE(ND,XV,TGREY)
	    YAXIS='T(10\u4\dK)'
	  ELSE
	    WRITE(T_OUT,*)'Error -Grey Temperature distribution not available'
	    WRITE(T_OUT,*)'Call GREY option first --- ',
	1                'has higher spatial resolution.'
	  END IF
!
	ELSE IF(XOPT .EQ. 'ROP')THEN
	  DO I=1,ND
	    YV(I)=LOG10(1.0D-10*ROSS_MEAN(I)/MASS_DENSITY(I))
	    ZV(I)=LOG10(1.0D+06*MASS_DENSITY(I)/T(I)**3)
	    WRITE(6,*)I,ZV(I),YV(I)
	  END DO
	  CALL DP_CURVE(ND,ZV,YV)
	  XAXIS='\gr/T\u3\d\d6\u'
	  YAXIS='\gk(ross)'
!
	ELSE IF(XOPT .EQ. 'THOP')THEN
	  DO I=1,ND
	    CHIROSS(I)=CLUMP_FAC(I)*ROSS_MEAN(I)*R(ND)*R(ND)/R(I)/R(I)
	  END DO
	  CALL TORSCL(TA,CHIROSS,R,TB,TC,ND,METHOD,TYPE_ATM)
	  DO I=1,ND
	    CHIROSS(I)=CLUMP_FAC(I)*ROSS_MEAN(I)
	  END DO
	  CALL TORSCL(TAUROSS,CHIROSS,R,TB,TC,ND,METHOD,TYPE_ATM)
	  I=1
	  DO WHILE(TAUROSS(I) .LT. 0.1)
	    I=I+1
	  END DO
	  J=1
	  DO WHILE(TAUROSS(J) .LT. 1.0)
	    J=J+1
	  END DO
	  K=1
	  DO WHILE(TAUROSS(K) .LT. 3.0)
	    K=K+1
	  END DO
!	  Q0=1.0D0; QINF=0.5D0; GAM_HOPF=1.0
          T2=1.0D+10*R(ND)/RAD_SUN()
          T2=1.0D-04*TEFF_SUN()*(LUM/T2**2)**0.25                !Units of 10^4 K
	  DO I=1,ND
	    T1=1.3333333D0*( (T(I)/T2)**4 )*TAUROSS(I)/TA(I)-TAUROSS(I)
	    WRITE(6,*)I,TAUROSS(I),TA(I),T(I),T2
	    YV(I)=T1
	  END DO
	  CALL DP_CURVE(ND,TA,YV)
!
	ELSE IF(XOPT .EQ. 'TRAT')THEN
!
! Plots the ration of the ACTUAL T distribution to the GREY value.
!
	  IF(GREY_COMP)THEN
	    DO I=1,ND
	      YV(I)=T(I)/TGREY(I)
	    END DO
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='T(10\u4\dK)'
	  ELSE
	    WRITE(T_OUT,*)'Error -Grey Temperature distribution not available'
	    WRITE(T_OUT,*)'Call GREY option first --- '
	  END IF
!
	ELSE IF(XOPT .EQ. 'TCMF')THEN
	  DO I=1,ND
	    YV(I)=CMFGEN_TGREY(I)
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='T(10\u4\dK)'
!
	ELSE IF(XOPT .EQ. 'T')THEN
	  DO I=1,ND
	    YV(I)=T(I)
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='T(10\u4\dK)'
!
! Compute dlnT/dlnP for comparison with adiabatic temperature gradient.
!
	ELSE IF(XOPT .EQ. 'DTDP')THEN
	  T1=1.0D+04*BOLTZMANN_CONSTANT()
	  TA(1:ND)=T1*(ED(1:ND)+POP_ATOM(1:ND))*T(1:ND)
	  ELEC=.TRUE.
	  CALL USR_OPTION(ELEC,'RP','T','Include radiation pressure')
	  IF(ELEC)THEN
	    T1=1.0D+16*STEFAN_BOLTZ()/SPEED_OF_LIGHT()
	    TA(1:ND)=TA(1:ND)+T1*(T(1:ND)**4)
	  END IF
	  TA(1:ND)=DLOG( TA(1:ND) )
	  TB(1:ND)=DLOG( T(1:ND) )
	  CALL DERIVCHI(TC,TB,TA,ND,'LINMON')
	  YV(1:ND)=TC(1:ND)
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='dlnT/dlnP'
!
	ELSE IF(XOPT .EQ. 'DERAD')THEN
	  CALL USR_OPTION(ELEC,'INTEG','F','Integrated luminosity')
	  IF(ELEC)THEN
	    YV(1:ND)=3.280D-03*dE_RAD_DECAY(1:ND)*R(1:ND)*R(1:ND)  !(4*PI*Dex(+30)/L(sun)
	    CALL LUM_FROM_ETA(YV,R,ND)
	    DO I=ND-1,1,-1
	      YV(I)=YV(I+1)+YV(I)
	    END DO
	    YAXIS='E(rad)(L\dsun\u)'
	  ELSE
	    YV(1:ND)=DLOG10(dE_RAD_DECAY(1:ND))
	    YAXIS='Log \ge(ergs\u \dcm\u-3 \ds\u-1\d)'
	  END IF
	  CALL DP_CURVE(ND,XV,YV)
!
	ELSE IF(XOPT .EQ. 'EK')THEN
	  CALL USR_OPTION(ELEC,'INTEG','T','Integrate thermal kinetic energy?')
	  YV(1:ND)=1.5D+04*(ED(1:ND)+POP_ATOM(1:ND))*T(1:ND)*BOLTZMANN_CONSTANT()
	  IF(ELEC)THEN
	    YV(1:ND)=3.280D-03*YV(1:ND)*R(1:ND)*R(1:ND)  !(4*PI*Dex(+30)/L(sun)
	    CALL LUM_FROM_ETA(YV,R,ND)
	    DO I=ND-1,1,-1
	      YV(I)=YV(I+1)+YV(I)
	    END DO
	    YAXIS='E(kinetic)(s.L\dsun\u)'
	  ELSE
	    YV(1:ND)=DLOG10(YV(1:ND))
	    YAXIS='Log EK(ergs\u \dcm\u-3)'
	  END IF
	  CALL DP_CURVE(ND,XV,YV)
!
	ELSE IF(XOPT .EQ. 'ERAD')THEN
	  CALL USR_OPTION(ELEC,'INTEG','T','Integrate radiative energy?')
	  T1=4.0D+16*STEFAN_BOLTZ()/SPEED_OF_LIGHT()
	  YV(1:ND)=T1*(T(1:ND)**4)
	  IF(ELEC)THEN
	    YV(1:ND)=4.0D+30*PI*YV(1:ND)*R(1:ND)*R(1:ND)  !(*PI*Dex(+30) !/L(sun)
	    CALL LUM_FROM_ETA(YV,R,ND)
	    DO I=ND-1,1,-1
	      YV(I)=YV(I+1)+YV(I)
	    END DO
	    YAXIS='Log E(rad)'
	    WRITE(6,'(A,ES14.4,A)')'Total radiation density ',YV(1),' s.Lsun'
	    YV(1:ND-1)=DLOG10(YV(1:ND-1))
	    CALL DP_CURVE(ND-1,XV,YV)
	  ELSE
	    YAXIS='Log E\drad\u(ergs\u \dcm\u-3)'
	    YV(1:ND)=DLOG10(YV(1:ND))
	    CALL DP_CURVE(ND,XV,YV)
	  END IF
!
	ELSE IF(XOPT .EQ. 'EI')THEN
	  CALL USR_OPTION(ELEC,'INTEG','T','Integrate intenal energy?')
          YV(1:ND)=0.0D0
!
	  DO ISPEC=1,NSPEC
	    T1=0.0D0			!Total energy of state
	    T2=0.0D0                    !Running ioization energy
!
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      DO I=1,ATM(ID)%NXzV_F
	        T1=T2+(ATM(ID)%EDGEXzV_F(1)-ATM(ID)%EDGEXzV_F(I))
	        DO J=1,ND
	          YV(J)=YV(J)+T1*ATM(ID)%XzV_F(I,J)
	        END DO
	      END DO
	      T2=T2+ATM(ID)%EDGEXzV_F(1)
	    END DO
	  END DO
	  YV(1:ND)=1.0D+04*HDKT*YV(1:ND)*BOLTZMANN_CONSTANT()
	  IF(ELEC)THEN
	    YV(1:ND)=3.280D-03*YV(1:ND)*R(1:ND)*R(1:ND)  !(4*PI*Dex(+30)/L(sun)
	    CALL LUM_FROM_ETA(YV,R,ND)
	    DO I=ND-1,1,-1
	      YV(I)=YV(I+1)+YV(I)
	    END DO
	    YAXIS='E(kinetic)(s.L\dsun\u)'
	  ELSE
	    YV(1:ND)=DLOG10(YV(1:ND))
	    YAXIS='Log EI(ergs\u \dcm\u-3)'
	  END IF
	  CALL DP_CURVE(ND,XV,YV)
!
	ELSE IF(XOPT .EQ. 'PGONP')THEN
	  TA(1:ND)=1.0D+04*BOLTZMANN_CONSTANT()*(ED(1:ND)+POP_ATOM(1:ND))*T(1:ND)
	  T1=4.0D+16*STEFAN_BOLTZ()/SPEED_OF_LIGHT()/3.0D0
	  TB(1:ND)=T1*T(1:ND)**4
	  YV(1:ND)=TA(1:ND)/TB(1:ND)
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='Pgas/P'
!
	ELSE IF(XOPT .EQ. 'DPDR')THEN
	  TA(1:ND)=1.0D+04*BOLTZMANN_CONSTANT()*(ED(1:ND)+POP_ATOM(1:ND))*T(1:ND)
	  DO I=2,ND-1
	    YV(I)=1.0D-10*(TA(I+1)-TA(I-1))/(R(I-1)-R(I+1))/MASS_DENSITY(I)
	    TC(I)=1.0D-10*TA(I)*LOG(TA(I+1)/TA(I-1))/LOG(R(I-1)/R(I+1))/R(I)/MASS_DENSITY(I)
!	    TC(I)=LOG(V(I+1)/V(I-1))/LOG(R(I+1)/R(I-1))
!	    TB(I)=1.0D-03*SQRT(BOLTZMANN_CONSTANT()*(POP_ATOM(I)+ED(I))*T(I)/MASS_DENSITY(I))
!	    TC(I)=TB(I)*TB(I)*TC(I)/R(I)
	  END DO
	  YV(1)=YV(2); YV(ND)=YV(ND-1)
	  TC(1)=TC(2); TC(ND)=TC(ND-1)
	  DO I=1,ND
	    TB(I)=1.0D-03*SQRT(BOLTZMANN_CONSTANT()*(POP_ATOM(I)+ED(I))*T(I)/MASS_DENSITY(I))
	    TB(I)=TB(I)*TB(I)*(SIGMA(I)+1.0D0)/R(I)
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  CALL DP_CURVE(ND,XV,TB)
	  CALL DP_CURVE(ND,XV,TC)
	  YAXIS='\gr\u-1\d|dPdR|'
!
	ELSE IF(XOPT .EQ. 'BEK')THEN
	  YV(1:ND)=0.5D+10*MASS_DENSITY(1:ND)*V(1:ND)*V(1:ND)
	  YV(1:ND)=4*PI*1.0D+30*YV(1:ND)*R(1:ND)*R(1:ND)
	  CALL LUM_FROM_ETA(YV,R,ND)
	  DO I=ND-1,1,-1
	    YV(I)=YV(I+1)+YV(I)
	  END DO	
	  EK_EJECTA=YV(1)
	  WRITE(6,*)RED_PEN
	  WRITE(6,'(A,ES10.3,A)')' Kinetic energy of ejecta is ',EK_EJECTA,' ergs'
	  YAXIS='Log EK(ergs\u \dcm\u-3)'
	  YV(1:ND-1)=DLOG10(YV(1:ND-1))
	  CALL DP_CURVE(ND-1,XV,YV)
!
	  T1=4.0D+30*PI/MASS_SUN()
	  DO I=1,ND
	    ZETA(I)=T1*MASS_DENSITY(I)*CLUMP_FAC(I)*R(I)*R(I)
	  END DO
	  CALL TORSCL(TA,ZETA,R,TB,TC,ND,METHOD,TYPE_ATM)
	  WRITE(6,'(A,F9.4,A)')' Mass of envelope (ejecta) is',TA(ND),' Msun'
	  T1=SQRT(2.0D0*EK_EJECTA/TA(ND)/MASS_SUN())/1.0D+05
	  WRITE(6,'(A,F10.2,A)')' SQRT(Mean square velocity) is ',T1,' km/s',DEF_PEN
!
!
!
	ELSE IF(XOPT .EQ. 'RONV')THEN
	  T1=1.0D+10/1.0D+05/24.0D0/3600.0D0
	  DO I=1,ND
	    YV(I)=V(I)/R(I)
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='t(days)'
!
	ELSE IF(XOPT .EQ. 'VEL')THEN
	  DO I=1,ND
	    YV(I)=V(I)
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='V(kms\u-1\d)'
!
	ELSE IF(XOPT .EQ. 'LOGV') THEN
	  CALL DLOGVEC(V,YV,ND)
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='Log V(kms\u-1\d) '
!
	ELSE IF(XOPT .EQ. 'VOL') THEN
	  CALL USR_OPTION(T1,'SN_AGE','1.0','SN age in days')
	  CALL USR_OPTION(T2,'dAGE','0.1','Fractional change in SN age')
	  T2=1.0D-05*T1*24.0D0*3600.0D0*T2
	  DO I=1,ND
	    TA(I)=(1.0D0+V(I)*T2/R(I))**2
	    TA(I)=TA(I)*(1.0D0+V(I)*T2/R(I)*(SIGMA(I)+1.0D0))
	    TA(I)=TA(I)**(1.0D0/3.0D0)
	  END DO
	  CALL DP_CURVE(ND,XV,TA)
	  YAXIS='Vf/Vi**0.33'
!
	ELSE IF(XOPT .EQ. 'RZERO') THEN
	  CALL USR_OPTION(T1,'SN_AGE','1.0','SN age in days')
	  T2=1.0D-05*T1*24.0D0*3600.0D0
	  DO I=1,ND
	    TA(I)=R(I)-V(I)*T2
	  END DO
	  CALL DP_CURVE(ND,XV,TA)
	  YAXIS='Rzero'
!
! Allows various velocity parameters to be varied to test their
! effect on the velocity law.
!
	ELSE IF(XOPT .EQ. 'TSTV')THEN
	  CALL USR_OPTION(SCLHT,'SCLHT','0.01D0','DEL R/R*')
	  CALL USR_OPTION(VCORE,'VCORE','1.0D0',' ')
	  CALL USR_OPTION(VPHOT,'VPHOT','100.0D0',' ')
	  CALL USR_OPTION(V_BETA1,'BETA1','1.0D0',' ')
	  CALL USR_OPTION(V_EPS1,'EPS1','1.0D0',' ')
	  CALL USR_OPTION(VINF1,'VINF1','1.0D0',' ')
	  CALL USR_OPTION(V_BETA2,'BETA2','1.0D0',' ')
	  CALL USR_OPTION(V_EPS2,'EPS2','1.0D0',' ')
	  DEFAULT=WR_STRING(VINF1)
	  CALL USR_OPTION(VINF2,'VINF2',DEFAULT,' ')
!
	  CALL USR_OPTION(TYPE,'PLOT','VEL','VEL, SIG, LOGV, LOGS')
!
	  CALL STARPCYG_TST(Z,DTAU,ZV,R(1),R(ND),
	1                   SCLHT,VCORE,VPHOT,VINF1,V_BETA1,V_EPS1,
	1                   VINF2,V_BETA2,V_EPS2,ND,TA,TB,TC,L_FALSE,LU_IN)
!
	  Z(1:ND)=DLOG10(Z(1:ND)/R(ND))
	  TYPE=UC( TRIM(TYPE) )
	  IF(TYPE .EQ. 'LOGV')THEN
	    CALL DLOGVEC(DTAU,YV,ND)
	    YAXIS='Log(V(kms\u-1\d))'
	  ELSE IF(TYPE .EQ. 'SIG')THEN
	    YV(1:ND)=ZV(1:ND)+1.0D0
	    YAXIS='\gs+1'
	  ELSE IF(TYPE .EQ. 'LOGS')THEN
	    YV(1:ND)=DLOG10(ZV(1:ND)+1.0D0)
	    YAXIS='LOG(\gs+1)'
	  ELSE
	    YV(1:ND)=DTAU(1:ND)
	    YAXIS='V(kms\u-1\d)'
	  END IF
	  CALL DP_CURVE(ND,Z,YV)
	  XAXSAV=XAXIS
	  XAXIS='Log(r/R\d*\u)'
!
! To convert from cm/s to km/s we cale the sound speed by a factor of 10^{-5}. However,
! T in units of 10^4 K, thus the scale factor is on 10^{-3}.
!
	ELSE IF(XOPT .EQ. 'ISOC')THEN
	  DO I=1,ND
	    YV(I)=1.0D-03*SQRT(BOLTZMANN_CONSTANT()*(POP_ATOM(I)+ED(I))*T(I)/MASS_DENSITY(I))
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='c(km/s)'
!
	ELSE IF(XOPT .EQ. 'SIGMA')THEN
	  DO I=1,ND
	    YV(I)=DLOG10(SIGMA(I)+1.0)
	  END DO
	  DO I=2,ND-1
	    TC(I)=DLOG10(LOG(V(I+1)/V(I-1))/LOG(R(I+1)/R(I-1)))
	  END DO
	  TC(1)=TC(2); TC(ND)=TC(ND-1)
	  CALL DP_CURVE(ND,XV,YV)
	  CALL DP_CURVE(ND,XV,TC)
	  YAXIS='Log(\gs+1)'
!
	ELSE IF(XOPT .EQ. 'FONR')THEN
	  IF(ROSS_MEAN(1) .NE. 0.0D0)THEN
	    DO I=1,ND
	      YV(I)=FLUX_MEAN(I)/ROSS_MEAN(I)
	    END DO
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='\gx(Flux)/\gx(Ross)'
	  ELSE
	    WRITE(6,*)'Error --- Rosseland mean opacity not defined'
	  END IF
!
	ELSE IF(XOPT .EQ. 'CAKT')THEN
	  CALL USR_HIDDEN(T1,'VTH','10.0','Thermal doppler velocity (km/s)')
	  DO I=1,ND
	    T2=(SIGMA(I)+1.0D0)*V(I)/R(I)
	    YV(I)=DLOG10(6.65D-15*ED(I)*T1/T2)
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='Log(t)'
!
! We subtract 1 since M(t), by definition, does not include the
! force on the elctrons. With the present definition, it does include
! bound-free and free contributins (only important at depth.)
!
	ELSE IF(XOPT .EQ. 'MT')THEN
	  T1=1.0D-30*LUM_SUN()*LUM/4.0D0/PI/C_CMS
	  T2=T1*6.65D-15
	  DO I=1,ND
	    YV(I)=FLUX_MEAN(I)/ED(I)/6.65D-15 - 1.0D0
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='M(t)'
!
! Sobolev radial optical depth scale. Assumes fg=1, POP_ATOM for the
! levelp population.
!
	ELSE IF(XOPT .EQ. 'TAUSOB')THEN
	  CALL USR_OPTION(FL,'LAM','1000','Wavelength in Ang')
	  FL=ANG_TO_HZ/FL
	  DO I=1,ND
	    T1=3.0D-10*OPLIN*POP_ATOM(I)*R(I)/V(I)/FL
	    YV(I)=LOG10(T1)
	    ZV(I)=LOG10(T1/(1.0D0+SIGMA(I)))
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  CALL DP_CURVE(ND,XV,ZV)
	  YAXIS='Log(\gt\dSob\u/fX)'
!
	ELSE IF(XOPT .EQ. 'YATOM')THEN
	  DO I=1,ND
	    YV(I)=LOG10(POP_ATOM(I))
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='Log(N\di\u)'
	ELSE IF(XOPT .EQ. 'YSPEC')THEN
	  FOUND=.FALSE.
	  ELEC=.FALSE.
	  FLAG=.FALSE.
!
	  WRITE(6,*)'Option plots fractional abundance (N/N[atom]), mass fraction, or species density.'
!
	  CALL USR_OPTION(ELEC,'FRAC','T','Fractional abundance')
	  IF(.NOT. ELEC)CALL USR_OPTION(FLAG,'MF','T','Mass fraction')
	  IF(ELEC)THEN
	    DO ISPEC=1,NSPEC
	      IF(XSPEC .EQ. SPECIES(ISPEC) .OR. (XSPEC .EQ. 'ALL' .AND.
	1              POPDUM(ND,ISPEC) .GT. 0.0D0))THEN
	        YV(1:ND)=LOG10(POPDUM(1:ND,ISPEC)/POP_ATOM(1:ND)+1.0D-100)
	        CALL DP_CURVE(ND,XV,YV)
	        FOUND=.TRUE.
	        IF(XSPEC .NE. 'ALL')EXIT
	      END IF
	    END DO
	    YAXIS='Log Fractional abundance (N\dX\u/N\dA\u)'
	  ELSE IF(FLAG)THEN
	    DO ISPEC=1,NSPEC
	      IF(XSPEC .EQ. SPECIES(ISPEC) .OR. (XSPEC .EQ. 'ALL' .AND.
	1              POPDUM(ND,ISPEC) .GT. 0.0D0))THEN
	        T1=AT_MASS(ISPEC)*ATOMIC_MASS_UNIT()
	        DO I=1,ND
	          YV(I)=T1*POPDUM(I,ISPEC)/MASS_DENSITY(I)
	        END DO
	        FOUND=.TRUE.
	        CALL DP_CURVE(ND,XV,YV)
	        IF(XSPEC .NE. 'ALL')EXIT
	      END IF
	    END DO
	    YAXIS='Mass fraction'
	  ELSE
	    DO ISPEC=1,NSPEC
	      IF(XSPEC .EQ. SPECIES(ISPEC) .OR. (XSPEC .EQ. 'ALL' .AND.
	1              POPDUM(ND,ISPEC) .GT. 0.0D0))THEN
	        YV(1:ND)=LOG10(POPDUM(1:ND,ISPEC)+1.0D-100)
	        CALL DP_CURVE(ND,XV,YV)
	        FOUND=.TRUE.
	        IF(XSPEC .NE. 'ALL')EXIT
	      END IF
	    END DO
	    YAXIS='Log Species density (cm\u-3\d)'
	  END IF
!
	  IF(XSPEC .EQ. 'ALL')THEN
	    CALL WR_SPEC_SUM_V2(ELEC,FLAG,XV,ND)
	  ELSE IF(.NOT. FOUND)THEN
	    WRITE(T_OUT,*)'Error --- unrecognized species'
	  END IF
!
	ELSE IF(XOPT .EQ. 'YION')THEN
	  DO I=1,ND
	    YV(I)=LOG10(POPION(I))
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='Log(N\di\u)'
	  WRITE(T_OUT,*)'Warning: this option plots ionized species only'
!
	ELSE IF(XOPT .EQ. 'CLUMP')THEN
	  
	  CALL USR_HIDDEN(ELEC,'RECIP','F','Plot 1 / fillingfactor')
	  IF(ELEC)THEN
	    YV(1:ND)=1.0D0/CLUMP_FAC(1:ND)
	    YAXIS='Recipricoal Filling Factor'
	  ELSE
	    YV(1:ND)=CLUMP_FAC(1:ND)
	    YAXIS='Filling Factor'
	  END IF
	  CALL DP_CURVE(ND,XV,YV)
	ELSE IF(XOPT .EQ. 'PROF')THEN
	  CALL PLT_PROFS()
!
!
!
	ELSE IF(XOPT .EQ. 'LNID')THEN
!
	  WRITE(T_OUT,'(A)')' '
	  WRITE(T_OUT,'(A)')' Create line ID file for stars with weak winds.'
	  WRITE(T_OUT,'(A)')' The list is created based on the line optical depth at a give Tau(e.s.)'
	  WRITE(T_OUT,'(A)')' Use XTAUC and TAUC option to find Tau(e.s) at Tau(cont)=2/3'
	  WRITE(T_OUT,'(A)')' Use PLNID for plane-parallel models'
	  WRITE(T_OUT,'(A)')' PLNID may also work better for O stars with weak winds'
	  WRITE(T_OUT,'(A)')' '
	  DO I=1,ND
	    ESEC(I)=6.65D-15*ED(I)
	  END DO
          CALL TORSCL(TA,ESEC,R,TB,TC,ND,METHOD,TYPE_ATM)
	  CALL USR_OPTION(TAU_VAL,'TAUES','0.67D0',
	1          'Electron scattering optical depth at whch line optical depth is evaluated')
	  CALL USR_OPTION(VTURB,'VTURB','10.0','Turbulent velcity for line profile (units km/s)?')
!
	  DO I=1,ND
	    J=I
	    IF(TA(I) .GT. TAU_VAL)THEN
	      EXIT
	    END IF
	  END DO
	  DPTH_INDX=J
!
	  DEFAULT='0.01'
	  CALL USR_OPTION(TAU_LIM,'TAU',DEFAULT,'We only output lines with Tau > ?')
	  CALL USR_OPTION(FLAG,'ADD','F','Add transition name to output file?')
!
	  WRITE(T_OUT,'(/,A,/)')' Values at depth just exceeeding requested TAU(e.s)'
          WRITE(T_OUT,'(1X,A,1P,E14.6)')'     R(I)/R*=',R(DPTH_INDX)/R(ND)
          WRITE(T_OUT,'(1X,A,1P,E14.6)')'        V(I)=',V(DPTH_INDX)
          WRITE(T_OUT,'(1X,A,1P,E14.6)')'        T(I)=',T(DPTH_INDX)
          WRITE(T_OUT,'(1X,A,1P,E14.6)')'       ED(I)=',ED(DPTH_INDX)
          WRITE(T_OUT,'(1X,A,1P,E14.6)')'TAU(e.s.)(I)=',TA(DPTH_INDX)
!
	  WRITE(T_OUT,'(A)')' '
	  DEFAULT=WR_STRING(LAM_ST)
	  CALL USR_OPTION(LAM_ST,'LAMST',DEFAULT,FREQ_INPUT)
	  DEFAULT=WR_STRING(LAM_EN)
	  CALL USR_OPTION(LAM_EN,'LAMEN',DEFAULT,FREQ_INPUT)
	  DEFAULT='LINE_ID'
	  CALL USR_OPTION(FILENAME,'FILE',DEFAULT,'Output file for line IDs')
	  CALL GEN_ASCI_OPEN(73,FILENAME,'UNKNOWN',' ',' ',IZERO,IOS)
!
! MNL_F (MNUP_F) denotes the lower (upper) level in the full atom.
!
	  CNT=0
!
! Determine index range encompassing desired wavelength range.
!
	  NU_ST=ANG_TO_HZ/LAM_ST
	  IF(NU_ST .LT. VEC_FREQ(1))THEN
	    NL=GET_INDX_DP(NU_ST,VEC_FREQ,N_LINE_FREQ)
	    IF(VEC_FREQ(NL) .GT. NU_ST)NL=NL+1
	  ELSE
	    NL=1
	  END IF
!
	  NU_EN=ANG_TO_HZ/LAM_EN
	  IF(NU_EN .GT. VEC_FREQ(N_LINE_FREQ))THEN
	    NUP=GET_INDX_DP(NU_EN,VEC_FREQ,N_LINE_FREQ)
	    IF(VEC_FREQ(NUP) .LT. NU_EN)NUP=NUP-1
	  ELSE
	    NUP=N_LINE_FREQ
	  END IF
!
! We use T4 for the interpolation to the requested optical depth.
!
	  T4=(TAU_VAL-TA(DPTH_INDX-1))/(TA(DPTH_INDX)-TA(DPTH_INDX-1))
	  DO LINE_INDX=NL,NUP
	    MNL_F=VEC_MNL_F(LINE_INDX)
	    MNUP_F=VEC_MNUP_F(LINE_INDX)
	    FL=VEC_FREQ(LINE_INDX)
!
	    IF(ANG_TO_HZ/FL .GT. LAM_ST .AND. ANG_TO_HZ/FL .LT. LAM_EN)THEN
	      DO ID=1,NUM_IONS
	        TAU_SOB=0.0D0
	        IF(VEC_SPEC(LINE_INDX) .EQ. ION_ID(ID) .AND.
	1           (XSPEC .EQ. ' ' .OR. XSPEC .EQ. UC(ION_ID(ID))) )THEN
	          DO I=1,DPTH_INDX
	            T1=ATM(ID)%W_XzV_F(MNUP_F,I)/ATM(ID)%W_XzV_F(MNL_F,I)
	            CHIL(I)=OPLIN*VEC_OSCIL(LINE_INDX)*( T1*ATM(ID)%XzV_F(MNL_F,I)-
	1               ATM(ID)%GXzV_F(MNL_F)*ATM(ID)%XzV_F(MNUP_F,I)/ATM(ID)%GXzV_F(MNUP_F) )
	            T2=12.86D0*SQRT( T(I)/AT_MASS(SPECIES_LNK(ID))+(VTURB/12.86D0)**2 )/C_KMS
	            CHIL(I)=MAX(1.0D-15*CHIL(I)/FL/T2/SQRT(PI),1.0D-10)
	          END DO
                  CALL TORSCL_V2(TA,CHIL,R,TB,TC,DPTH_INDX,'ZERO',TYPE_ATM,L_FALSE)
	          TAU_SOB=T4*TA(DPTH_INDX)+(1.0D0-T4)*TA(DPTH_INDX-1)
	          STRING=' ';  IF(FLAG)STRING=VEC_TRANS_NAME(LINE_INDX)
	          IF(TAU_SOB .GT. TAU_LIM)THEN
	            WRITE(73,'(A,3ES14.5,3X,F3.0,I6,5X,A)')ION_ID(ID),ANG_TO_HZ/FL,
	1                                TAU_SOB,ANG_TO_HZ/FL,1.0,2,TRIM(STRING)
	          END IF
	        END IF
	      END DO
	    END IF
	  END DO
	  CLOSE(UNIT=73)
! 
!
! Designed for use with plane-parallel models. We compute the central intesnity of
! each line to deduce line-identifcations.
!
	ELSE IF(XOPT .EQ. 'PLNID')THEN
!
	  WRITE(T_OUT,'(A)')' '
	  WRITE(T_OUT,'(A)')' Create line ID file for plane-paraillel models'
	  WRITE(T_OUT,'(A)')' The list is created based on an estimate of the central line intensity'
	  WRITE(T_OUT,'(A)')' The procedure may also work for most O stars -- especially those with weak winds'
	  WRITE(T_OUT,'(A)')' '
	  CALL USR_OPTION(VTURB,'VTURB','10.0','Turbulent velcity for line profile (units km/s)?')
!
	  DEFAULT='0.01'
	  CALL USR_OPTION(DEPTH_LIM,'DEPTH',DEFAULT,
	1        'Output lines whose central depth from the continuum is > LIMIT')
	  CALL USR_OPTION(FLAG,'ADD','F','Add transition name to output file?')
!
	  WRITE(T_OUT,'(A)')' '
	  DEFAULT=WR_STRING(LAM_ST)
	  CALL USR_OPTION(LAM_ST,'LAMST',DEFAULT,FREQ_INPUT)
	  DEFAULT=WR_STRING(LAM_EN)
	  CALL USR_OPTION(LAM_EN,'LAMEN',DEFAULT,FREQ_INPUT)
	  DEFAULT='LINE_ID'
	  CALL USR_OPTION(FILENAME,'FILE',DEFAULT,'Output file for line IDs')
	  CALL GEN_ASCI_OPEN(73,FILENAME,'UNKNOWN',' ',' ',IZERO,IOS)
	  CALL USR_HIDDEN(FLUSH_FILE,'FLUSH','F','Flush outout immediately?')
!
! MNL_F (MNUP_F) denotes the lower (upper) level in the full atom.
!
! Determine index range encompassing desired wavelength range.
!
	  NU_ST=ANG_TO_HZ/LAM_ST
	  IF(NU_ST .LT. VEC_FREQ(1))THEN
	    NL=GET_INDX_DP(NU_ST,VEC_FREQ,N_LINE_FREQ)
	    IF(VEC_FREQ(NL) .GT. NU_ST)NL=NL+1
	  ELSE
	    NL=1
	  END IF
!
	  NU_EN=ANG_TO_HZ/LAM_EN
	  IF(NU_EN .GT. VEC_FREQ(N_LINE_FREQ))THEN
	    NUP=GET_INDX_DP(NU_EN,VEC_FREQ,N_LINE_FREQ)
	    IF(VEC_FREQ(NUP) .LT. NU_EN)NUP=NUP-1
	  ELSE
	    NUP=N_LINE_FREQ
	  END IF
	  WRITE(6,*)' '
	  WRITE(6,*)'Number of line transitions in the interval is',NUP-NL+1
!
	  FL_SAVE=0.0D0
	  DO LINE_INDX=NL,NUP
	    MNL_F=VEC_MNL_F(LINE_INDX)
	    MNUP_F=VEC_MNUP_F(LINE_INDX)
	    FL=VEC_FREQ(LINE_INDX)
!
! Since the opacities are sequentially computed, we dont need to compute them 
! at all frequencies. At present we simply compute them every 1000 km/s.
! This saves considerable computing time.
!
	    IF( 2.998D+05*ABS(FL_SAVE-FL)/FL .GT. 1000.0D0)THEN
	      FL_SAVE=FL
	      INCLUDE 'OPACITIES.INC'
!
! Compute DBB and DDBBDT for diffusion approximation. DBB=dB/dR
! and DDBBDT= dB/dTR .
!
	      T1=HDKT*FL/T(ND)
	      T2=1.0D0-EMHNUKT(ND)
	      DBB=TWOHCSQ*( FL**3 )*T1*DTDR/T(ND)*EMHNUKT(ND)/(T2**2)
	      IF(.NOT. DIF)DBB=0.0D0
!
! Adjust the opacities and emissivities for the influence of clumping.
!
	      DO I=1,ND
	        ETA(I)=ETA(I)*CLUMP_FAC(I)
	        CHI(I)=CHI(I)*CLUMP_FAC(I)
	        ESEC(I)=ESEC(I)*CLUMP_FAC(I)
	        CHI_RAY(I)=CHI_RAY(I)*CLUMP_FAC(I)
	        CHI_SCAT(I)=CHI_SCAT(I)*CLUMP_FAC(I)
	      END DO
	    END IF
!
! ID must be set here as it modified by OPACITIES.INC
!
	    ID=VEC_ION_INDX(LINE_INDX)
	    GLDGU=ATM(ID)%GXzV_F(MNL_F)/ATM(ID)%GXzV_F(MNUP_F)
	    DO I=1,ND
	      T1=ATM(ID)%W_XzV_F(MNUP_F,I)/ATM(ID)%W_XzV_F(MNL_F,I)
	      CHIL(I)=OPLIN*ATM(ID)%AXzV_F(MNL_F,MNUP_F)*(T1*ATM(ID)%XzV_F(MNL_F,I)-GLDGU*ATM(ID)%XzV_F(MNUP_F,I) )
	      ETAL(I)=EMLIN*FREQ*ATM(ID)%AXzV_F(MNUP_F,MNL_F)*ATM(ID)%XzV_F(MNUP_F,I)
	      T2=12.86D0*SQRT( T(I)/AT_MASS(SPECIES_LNK(ID))+(VTURB/12.86D0)**2 )/C_KMS
	      T2=1.0D-15/FL/T2/SQRT(PI)
	      CHIL(I)=MAX(T2*CHIL(I)*CLUMP_FAC(I),1.0D-10)
	      ETAL(I)=ETAL(I)*T2*CLUMP_FAC(I)
	    END DO
	    CALL GET_FLUX_DEFICIT(T2,R,ETA,CHI,CHI_RAY,CHI_SCAT,ESEC,ETAL,CHIL,FL,DBB,ND)
!
	    IF(ABS(T2).GT. DEPTH_LIM)THEN
	      STRING=' ';  IF(FLAG)STRING=VEC_TRANS_NAME(LINE_INDX)
	      T1=ANG_TO_HZ/FL
	      IF(.NOT. ATM(ID)%OBSERVED_LEVEL(MNL_F) .OR. .NOT. ATM(ID)%OBSERVED_LEVEL(MNUP_F))T1=-T1
	      WRITE(73,'(A,3ES14.5,3X,F3.0,I6,5X,A)')ION_ID(ID),T1,
	1                               T2,ANG_TO_HZ/FL,1.0,2,TRIM(STRING)
	      IF(FLUSH_FILE)FLUSH(UNIT=73)
	    END IF
	    IF(MOD(LINE_INDX-NL,MAX(10,(NUP-NL)/10)) .EQ. 0)THEN
	      WRITE(6,'(A,I7,A,I7,A)')' Done',LINE_INDX-NL+1,' of ',NUP-NL+1,' lines'
	    END IF
	  END DO
	  CLOSE(UNIT=73)
!
! Section to allow various line optical deph and wavelength distributions to be examined.
!
! 1. POW: Plots the # of lines in Log(Tau) or Log(Line strength) space. Each uniform logarithmic bin
!         is 0.25 wide in LOG10 space. For Log(tau) space, the line optical depth is normalized by
!         the elctron scatering opacity.
!
! 2. DIST: PLots the the ilogaritmic Sobolev optical depth (radial, or angle averaged) for every 
!          line in a given wavelenth interval as a function of the line wavelngth. Lines with
!          a negative optical depth are set to -20.
!
! 3. NV:   Plots the  number of lines with optical depth > TAU_MIN in a given velocity interval.
! 
	ELSE IF(XOPT .EQ. 'DIST' 
	1        .OR. XOPT .EQ. 'NV' 
	1        .OR. XOPT .EQ. 'POW')THEN
!
! Used to count the number lines as a function of optical depth over the entire
! wavelngth interval.
!
	  TAU_GRT_LOGX(:)=0.0D0
!
	  IF(DPTH_INDX .LT. 1 .OR. DPTH_INDX .GT. ND)DPTH_INDX=ND/2
	  DEFAULT=WR_STRING(DPTH_INDX)
	  CALL USR_OPTION(DPTH_INDX,'DPTH',DEFAULT,'Depth for line plot')
!
! Ouput summary of model data at depth point
!
          WRITE(T_OUT,'(1X,A,1P,E14.6)')'    R(I)/R*=',R(DPTH_INDX)/R(ND)
          WRITE(T_OUT,'(1X,A,1P,E14.6)')'       V(I)=',V(DPTH_INDX)
          WRITE(T_OUT,'(1X,A,1P,E14.6)')'       T(I)=',T(DPTH_INDX)
          WRITE(T_OUT,'(1X,A,1P,E14.6)')'      ED(I)=',ED(DPTH_INDX)
          TA(1:ND)=5.65D-15*ED(1:ND)
          CALL TORSCL(DTAU,TA,R,TB,TC,ND,METHOD,TYPE_ATM)
          WRITE(T_OUT,'(1X,A,1P,E14.6)')'  TAU_ES(I)=',DTAU(DPTH_INDX)
!
	  IF(XOPT .EQ. 'POW')THEN
	    CALL USR_OPTION(DO_TAU,'TAU','T','Plot versus TAU or Line strength')
	  ELSE IF(XOPT .EQ. 'DIST')THEN
	    CALL USR_OPTION(DO_TAU,'SRCE','F','Plot source function instead of tau')
	  END IF
	  CALL USR_HIDDEN(WEIGHT_NV,'WGT','T','Weight Tau by [1.0D0-E(-Tau)]?')
	  CALL USR_HIDDEN(MEAN_TAU,'MEAN','F','Use MEAN instead of radial Tau?')
!
	  DEFAULT=WR_STRING(LAM_ST)
	  CALL USR_OPTION(LAM_ST,'LAMST',DEFAULT,FREQ_INPUT)
	  DEFAULT=WR_STRING(LAM_EN)
	  CALL USR_OPTION(LAM_EN,'LAMEN',DEFAULT,FREQ_INPUT)
!
! This makes sure LAM_ST,LAM_END are in A.
!
	  IF(KEV_INPUT)THEN
	    LAM_ST=ANG_TO_HZ/(LAM_ST*KEV_TO_HZ)
	    LAM_EN=ANG_TO_HZ/(LAM_EN*KEV_TO_HZ)
	  ELSE IF(HZ_INPUT)THEN
	    LAM_ST=ANG_TO_HZ/LAM_ST
	    LAM_EN=ANG_TO_HZ/LAM_EN
	  END IF
	  T1=MIN(LAM_ST,LAM_EN); LAM_EN=MAX(LAM_ST,LAM_EN); LAM_ST=T1
!
! TAU_CONSTANT is used to compute the optical depth. Two choices are
! availabe:
!    (1) Radial optical depth. Valid in outer region.
!    (2) Angle avearged optical depth. Might be more appropriated at depth.
!           For an O star it would be better to average over the disk of the
!           star. For a W-R star the disk is not well defined.
! 
	  I=DPTH_INDX
	  TAU_CONSTANT=1.0D-15*C_KMS*R(I)/V(I)
	  IF(MEAN_TAU)THEN
	    T1=SIGMA(DPTH_INDX)
	    T2=SQRT(ABS(T1))
	    IF(T2 .LT. 0.01)THEN
	      TAU_CONSTANT=TAU_CONSTANT*(1.0D0-T1/3.0D0+T1*T1/5.0D0)
	    ELSE IF(T1 .LT. 0)THEN
	      TAU_CONSTANT=0.5D0*TAU_CONSTANT*LOG( (1.0D0+T2)/(1.0D0-T2) )/T2
	    ELSE
	      TAU_CONSTANT=TAU_CONSTANT*ATAN(T2)/T2
	    END IF
	  ELSE
	    TAU_CONSTANT=TAU_CONSTANT/(1.0D0+SIGMA(I))
	  END IF
!
! NB: DELTA_TAU is the spacing in LOG10(Tau)
!
	  IF(XOPT .EQ. 'POW')THEN
	    DELTA_TAU=0.25D0
	    TAU_BEG=-5.0D0
	    NBINS=(6.0D0-TAU_BEG)/DELTA_TAU+1
	    DO I=1,NBINS
 	      ZV(I)=TAU_BEG+(I-1)*DELTA_TAU
	    END DO
	    YV(:)=0.0D0
	  ELSE IF(XOPT .EQ. 'NV')THEN
	    DELTA_TAU=V(DPTH_INDX)/C_KMS
	    NBINS=LOG(LAM_EN/LAM_ST)/LOG(1.0D0+DELTA_TAU)+1
	    DO I=1,NBINS
	      ZV(I)=LAM_ST*(1.0D0+DELTA_TAU)**(I-1)
	    END DO
	    YV(:)=0.0D0
	    CALL USR_OPTION(TAU_MIN,'TAU_MIN','1.0D0','Minimum Tau')
	    WRITE(T_OUT,'(1X,A,1P,E11.4)')'Velocity is:',V(DPTH_INDX),' km/s'
	  END IF
!
! MNL_F (MNUP_F) denotes the lower (upper) level in the full atom.
!
	  CNT=0
!
! Determine index range encompassing desired wavelength range.
!
	  NU_ST=ANG_TO_HZ/LAM_ST
	  IF(NU_ST .LT. VEC_FREQ(1))THEN
	    NL=GET_INDX_DP(NU_ST,VEC_FREQ,N_LINE_FREQ)
	    IF(VEC_FREQ(NL) .GT. NU_ST)NL=NL+1
	  ELSE
	    NL=1
	  END IF
!
	  NU_EN=ANG_TO_HZ/LAM_EN
	  IF(NU_EN .GT. VEC_FREQ(N_LINE_FREQ))THEN
	    NUP=GET_INDX_DP(NU_EN,VEC_FREQ,N_LINE_FREQ)
	    IF(VEC_FREQ(NUP) .LT. NU_EN)NUP=NUP-1
	  ELSE
	    NUP=N_LINE_FREQ
	  END IF
!
	  DO LINE_INDX=NL,NUP
	    MNL_F=VEC_MNL_F(LINE_INDX)
	    MNUP_F=VEC_MNUP_F(LINE_INDX)
	    FL=VEC_FREQ(LINE_INDX)
	    I=DPTH_INDX
!
	    IF(ANG_TO_HZ/FL .GT. LAM_ST .AND. ANG_TO_HZ/FL .LT. LAM_EN)THEN
	      DONE_LINE=.FALSE.
	      DO ID=1,NUM_IONS
	        IF(VEC_SPEC(LINE_INDX) .EQ. ION_ID(ID) .AND.
	1           (XSPEC .EQ. ' ' .OR. XSPEC .EQ. UC(ION_ID(ID))) )THEN
	          T1=ATM(ID)%W_XzV_F(MNUP_F,I)/ATM(ID)%W_XzV_F(MNL_F,I)
	          CHIL(I)=OPLIN*VEC_OSCIL(LINE_INDX)*(
	1            T1*ATM(ID)%XzV_F(MNL_F,I)-
	1            ATM(ID)%GXzV_F(MNL_F)*ATM(ID)%XzV_F(MNUP_F,I)/
	1            ATM(ID)%GXzV_F(MNUP_F) )
	          ETAL(I)=EMLIN*FL*ATM(ID)%AXzV_F(MNUP_F,MNL_F)*ATM(ID)%XzV_F(MNUP_F,I)
	          DONE_LINE=.TRUE.
	        END IF
	      END DO
!
	      IF(DONE_LINE)THEN
	        IF(XOPT .EQ. 'DIST' .OR. XOPT .EQ. 'NV' .OR. DO_TAU)THEN
	          TAU_SOB=TAU_CONSTANT*CHIL(I)/FL
	          DO I=ITAU_GRT_LIM,-ITAU_GRT_LIM,-1
	             IF(TAU_SOB .GT. 10.0**I)THEN
	               TAU_GRT_LOGX(I)=TAU_GRT_LOGX(I)+1
	               EXIT
	             END IF
	          END DO
	          IF(XOPT .EQ. 'POW' .AND. TAU_SOB .GT. 0)THEN
                    I=(DLOG10(TAU_SOB)-TAU_BEG)/DELTA_TAU+1
	            IF(I .GT. 0 .AND. I .LE. NBINS)YV(I)=YV(I)+1
	          END IF
	        ELSE
!
! For simplicty, we choos a Doppler velocity of 10 km/s and ignore the
! factor of PI.
!
	          TAU_SOB=1.0D-15*CHIL(I)/(FL/2.998D+04)/(6.65D-15*ED(I))
	          IF(TAU_SOB .GT. 0)THEN
                    I=(DLOG10(TAU_SOB)-TAU_BEG+0.5D0*DELTA_TAU)/DELTA_TAU+1
	            IF(I .GT. 0 .AND. I .LE. NBINS)YV(I)=YV(I)+1
	            IF(I .GE. NBINS)YV(NBINS)=YV(NBINS)+1
	          END IF
	          DO I=ITAU_GRT_LIM,-ITAU_GRT_LIM,-1
	             IF(TAU_SOB .GT. 10.0**I)THEN
	               TAU_GRT_LOGX(I)=TAU_GRT_LOGX(I)+1
	               EXIT
	             END IF
	          END DO
	        END IF
!
	        CNT=CNT+1
	        IF(XOPT .EQ. 'DIST' .AND. DO_TAU)THEN
	          ZV(2*CNT-1)=0.2998E+04/FL		!Angstroms
	          ZV(2*CNT)=0.2998E+04/FL
	          YV(2*CNT-1)=-30.0D0
	          T1=ETAL(DPTH_INDX)/CHIL(DPTH_INDX)
	          IF(T1 .GT. 1.0D-30)THEN
	            YV(2*CNT)=DLOG10(T1)
	          ELSE IF(T1 .GE. 0.0D0)THEN
	            YV(2*CNT)=-30.0D0
	          ELSE IF(T1 .GT. -1.0D-30)THEN
	            YV(2*CNT)=-30.0D0
	          ELSE
	            YV(2*CNT)=-60.0D0-DLOG10(-T1)
	          END IF
	          TAU_SOB=1.0D-15*CHIL(DPTH_INDX)/(FL/2.998D+04)/(6.65D-15*ED(DPTH_INDX))
	          WRITE(66,'(6ES14.4)')FL,0.29979D+04/FL,TAU_SOB,T1,YV(2*CNT-1),YV(2*CNT)
	        ELSE IF(XOPT .EQ. 'DIST')THEN
	          ZV(2*CNT-1)=0.2998E+04/FL		!Angstroms
	          ZV(2*CNT)=0.2998E+04/FL
	          YV(2*CNT-1)=-10.0D0
	          IF(TAU_SOB .GT. 1.0D-10)THEN
	            YV(2*CNT)=DLOG10(TAU_SOB)
	          ELSE IF(TAU_SOB .GE. 0.0D0)THEN
	            YV(2*CNT)=-10.0D0
	          ELSE IF(TAU_SOB .GT. -1.0D-10)THEN
	            YV(2*CNT)=-10.0D0
	          ELSE
	            YV(2*CNT)=-10.0D0-DLOG10(-1.0D+10*TAU_SOB)
	          END IF
	        ELSE IF(XOPT .EQ. 'NV')THEN
	          IF(TAU_SOB .GT. TAU_MIN)THEN
	            T1=C_KMS/FL/100.0D0			!Angstroms
	            T2=LAM_ST*(1.0D0-0.5D0*DELTA_TAU)
	             I=LOG(T1/T2)/LOG(1.0D0+DELTA_TAU)+1
	            T1=0.0D0; T1=MAX(T1,TAU_SOB)
	            IF(WEIGHT_NV)THEN
	              YV(I)=YV(I)+(1.0-DEXP(-T1))
	            ELSE
	              YV(I)=YV(I)+1
	            END IF
	          END IF
	        END IF			!Specifi option
	      END IF			!Line included (e.g. for species)
	    END IF			!Line between LAM_ST and LAM_END
	  END DO			!Loop overlines
!
	  WRITE(T_OUT,'(1X,A,1PE12.4,A)')'ED:   ',ED(DPTH_INDX),'/cm^3'
	  WRITE(T_OUT,'(1X,A,1PE12.4,A)')'Vel:  ',V(DPTH_INDX),'km/s'
	  IF(XOPT .EQ. 'POW' .AND. .NOT. DO_TAU)THEN
	    WRITE(T_OUT,'(1X,A,F9.4)')'/\Line_strength=',DELTA_TAU
	    DO I=ITAU_GRT_LIM,-ITAU_GRT_LIM,-1
	      T1=10.0D0**I
	      WRITE(T_OUT,'(1X,A,I2,A,I9)')'CHIL(Vo)/CHI_ES > 10**(',I,'):',TAU_GRT_LOGX(I)
	    END DO
	  ELSE
	    WRITE(T_OUT,'(1X,A,F9.4)')'/\Tau=',DELTA_TAU
	    WRITE(T_OUT,'(A)')' '
	    WRITE(T_OUT,'(1X,A,I7)')' Total number of lines in interval is',NUP-NL+1
	    WRITE(T_OUT,'(1X,A)')' Number of lines in each decade of optical depth (NOT cummualtive).'
	    WRITE(T_OUT,'(A)')' '
	    DO I=ITAU_GRT_LIM,-ITAU_GRT_LIM,-1
	      T1=10.0D0**I
	      WRITE(T_OUT,'(1X,A,I2,A,I9)')'Tau > 10**(',I,'):',TAU_GRT_LOGX(I)
	    END DO
	  END IF
!
	  IF(XOPT .EQ. 'DIST')THEN
	    I=2*CNT
	    IF(HZ_INPUT)THEN
	      DO J=1,I
	         ZV(J)=ANG_TO_HZ/ZV(J)
	      END DO
	    ELSE IF(KEV_INPUT)THEN
	      DO J=1,I
	         ZV(J)=ANG_TO_HZ/ZV(J)/KEV_TO_HZ
	      END DO
	    END IF
	    CALL DP_CURVE(I,ZV,YV)
	    CALL GRAMON_PGPLOT('\gl(\A)','Log \gt\dSob\u',' ',' ')
	  ELSE IF(XOPT .EQ. 'NV')THEN
	    CALL DP_CURVE(NBINS,ZV,YV)
	    CALL GRAMON_PGPLOT('\gl(\A)','N/V',' ',' ')
	  ELSE
!
! We divide by DELTA_TAU so that dN/dLOG_TAU * DELTA_TAU gives the number
! of lines in that bin. The bin is cenetered on the log(tau)=ZV(I).
!
	    DO I=1,NBINS
	      IF(YV(I) .NE. 0)THEN
	        YV(I)=DLOG10(YV(I)/DELTA_TAU)
	      ELSE
	        YV(I)=-10.0D0
	      END IF
	    END DO
	    XAXIS='Log \gt'
	    YAXIS='Log(dN/dLog\gt)'
	    CALL DP_CURVE(NBINS,ZV,YV)
!
! Draw the canonical slope of 2/3 on the plot.
!
            I=-TAU_BEG/DELTA_TAU+1
	    T1=YV(I)
	    DO I=1,NBINS
	      YV(I)=T1-0.33333D0*ZV(I)
	    END DO
	    CALL DP_CURVE(NBINS,ZV,YV)
	  END IF
!
! 
!
	ELSE IF(XOPT .EQ. 'GF')THEN
!
	   I=0
	   DO LINE_INDX=1,N_LINE_FREQ
	     MNL_F=VEC_MNL_F(LINE_INDX)
	     MNUP_F=VEC_MNUP_F(LINE_INDX)
	     FL=VEC_FREQ(LINE_INDX)
!
	     T1=0
	     IF( XSPEC .EQ. UC(VEC_SPEC(LINE_INDX)) )THEN
	       DO ID=1,NUM_IONS
	         IF(XSPEC .EQ. UC(ION_ID(ID)))THEN
	           T1=ABS(VEC_OSCIL(LINE_INDX)*ATM(ID)%GXzV_F(MNL_F))
	         END IF
	       END DO
	     END IF
!
	     IF(T1 .GT. 0)THEN
	       I=I+1
	       ZV(2*I-1)=0.2998E+04/FL			!Angstroms
	       ZV(2*I)=0.2998E+04/FL
	       YV(2*I-1)=-10.0D0
	       YV(2*I)=DLOG10(T1)
	     END IF
	   END DO
!
	   IF(I .EQ. 0)THEN
	     WRITE(T_OUT,*)'Error no matching ioization stage'
	   ELSE
             WRITE(T_OUT,*)'Number of transitions found is',I
	     I=2*I
	     CALL DP_CURVE(I,ZV,YV)
	   END IF
!
! 
	ELSE IF(XOPT .EQ. 'CJBB')THEN
!
! RJ has previously been computed.
!
	  WRITE(6,*)' '
	  WRITE(6,*)'Option plots J/B (def), J/S or J versus X(def) or Log(Tau)'
!
	  CALL USR_OPTION(NORAD,'TPLT','F','Plot against Tau (alt is current X axis)?')
	  CALL USR_OPTION(JONS,'JONS','F','Plot J/S?')
	  JONLY=.FALSE.
	  IF(.NOT. JONS)THEN
	    CALL USR_OPTION(JONLY,'JONLY','F','Plot J only?')
	  END IF
	  IF(.NOT. JONS .AND. .NOT. JONLY)THEN
	    CALL USR_OPTION(SONLY,'SONLY','F','Plot S only?')
	  END IF
!
	  IF(JONS)THEN
	    DO I=1,ND
	      YV(I)=DLOG10( RJ(I)*CHI(I)/ETA_WITH_ES(I) )
	    END DO
	    YAXIS='J\dc\u/S'
	  ELSE IF(JONLY)THEN
	    DO I=1,ND
	      YV(I)=DLOG10(RJ(I))
	    END DO
	    YAXIS='J\dc\u'
	  ELSE IF(SONLY)THEN
	    DO I=1,ND
	      YV(I)=DLOG10(ETA_WITH_ES(I)/CHI(I))
	    END DO
	    YAXIS='S'
	  ELSE
	    DO I=1,ND
	      YV(I)=DLOG10( RJ(I)*
	1          ( EXP(HDKT*FREQ/T(I))-1.0D0 )/TWOHCSQ/(FREQ**3)  )
	    END DO
	    YAXIS='J\dc\u/B(T\de\u)'
	  END IF
!
	  IF(NORAD)THEN
	    IF(.NOT. ELEC)THEN
	      DO I=1,ND
	       CHI(I)=CHI(I)-ESEC(I)
	      END DO
	    END IF
	    CALL TORSCL(TA,CHI,R,TB,TC,ND,METHOD,TYPE_ATM)
	    DO I=1,ND
	      ZV(I)=DLOG10(TA(I))
	    END DO
	    CALL DP_CURVE(ND,ZV,YV)
	    XAXSAV=XAXIS
	    XAXIS='Log(\gt)'
	  ELSE
	    CALL DP_CURVE(ND,XV,YV)
	    WRITE(T_OUT,*)YV(ND-4),RJ(ND-4),CHI(ND-4),ETA(ND-4),ESEC(ND-4)
	    WRITE(T_OUT,*)YV(ND),RJ(ND),CHI(ND),ETA(ND),ESEC(ND)
	    OPEN(UNIT=37,FILE='RJ_COMP',STATUS='UNKNOWN')
	      WRITE(37,'(5X,A4,6X)')'  YV','  RJ',' CHI',' ETA','ESEC'
	      DO I=1,ND
	        WRITE(37,'(5E15.6)')YV(I),RJ(I),CHI(I),ETA(I),ESEC(I)
	      END DO
	    CLOSE(UNIT=37)
	  END IF
!
! 
! **************************************************************************
! **************************************************************************
!
! To enable a comparison of J (continuum) compared with that computed on a
! finer depth grid.
!
! **************************************************************************
! **************************************************************************
!
	ELSE IF(XOPT .EQ. 'JEXT'
	1  .OR. XOPT .EQ. 'EWEXT')THEN
!
! Compute J. DBB and S1 have been previously evaluated. VT is used
! for dCHI_dR.
!
	  CALL NEWJSOLD(TA,TB,TC,XM,WM,FB,RJ,DTAU,R,Z,P,
	1               ZETA,THETA,CHI,VT,JQW,
	1               THICK,DIF,DBB,IC,NC,ND,NP,METHOD)
!
	  S1=(ETA(1)+RJ(1)*CHI_SCAT(1))/CHI(1)
	  DO I=1,ND
	    SOURCE(I)=(ETA(I)+CHI_SCAT(I)*RJ(I))/CHI(I)
	  END DO
	  CALL NORDFLUX(TA,TB,TC,XM,DTAU,R,Z,P,SOURCE,CHI,dCHIdR,HQW,VT,
	1               S1,THICK,DIF,DBB,IC,NC,ND,NP,METHOD)
!
! Compute observed flux in Janskys for an object at 1 kpc .
!	(const=dex(23)*2*pi*dex(20)/(3.0856dex(21))**2 )
!
	  T1=6.599341D0*VT(1)*2.0D0		!2 DUE TO 0.5U
          WRITE(T_OUT,'('' Observed Flux is'',1Pe12.4,''Jy'')')T1
          WRITE(LU_LOG,'('' Observed Flux is'',1Pe12.4,''Jy'')')T1
!
	  CALL EXTEND_OPAC(CHIEXT,ETAEXT,ESECEXT,RJEXT,COEF,INDX,NDX
	1                     ,CHI,ETA,CHI_SCAT,RJ,ND)
!
	  DO I=1,NDX
	    ZETAEXT(I)=ETAEXT(I)/CHIEXT(I)
	    THETAEXT(I)=ESECEXT(I)/CHIEXT(I)
	    SOURCEEXT(I)=ZETAEXT(I)+THETAEXT(I)*RJEXT(I)
	    FOLD(I)=1.0D0/3.0D0
	  END DO
	  S1=SOURCE(1)
!
	  INACCURATE=.TRUE.
	  J=0				!Loop counter
	  DO WHILE(INACCURATE)
	    J=J+1
	    WRITE(T_OUT,'('' Beginning '',I3,''th loop'')')J
	    CALL FQCOMP_V2(TA,TB,TC,XM,DTAU,REXT,Z,PEXT,Q,FEXT,
	1          SOURCEEXT,CHIEXT,ETAEXT,JQWEXT,KQWEXT,DBB,HBC,
	1          INBC,IC,S1,THICK,INNER_BND_METH,NCX,NDX,NPX,METHOD)
	    T1=DBB
	    IF(INNER_BND_METH .EQ. 'ZERO_FLUX' .OR. INNER_BND_METH .EQ. 'HOLLOW')T1=0.0D0
	    CALL JFEAUNEW(TA,TB,TC,DTAU,REXT,RJEXT,Q,FEXT,
	1          ZETAEXT,THETAEXT,CHIEXT,T1,IC,HBC,
	1          INBC,THICK,DIF,NDX,METHOD)
	    S1=ZETAEXT(1)+THETAEXT(1)*RJEXT(1)
	    INACCURATE=.FALSE.
	    T1=0.0D0
	    DO I=1,NDX
	     T1=MAX(ABS(FOLD(I)-FEXT(I)),T1)
	     FOLD(I)=FEXT(I)
	     SOURCEEXT(I)=ZETAEXT(I)+THETAEXT(I)*RJEXT(I)
	    END DO
	    IF(T1 .GT. 1.0E-05)INACCURATE=.TRUE.
	    WRITE(T_OUT,'('' Maximum fractional change is'',1PE11.4)')T1
	  END DO
!
	  IF(XOPT .EQ. 'JEXT')THEN
!
! Will use AV for the new RJ
!
	    CALL UNGRID(AV,ND,RJEXT,NDX,GRID)
!
	    S1=(ETA(1)+AV(1)*CHI_SCAT(1))/CHI(1)
	    DO I=1,ND
	      SOURCE(I)=(ETA(I)+CHI_SCAT(I)*AV(I))/CHI(I)
	    END DO
	    CALL NORDFLUX(TA,TB,TC,XM,DTAU,R,Z,P,SOURCE,CHI,dCHIdR,HQW,VT,
	1                 S1,THICK,DIF,DBB,IC,NC,ND,NP,METHOD)
!
! Compute observed flux in Janskys for an object at 1 kpc .
!	(const=dex(23)*2*pi*dex(20)/(3.0856dex(21))**2 )
!
	    T1=6.599341D0*VT(1)*2.0D0		!2 DUE TO 0.5U
            WRITE(T_OUT,'('' Observed Flux[Ext] is'',1Pe12.4,''Jy'')')T1
            WRITE(LU_LOG,'('' Observed Flux[Ext] is'',1Pe12.4,''Jy'')')T1
!
! With this definition for YV, can never get an error worse than 200%.
! An error of 100% implies that RJ and RJEXT differ by a factor of 3.
!
	    DO I=1,ND
	      YV(I)=200.0D0*(AV(I)-RJ(I))/(AV(I)+RJ(I))
	    END DO
	    YAXIS='% error'
!
	    CALL USR_HIDDEN(ELEC,'LOGE','F','Log of error?')
	    IF(ELEC)THEN
	      DO I=1,ND
	        YV(I)=LOG10(YV(I))
	      END DO
	      YAXIS='Log(% error)'
	    END IF
!
	    CALL DP_CURVE(ND,XV,YV)
	  ELSE IF(XOPT .EQ. 'EWEXT')THEN
!
	    DO I=1,ND
	      SOURCE(I)=ZETA(I)+THETA(I)*RJEXT(GRID(I))
	    END DO
!
! TA is used for the line flux. Integral of TA dlog(r) is
! the line EW.
!
	    CALL SOBEW_GRAD_V2(SOURCE,CHI,CHI_SCAT,CHIL,ETAL,
	1              FORCE_MULT,LUM,
	1              V,SIGMA,R,P,JQW,HQW,TA,T1,S1,
	1              FREQ,'HOLLOW',DBB,IC,THICK,.FALSE.,NC,NP,ND,METHOD)
!
	    CALL DP_CURVE(ND,XV,FORCE_MULT)
	    YAXIS='Force Multiplier'
!
	    OPEN(UNIT=18,FILE='DSOB_FORCE_MULT',STATUS='UNKNOWN')
	      WRITE(18,'(3X,A1,10X,A1,15X,A1,13X,A1)')'I','R','V','M'
              DO I=1,ND
                WRITE(18,'(1X,I3,3X,3ES14.6)')I, R(I),V(I),FORCE_MULT(I)
              END DO
            CLOSE(UNIT=18)
!
! S1 is the continuum flux in Jy for an object at 1kpc.
! T1 is the line equivalent width in Angstroms.
! T3 is the line flux in ergs/cm^2/s
!
	    T2=LAMVACAIR(FREQ)		!Wavelength(Angstroms)
	    T3=T1*S1*1.0D-23*FREQ*1.0D+15/T2
	    WRITE(LU_NET,40008)T1,S1,T2,T3
	    WRITE(T_OUT,'(A)')RED_PEN
	    WRITE(T_OUT,40008)T1,S1,T2,T3
40008	    FORMAT(1X,'EW =',ES10.3,' Ang',5X,'I =',ES10.3,' Jy',5X,
	1             'Lambda =',ES11.4,' Ang',5X,'Line flux=',ES10.3,' ergs/cm^2/s')
	    WRITE(T_OUT,'(A)')DEF_PEN
	  END IF
!
! 
! **************************************************************************
! **************************************************************************
!
! Comparison of line source function with Jc and B.
!
! **************************************************************************
! **************************************************************************
!
	ELSE IF(XOPT .EQ. 'SRCEBB')THEN
!
! Assumes line opacity and emissivity have been computed in the set up.
!
	  DO I=1,ND
	    IF(CHIL(I).LT. 1.0D-30)THEN
	      YV(I)=-30
	    ELSE
	      YV(I)=DLOG10(  ETAL(I)/CHIL(I)*
	1          ( EXP(HDKT*FREQ/T(I))-1.0D0 )/TWOHCSQ/(FREQ**3)  )
	    END IF
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='S/B(T\de\u)'
!
	ELSE IF(XOPT .EQ. 'SRCEJC')THEN
!
! RJ has previously been computed, as have the Line opacity and emissivity.
!
	  DO I=1,ND
	    IF(CHIL(I).LT. 1.0D-30)THEN
	      YV(I)=-30
	    ELSE
	      YV(I)=DLOG10( ETAL(I)/CHIL(I)/RJ(I) )
	    END IF
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='S/J\dc\u'
!
	ELSE IF(XOPT .EQ. 'SRCE')THEN
!
! Assumes line opacity and emissivity have been computed in the set up.
!
	  CALL USR_OPTION(ELEC,'LIN','T','Plot S instead of Log S')
	  IF(ELEC)THEN
	    DO I=1,ND
	      IF( ABS(CHIL(I)) .LT. 1.0D-50)THEN
	         YV(I)=-30
	      ELSE
	        YV(I)=ETAL(I)/CHIL(I)
	      END IF
	    END DO
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='S'
	  ELSE
	    DO I=1,ND
	      IF(CHIL(I).LT. 1.0D-30)THEN
	         YV(I)=-30
	      ELSE
	        YV(I)=DLOG10( ETAL(I)/CHIL(I) )
	      END IF
	    END DO
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='Log S'
	  END IF
! 
!

	ELSE IF(XOPT .EQ. 'DION')THEN
	  TYPE=' '
	  DO ID=1,NUM_IONS
	    IF(XSPEC .EQ. UC(TRIM(ION_ID(ID))) ) THEN
	      CALL DLOGVEC(ATM(ID)%DXzV_F,YV,ND)
	      CALL DP_CURVE(ND,XV,YV)	
	      IF(ATM(ID)%ZXzV .LE. 9)WRITE(TYPE,'(I1)')NINT(ATM(ID)%ZXzV)
	      IF(ATM(ID)%ZXzV .GE. 10)WRITE(TYPE,'(I2)')NINT(ATM(ID)%ZXzV)
	      YAXIS='Log(H\u'//TRIM(TYPE)//'\d+)'
	      YV(1:ND)=ATM(ID)%DXzV(1:ND)
	    END IF
	  END DO
!
! 
!
	ELSE IF(XOPT .EQ. 'SCL')THEN
	  CALL USR_OPTION(T1,'SCL_FAC',' ',
	1   'Factor to multiply POPULATIONS by (<1)')
	  CALL USR_OPTION(K,'DEPTH',' ',
	1   'Depth index below which [for 1 to I] populations are reduced.')
!
	  DO ISPEC=1,NSPEC
	    IF(XSPEC .EQ. SPECIES(ISPEC))THEN
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	        IF(ATM(ID)%XzV_PRES)THEN
	          DO I=1,K
	            DO J=1,ATM(ID)%NXzV_F
	              ATM(ID)%XzV_F(J,I)=ATM(ID)%XzV_F(J,I)*T1
	            END DO
	            ATM(ID)%DXzV_F(I)=ATM(ID)%DXzV_F(I)*T1
	            IF(.NOT. ATM(ID+1)%XzV_PRES)ATM(ID+1)%XzV_F(1,I)=
	1                ATM(ID)%XzV_F(1,I)*T1
	          END DO
	        END IF
	      END DO
	    END IF
	  END DO
! 
!
	ELSE IF(XOPT .EQ. 'AV')THEN
	  CALL USR_OPTION(K,'DEPTH',' ','Depth index for averaging')
	  DO ID=1,NUM_IONS
            IF(XSPEC .EQ. UC(ION_ID(ID)) .AND. ATM(ID)%XzV_PRES)THEN
	      DO J=1,ATM(ID)%NXzV_F
	        ATM(ID)%XzV_F(J,K)=SQRT(ATM(ID)%XzV_F(J,K-1)*ATM(ID)%XzV_F(J,K+1) )
	      END DO
	      IF(ID .NE. 1)THEN
	        IF(ATM(ID-1)%XzV_PRES)THEN
                  ATM(ID-1)%DXzV_F(K)=SQRT(ATM(ID-1)%DXzV_F(K-1)*ATM(ID-1)%DXzV_F(K+1))
	        END IF
	      END IF
	      IF(ID .NE. NUM_IONS)THEN
	        IF(.NOT. ATM(ID+1)%XzV_PRES)THEN
                  ATM(ID)%DXzV_F(K)=SQRT(ATM(ID)%DXzV_F(K-1)*ATM(ID)%DXzV_F(K+1))
	        END IF
	      END IF
	      EXIT
	    END IF
	  END DO
!
	ELSE IF(XOPT .EQ. 'REP')THEN
	  CALL USR_OPTION(K,'DEPTH',' ','Depth to replace')
	  I=K+1
	  CALL USR_OPTION(I,'DEPTH',' ','Replacement depth')
	  DO ID=1,NUM_IONS
            IF(XSPEC .EQ. UC(ION_ID(ID)) .AND. ATM(ID)%XzV_PRES)THEN
	      DO J=1,ATM(ID)%NXzV_F
	        ATM(ID)%XzV_F(J,K)=ATM(ID)%XzV_F(J,I)*MASS_DENSITY(K)/MASS_DENSITY(I)
	      END DO
	      IF(ID .NE. 1)THEN
	        IF(ATM(ID-1)%XzV_PRES)THEN
                  ATM(ID-1)%DXzV_F(K)=ATM(ID-1)%DXzV_F(I)*MASS_DENSITY(K)/MASS_DENSITY(I)
	        END IF
	      END IF
	      IF(ID .NE. NUM_IONS)THEN
	        IF(.NOT. ATM(ID+1)%XzV_PRES)THEN
                  ATM(ID)%DXzV_F(K)=ATM(ID)%DXzV_F(I)*MASS_DENSITY(K)/MASS_DENSITY(I)
	        END IF
	      END IF
	      EXIT
	    END IF
	  END DO
! 
!
	ELSE IF (XOPT .EQ. 'FIXT')THEN
!
	  CALL USR_OPTION(ELEC,'DT','F','Set T values at certain depths?')
	  IF(ELEC)THEN
	    CALL USR_OPTION(LEV,ITEN,ITWO,'OWIN','0,0,0,0,0,0,0,0,0,0',
	1        'Depth section to be changed [2: lower and upper limit]')
            DEFAULT='1.0'
	    DO I=LEV(1),LEV(2)
	      CALL USR_OPTION(T1,'TEMP',DEFAULT,'Revised temperature')
	      T(I)=T1
	      DEFAULT=WR_STRING(T1)
	    END DO
	  ELSE
	    CALL USR_HIDDEN(T1,'SCALE','1.0',
	1      'Optical Depth below which T is held approximately constant')
!
! As we don't require the old T, we can overite it straight away.
!
	    DO I=1,ND
	      IF(T1 .LE. 0)THEN
	        T1=1.0D0
	      ELSE
	        T1=TAUROSS(I)/T1
	        T1=1.0D0-EXP(-T1)
	      END IF
	      T(I)=T1*TGREY(I)+(1.0-T1)*T(I)
	    END DO
!

	    DO ID=1,NUM_IONS
	      CALL UPDATE_POPS(ATM(ID)%XzV_F,ATM(ID)%XzVLTE_F,
	1                    ATM(ID)%NXzV_F,ND,T,TA)
	    END DO
	    DO I=1,ND
	      T(I)=TA(I)
	    END DO
	  END IF
!
! Compute LTE populations
!
!            INCLUDE 'EVAL_LTE_FULL.INC'
!
	  ROSS=.FALSE.
	  GREY_COMP=.FALSE.
!
! 
	ELSE IF(XOPT .EQ. 'XCOOL')THEN
!
	   WRITE(6,*)' '
	   WRITE(6,*)' Compares X-ray cooling time and flow times assuming smooth flow'
	   WRITE(6,*)' '
!
! The factor 0.375 = 1.5/4 (3/2 kT and 4 for the shock compression).
!
	   T2=BOLTZMANN_CONSTANT()
	    DO I=1,ND
	      TA(I)=0.375D0*T2*(ED(I)+POP_ATOM(I))*1.0D+28/ED(I)/POP_ATOM(I)/CLUMP_FAC(I)
	      TB(I)=1.0D+05*R(I)/V(I)
	    END DO
	    CALL DP_CURVE(ND,XV,TA)
	    CALL DP_CURVE(ND,XV,TB)
	    YAXIS='t\dX\u(T/10\u6\d)(\gL/10\u22\d) s; t(flow)'
	ELSE IF(XOPT .EQ.'XRAY')THEN
	  XRAYS=.NOT. XRAYS
	  IF(XRAYS)THEN
	    WRITE(T_OUT,*)'Xray opacities now included'
	    CALL USR_HIDDEN(FILL_FAC_XRAYS,'FILL','1.0',
	1	 'X-ray filling factor')
	    CALL USR_HIDDEN(T_SHOCK,'TSHOCK','300',
	1	 'Shock temperature of region producing X-rays')
	    CALL USR_HIDDEN(V_SHOCK,'VSHOCK','200 ',
	1	 'Minimum velocity for X-ray production')
	  ELSE
	    WRITE(T_OUT,*)'Xray opacities NOT included'
	  END IF
!
! Switch on/of level dissolution.
!
	ELSE IF(XOPT .EQ. 'DIS')THEN
	  LEVEL_DISSOLUTION= .NOT. LEVEL_DISSOLUTION
	  CALL COMP_LEV_DIS_BLK(ED,POPION,T,LEVEL_DISSOLUTION,ND)
	  INCLUDE 'EVAL_LTE_FULL.INC'
	  IF(LEVEL_DISSOLUTION)THEN
	    WRITE(T_OUT,*)'Effect of level dissolution on opacities '//
	1	 'now included'
	  ELSE
	    WRITE(T_OUT,*)'Effect of level dissolution on opacities '//
	1	 'NOT included'
	  END IF
	  WRITE(T_OUT,*)'Warning --- level dissolution option should match' ,
	1    'that of main code to avoid inconsisteancies.'
!
! 
! Computes the optical depth of lines of a given ionization stage and species at a 
! given depth in the atmosphere. Either the radial or tangential Sobolev optical
! depth is used.
!
	ELSE IF(XOPT .EQ. 'ARAT')THEN
	  WRITE(6,'(A)')RED_PEN
	  WRITE(6,'(A)')'Plotting line optical depth versus wavelength'
	  WRITE(6,'(A)')'Lines with different A/SUM ranges ar eplotted in different colors'
	  WRITE(6,'(A)')BLACK_PEN
!
	  XAXSAV=XAXIS
          CALL BRANCH_RAT(OMEGA_F,XV,YV,XAXIS,YAXIS,XSPEC,N_MAX,N_PLT_MAX,ND)
!
! Computes the optical depth of lines of a given ionization stage and species at a 
! given depth in the atmosphere. Either the radial or tangential Sobolev optical
! depth is used.
!
	ELSE IF(XOPT .EQ. 'LTAU')THEN
!
	  I=ND/2
	  DEFAULT=WR_STRING(I)
	  VALID_VALUE=.FALSE.
	  DO WHILE(.NOT. VALID_VALUE)
	    CALL USR_OPTION(I,'DEPTH',DEFAULT,'Depth for plotting TAUL')
	    IF(I .GE. 1 .AND. I .LE. ND)THEN
	       WRITE(T_OUT,'(A,I4,A,ES12.4)')'     Radius at depth',I,' is',R(I)
	       WRITE(T_OUT,'(A,I4,A,ES12.4)')'   Velocity at depth',I,' is',V(I)
	       WRITE(T_OUT,'(A,I4,A,ES12.4)')'Temperature at depth',I,' is',T(I)
	       VALID_VALUE=.TRUE.
	    END IF
	  END DO
	  CALL USR_HIDDEN(FLAG,'LOGX','F','Logarithmic in Wavelength?')
	  CALL USR_OPTION(RADIAL_TAU,'RD_TAU','TRUE',
	1      'Use radial (alt. is TANGENTIAL) direction to evaluate the Sobolev optical depth')
	  CALL USR_OPTION(LINE_STRENGTH,'LS','F','Plot line strength')
!
	  IF(LINE_STRENGTH)THEN
!
! Vdop=10 km/s and ignoring SQRT(PI)
!
	    WRITE(6,*)'Opacity at line center(for Vdop=10km/s) * SQRT(PI) normalized be e.s. opacity' 
	    TAU_CONSTANT=1.0D-15*OPLIN*2.998D+04/(6.65D-15*ED(I))
	  ELSE IF(V(I) .LE. 20.0D0)THEN
	    WRITE(6,*)'Using modified static optical depth'
	    TAU_CONSTANT=OPLIN*1.6914D-11		!Assumes Vdop=10 km/s
	  ELSE
	    TAU_CONSTANT=OPLIN*R(I)*2.998E-10/V(I)
	    IF(RADIAL_TAU)TAU_CONSTANT=TAU_CONSTANT/(1.0D0+SIGMA(I))
	  END IF
!
! Compute line opacity and emissivity.
!
	  J=0
	  FOUND=.FALSE.
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES .AND. (XSPEC .EQ. UC(ION_ID(ID)) .OR. XSPEC .EQ. 'ALL') )THEN
	      DO NL=1,ATM(ID)%NXzV_F
	        DO NUP=NL+1,ATM(ID)%NXzV_F
	          WRITE(25,*)NL,NUP,ATM(ID)%AXzV_F(NL,NUP)
	          IF(ATM(ID)%AXzV_F(NL,NUP) .NE. 0)THEN
	            J=J+1
	            T1=ATM(ID)%W_XzV_F(NUP,I)/ATM(ID)%W_XzV_F(NL,I)
	            FREQ=ATM(ID)%EDGEXzV_F(NL)-ATM(ID)%EDGEXzV_F(NUP)
	            GLDGU=ATM(ID)%GXzV_F(NL)/ATM(ID)%GXzV_F(NUP)
	            IF(FLAG)THEN
                      XV(J)=LOG10( ANG_TO_HZ/FREQ )
	            ELSE
                      XV(J)=ANG_TO_HZ/FREQ
	            END IF
	            T2=ATM(ID)%AXzV_F(NL,NUP)*(T1*ATM(ID)%XzV_F(NL,I)-GLDGU*ATM(ID)%XzV_F(NUP,I))
	            WRITE(25,*)NL,NUP,T1,T2
	            WRITE(25,*)J,XV(J),T2
	            IF(T2 .NE. 0)THEN
	              IF(V(I) .LE. 20.0D0 .AND. .NOT. LINE_STRENGTH)THEN
	                YV(J)=0.0D0
	                T3=T2
	                DO K=I,2,-1
	                  T2=T3
	                  T1=ATM(ID)%W_XzV_F(NUP,K)/ATM(ID)%W_XzV_F(NL,K)
	                  T3=T1*ATM(ID)%XzV_F(NL,K)-GLDGU*ATM(ID)%XzV_F(NUP,K)
	                  YV(J)=YV(J)+(R(K-1)-R(K))*(T2+T3)
	                  IF(V(K+1) .GT. 20.0D0)EXIT
	                END DO
	                YV(J)=LOG10( 0.5D0*ATM(ID)%AXzV_F(NL,NUP)*ABS(YV(J))*TAU_CONSTANT/FREQ)
	              ELSE
	                YV(J)=LOG10( ABS(T2)*TAU_CONSTANT/FREQ )
	              END IF
	            ELSE
	              J=J-1
	            END IF
	            WRITE(25,*)J,XV(J),T2
	            FOUND=.TRUE.
	          END IF
	        END DO
	      END DO
	    END IF
	  END DO
	  IF(.NOT. FOUND)THEN
	    WRITE(T_OUT,*)' Invalid population type or species unavailable.'
	    GOTO 1
	  END IF
!
	  XAXSAV=XAXIS
	  IF(FLAG)THEN
	    XAXIS='Log(\gl(\A))'
	  ELSE
	    XAXIS='\gl(\gV)'
	  END IF
	  YAXIS='Log(\gt)'
!
	  IF(J .NE. 0)WRITE(T_OUT,*)J,' lines plotted'
	  IF(J  .NE. 0)CALL DP_CURVE(J,XV,YV)
!
!`
	ELSE IF(XOPT .EQ. 'VLTAU')THEN
!
	  DEFAULT=WR_STRING(LAM_ST)
	  CALL USR_OPTION(LAM_ST,'LAMST',DEFAULT,FREQ_INPUT)
	  DEFAULT=WR_STRING(LAM_EN)
	  CALL USR_OPTION(LAM_EN,'LAMEN',DEFAULT,FREQ_INPUT)
	  NU_ST=ANG_TO_HZ/LAM_ST
	  NU_EN=ANG_TO_HZ/LAM_EN
	  T1=NU_ST
	  NU_ST=MIN(NU_ST,NU_EN)
	  NU_EN=MAX(T1,NU_EN)
!
	  DO I=1,ND
	    TA(I)=OPLIN*R(I)*2.998D-10/V(I)
	  END DO
!
	  J=0
	  FOUND=.FALSE.
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES .AND. (XSPEC .EQ. UC(ION_ID(ID)) .OR. XSPEC .EQ. 'ALL') )THEN
	      DO NL=1,ATM(ID)%NXzV_F
	        DO NUP=NL+1,ATM(ID)%NXzV_F
	          FREQ=ATM(ID)%EDGEXzV_F(NL)-ATM(ID)%EDGEXzV_F(NUP)
	          GLDGU=ATM(ID)%GXzV_F(NL)/ATM(ID)%GXzV_F(NUP)
	          IF(ATM(ID)%AXzV_F(NL,NUP) .NE. 0 .AND. FREQ .GT. NU_ST .AND. FREQ .LT. NU_EN)THEN
	            J=J+1
	            IF(FLAG)THEN
                      XV(J)=LOG10( ANG_TO_HZ/FREQ )
	            ELSE
                      XV(J)=ANG_TO_HZ/FREQ
	            END IF
	            ELEC=.FALSE.
	            DO I=1,ND
	              T1=ATM(ID)%W_XzV_F(NUP,I)/ATM(ID)%W_XzV_F(NL,I)
	              T2=ATM(ID)%AXzV_F(NL,NUP)*(T1*ATM(ID)%XzV_F(NL,I)-GLDGU*ATM(ID)%XzV_F(NUP,I))
	              T2=T2*TA(I)/FREQ
	              IF(ELEC .AND. T2 .LT. 1.0D0)THEN
	                ZV(J)=(YV(J)-V(I))/2
	                YV(J)=(YV(J)+V(I))/2
	                EXIT
	              ELSE IF(T2 .GT. 1.0D0 .AND. .NOT. ELEC)THEN
	                ELEC=.TRUE.
	                YV(J)=V(I)
	              ELSE IF(ELEC .AND. I .EQ. ND)THEN
	                ZV(J)=(YV(J)-V(I))/2
	                YV(J)=(YV(J)+V(I))/2
	              ELSE IF(I .EQ. ND)THEN
	                J=J-1
	              END IF
	            END DO
	          END IF
	        END DO
	      END DO
	    END IF
	  END DO
!
	  XAXSAV=XAXIS
	  IF(FLAG)THEN
	    XAXIS='Log(\gl(\A))'
	  ELSE
	    XAXIS='\gl(\gV)'
	  END IF
	  YAXIS='V(km/s)'
	  IF(J .NE. 0)WRITE(T_OUT,*)J,' lines plotted'
	  IF(J  .NE. 0)CALL DP_CURVE_AND_ER(J,XV,YV,ZV,' ')
!
!
! This section creates an image of the partial opacities (i.e., the opacities due
! to each ionizaion stage) at a given frequency as a function of depth.
! 
	ELSE IF(XOPT .EQ. 'DCHI')THEN
	  CALL USR_OPTION(FREQ,'LAM','0.0',FREQ_INPUT)
	  IF(FREQ .EQ. 0)GOTO 1
	  IF(KEV_INPUT)THEN
	    FREQ=FREQ*KEV_TO_HZ
	  ELSE IF(ANG_INPUT)THEN
	    FREQ=ANG_TO_HZ/FREQ
	  END IF
	  FL=FREQ
!
	  INCLUDE 'PAR_OPACITIES.INC'
	  CALL USR_OPTION(ELEC,'ELEC','F','Include electron scattering?')
	  IF(ELEC)THEN
	    CALL ESOPAC(CHI,ED,ND)
	  ELSE
	    CHI(1:ND)=0.0D0
	  END IF 
!
	  DO ID=1,NUM_IONS
	    CHI(1:ND)=CHI(1:ND)+CHI_PAR(1:ND,ID)
	  END DO
	  WRITE(6,*)'Computed the total opacity'
!
	  DO ID=1,NUM_IONS
	    CHI_PAR(:,ID)=CHI_PAR(:,ID)/CHI(1:ND)
            YMAPV(ID)=ID
	  END DO
	  CALL USR_OPTION(ELEC,'Depth','T','Plot against depth index')
	  IF(ELEC)THEN
	    DO I=1,ND
	      XMAPV(I)=I
	    END DO
	  ELSE
	    DO I=1,ND
	      WRITE(6,*)'I=',I,ED(I)
	      XMAPV(I)=LOG10(ED(I))
	    END DO
	  END IF
!
! Remove rows containing all zero's.
!
	  ETA_PAR(:,:)=CHI_PAR(:,:)
	  J=0
	  LAST_NON_ZERO=.FALSE.
	  DO ID=1,NUM_IONS
	    IF(SUM(ETA_PAR(:,ID)) .NE. 0)THEN
	      J=J+1
              CHI_PAR(:,J)=ETA_PAR(:,ID)
	      LAST_NON_ZERO=.TRUE.
	    ELSE IF(LAST_NON_ZERO)THEN
	      J=J+1
              CHI_PAR(:,J)=ETA_PAR(:,ID)
	      LAST_NON_ZERO=.FALSE.
	      WRITE(6,*)J,ION_ID(ID)
	    ELSE
	    END IF
	  END DO
!
! ETA_PAR is passed, but should not be imaged.
!
	  CALL MAP_PLOT(CHI_PAR,ETA_PAR,XMAPV,YMAPV,ND,J,
	1               " "," "," "," ")
!
!
!
! This section does one of two things:
!           1. Creates an image of the fractional contribution to the opacity
!              and emissivity of a single species as a function of depth and
!              wavelength.
!           2. When no species is specified, the maximum fractional contribution
!               of each species at ANY depth is ouput to the terminal
!           3. Plot species dependent contributions at a given depth.
!           
	ELSE IF(XOPT .EQ. 'MCHI')THEN
	  CALL USR_OPTION(LAM_ST,'LAM_ST','50.0',FREQ_INPUT)
	  CALL USR_OPTION(LAM_EN,'LAM_END','10000.0',FREQ_INPUT)
	  CALL USR_OPTION(NLAM,'NLAM','1000.0','# of wavlengths')
	  CALL USR_OPTION(ELEC,'ELEC','F','Include electron scattering?')
	  CALL USR_OPTION(DO_KAP,'DO_KAP','F','Plot species dependent contributions?')
	  IF(DO_KAP)THEN
	    CALL USR_OPTION(DPTH_INDX,'DPTH','30','Depth index?')
	  END IF
!
	  LOC_ION_ID='$$'
	  DO WHILE(LOC_ION_ID .EQ. '$$' .AND. .NOT. DO_KAP)
	    CALL USR_OPTION(LOC_ION_ID,'ION',' ','Ion ID (eg HI)')
	    DO ID=1,NUM_IONS
	      IF(LOC_ION_ID .EQ. ION_ID(ID) .OR. LOC_ION_ID .EQ. ' ')THEN
	         EXIT
	      END IF
	      IF(ID .EQ. NUM_IONS)THEN
	        WRITE(6,*)'Invalid ION_ID',LOC_ION_ID
	        LOC_ION_ID='$$'
	      END IF
	    END DO
	  END DO
!
	  IF(KEV_INPUT)THEN
	    NU_ST=LAM_ST*KEV_TO_HZ
	    NU_EN=LAM_EN*KEV_TO_HZ
	  ELSE IF(ANG_INPUT)THEN
	    NU_EN=ANG_TO_HZ/LAM_ST
	    NU_ST=ANG_TO_HZ/LAM_EN
	  ELSE
	    NU_ST=LAM_ST
	    NU_EN=LAM_EN
	  END IF 
	  DEL_NU=10.0D0**( LOG10(NU_EN/NU_ST)/(NLAM-1) )
!
	  IF(LOC_ION_ID .NE. ' ')THEN
	    ALLOCATE (CHI_LAM(NLAM,ND))
	    ALLOCATE (ETA_LAM(NLAM,ND))
	    ALLOCATE (LAM_VEC(NLAM));     LAM_VEC(:)=0.0D0
	  ELSE IF(DO_KAP)THEN
	    ALLOCATE (CHI_LAM(NLAM,NUM_IONS))
	    ALLOCATE (ETA_LAM(NLAM,NUM_IONS))
	    ALLOCATE (LAM_VEC(NLAM));     LAM_VEC(:)=0.0D0
	  ELSE
	    ALLOCATE (CHI_TOT_LAM(NUM_IONS)); CHI_TOT_LAM(:)=0.0D0
	    ALLOCATE (LAM_VEC(NUM_IONS));     LAM_VEC(:)=0.0D0
	  END IF
!
	  LAST_NON_ZERO=.FALSE.
	  NION=0
	  DO ML=1,NLAM
	    WRITE(6,*)'ML=',ML
	    FREQ=NU_ST*(DEL_NU**(ML-1))
	    FL=FREQ
!
	    INCLUDE 'PAR_OPACITIES.INC'
	    IF(ELEC)THEN
	      CALL ESOPAC(CHI,ED,ND)
	    ELSE
	      CHI(1:ND)=0.0D0
	    END IF
	    ETA(1:ND)=0.0D0
!
	    DO ID=1,NUM_IONS
	      CHI(1:ND)=CHI(1:ND)+CHI_PAR(1:ND,ID)
	      ETA(1:ND)=ETA(1:ND)+ETA_PAR(1:ND,ID)
	    END DO
!
! NB: First dimension of CHI_PAR is ND.
!
	    IF(.NOT. DO_KAP)THEN
	      DO ID=1,NUM_IONS
	        CHI_PAR(:,ID)=CHI_PAR(:,ID)/CHI(1:ND)
	        ETA_PAR(:,ID)=ETA_PAR(:,ID)/ETA(1:ND)
                YMAPV(ID)=ID
	      END DO
	    END IF
!
	    IF(LOC_ION_ID .EQ. ' ')THEN
	      DO ID=1,NUM_IONS
	        T1=MAXVAL(CHI_PAR(:,ID))
	        T2=MAXVAL(ETA_PAR(:,ID))
	        CHI_TOT_LAM(ID)=MAX(T1,CHI_TOT_LAM(ID))
	        LAM_VEC(ID)=MAX(T2,LAM_VEC(ID))
	      END DO
	    ELSE IF(DO_KAP)THEN
	       LAM_VEC(ML)=LOG10(ANG_TO_HZ/FREQ)
	       DO ID=1,NUM_IONS
	         CHI_LAM(ML,ID)=CHI_PAR(DPTH_INDX,ID)
	         ETA_LAM(ML,ID)=ETA_PAR(DPTH_INDX,ID)
	       END DO
	    ELSE
	       LAM_VEC(ML)=LOG10(ANG_TO_HZ/FREQ)
	       DO ID=1,NUM_IONS
	         IF(LOC_ION_ID .EQ. ION_ID(ID))THEN
	           CHI_LAM(ML,:)=CHI_PAR(:,ID)
	           ETA_LAM(ML,:)=ETA_PAR(:,ID)
	         END IF
	      END DO
	    END IF
	  END DO
!
	  IF(LOC_ION_ID .EQ. ' ')THEN
	    WRITE(6,'(2X,A3,2X,A7,A10,A10)')'ID','  ION  ','  % CHI   ','  % ETA'
	    DO ID=1,NUM_IONS
	      IF(CHI_TOT_LAM(ID) .NE. 0)THEN
	         ISPEC=SPECIES_LNK(ID)
	         WRITE(6,'(2X,I3,2X,A7,2F10.3)')
	1          ID,TRIM(ION_ID(ID)),100.0D0*CHI_TOT_LAM(ID),100.0D0*LAM_VEC(ID)
	      END IF
	    END DO
	    DEALLOCATE (CHI_TOT_LAM)
	    DEALLOCATE (LAM_VEC)
	  ELSE IF(DO_KAP)THEN
	    DO ID=1,NUM_IONS
	      CHI_LAM(:,ID)=1.0D-10*CHI_LAM(:,ID)/MASS_DENSITY(DPTH_INDX)
	    END DO
	    DO ISPEC=1,NSPEC
	      YV(1:NLAM)=0.0D0
	      ELEC=.FALSE.
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	        YV(1:NLAM)=YV(1:NLAM)+CHI_LAM(1:NLAM,ID)
	        ELEC=.TRUE.
	      END DO
	      IF(ELEC)THEN
	        CALL DP_CURVE(NLAM,LAM_VEC,YV)
	        WRITE(6,*)ISPEC,SPECIES(ISPEC)
	      END IF
	    END DO
	  ELSE
	    DO I=1,ND
              XMAPV(I)=I
	    END DO
	    CALL MAP_PLOT(CHI_LAM,ETA_LAM,LAM_VEC,XMAPV,NLAM,ND,
	1               " "," "," "," ")
!
	    DEALLOCATE (CHI_LAM)
	    DEALLOCATE (ETA_LAM)
	    DEALLOCATE (LAM_VEC)
	  END IF
!
!
! This section creates an image of the fractional contribution to the opacity
!              as a function of species and wavelength.
!
	ELSE IF(XOPT .EQ. 'LCHI' .OR. XOPT .EQ. 'PPHOT')THEN
	  CALL USR_OPTION(LAM_ST,'LAM_ST','50.0',FREQ_INPUT)
	  CALL USR_OPTION(LAM_EN,'LAM_END','10000.0',FREQ_INPUT)
	  CALL USR_OPTION(K,'DEPTH','20.0','Depth for opacities')
	  CALL USR_OPTION(NLAM,'NLAM','1000.0','# of wavlengths')
	  CALL USR_OPTION(ELEC,'ELEC','F','Include electron scattering?')
!
	  IF(KEV_INPUT)THEN
	    NU_ST=LAM_ST*KEV_TO_HZ
	    NU_EN=LAM_EN*KEV_TO_HZ
	  ELSE IF(ANG_INPUT)THEN
	    NU_EN=ANG_TO_HZ/LAM_ST
	    NU_ST=ANG_TO_HZ/LAM_EN
	  ELSE
	   NU_ST=LAM_ST
	   NU_EN=LAM_EN
	  END IF
!
	  LAST_NON_ZERO=.FALSE.
	  NION=0
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES)THEN
	      NION=NION+1
	      LAST_NON_ZERO=.TRUE.
	    ELSE IF(LAST_NON_ZERO)THEN
	      NION=NION+1
	      LAST_NON_ZERO=.FALSE.
	    ELSE
	    END IF
	  END DO
!
	  ALLOCATE (CHI_LAM(NLAM,NION))
	  ALLOCATE (ETA_LAM(NLAM,NION))
	  ALLOCATE (CHI_TOT_LAM(NLAM))
	  ALLOCATE (ETA_TOT_LAM(NLAM))
	  ALLOCATE (LAM_VEC(NLAM))
	  DEL_NU=EXP(LOG(NU_EN/NU_ST)/(NLAM-1))
!
	  ETA_LAM(:,:)=0.0D0
	  CHI_LAM(:,:)=0.0D0
	  CHI_TOT_LAM(:)=0.0D0
	  ETA_TOT_LAM(:)=0.0D0
	  WRITE(6,'(A3)',ADVANCE='NO')'ML='
	  DO ML=1,NLAM
	    IF(MOD(ML,200) .EQ. 0)THEN
	        WRITE(6,*)' '
	        WRITE(6,'(A3)',ADVANCE='NO')'ML='
	    END IF
	    IF(ML .EQ. NLAM)WRITE(6,'(A)')' '
	    IF(MOD(ML,20) .EQ. 0)WRITE(6,'(1X,I5)',ADVANCE='NO')ML
	    FREQ=NU_ST*(DEL_NU**(ML-1))
	    LAM_VEC(ML)=LOG10(ANG_TO_HZ/FREQ)
	    INCLUDE 'PAR_OPACITIES.INC'
	    J=0; CNT=0
	    DO ID=1,NUM_IONS
	      IF(ATM(ID)%XzV_PRES)THEN
	        J=J+1
                CHI_LAM(ML,J)=CHI_PAR(K,ID)
                ETA_LAM(ML,J)=ETA_PAR(K,ID)
	        LAST_NON_ZERO=.TRUE.
	        CHI_TOT_LAM(ML)=CHI_TOT_LAM(ML)+CHI_LAM(ML,J)
	        ETA_TOT_LAM(ML)=ETA_TOT_LAM(ML)+ETA_LAM(ML,J)
	        IF(XOPT .NE. 'PPHOT' .AND. ML .EQ. NLAM)THEN
	          CNT=CNT+1
	          IF(MOD(CNT-1,5) .EQ. 0)THEN
	            WRITE(6,'(4X,A5,1X,I3)')TRIM(ION_ID(ID)),J
	          ELSE
	            WRITE(6,'(4X,A5,1X,I3)',ADVANCE='NO')TRIM(ION_ID(ID)),J
	          END IF
	        END IF
	      ELSE IF(LAST_NON_ZERO)THEN
	        J=J+1
                CHI_LAM(ML,J)=0.0D0
                ETA_LAM(ML,J)=0.0D0
	        LAST_NON_ZERO=.FALSE.
	      ELSE
	      END IF
	    END DO
	  END DO
	  IF(XOPT .NE. 'PPHOT')WRITE(6,'(A,/)')' '
	  NION=J
!
	  CALL ESOPAC(ESEC,ED,ND)
	  IF(ELEC)THEN
	    CHI_TOT_LAM(1:NLAM)=CHI_TOT_LAM(1:NLAM)+ESEC(K)
	  END IF 
!
	  IF(XOPT .EQ. 'PPHOT')THEN
	    WRITE(6,*)'Total opacity is plotted in red'
	    XV(1:NLAM)=10**(LAM_VEC(1:NLAM))
	    YV(1:NLAM)=CHI_TOT_LAM(1:NLAM)/ESEC(K)
	    CALL DP_CURVE(NLAM,XV,YV)
	    YAXIS='Opacity/(e.s. opacity)'
	    DO ID=1,NUM_IONS
	      IF(ATM(ID)%XzV_PRES .AND. XSPEC .EQ. UC(ION_ID(ID)) )THEN
	        YV(1:NLAM)=CHI_LAM(1:NLAM,ID)
	        CALL DP_CURVE(NLAM,XV,YV)
	        EXIT
	      END IF
	    END DO
	    DO ISPEC=1,NSPEC
	      IF(XSPEC .EQ. UC(SPECIES(ISPEC)) )THEN
	        DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	          YV(1:NLAM)=CHI_LAM(1:NLAM,ID)
	          CALL DP_CURVE(NLAM,XV,YV)
	        END DO
	      END IF
	    END DO
	  ELSE
	    DO ID=1,NION
	      DO ML=1,NLAM
	        CHI_LAM(ML,ID)=CHI_LAM(ML,ID)/CHI_TOT_LAM(ML)
	        ETA_LAM(ML,ID)=ETA_LAM(ML,ID)/ETA_TOT_LAM(ML)
	      END DO
	    END DO
!
	    DO ID=1,NION
              YMAPV(ID)=ID
	    END DO
!
	    CALL MAP_PLOT(CHI_LAM,ETA_LAM,LAM_VEC,YMAPV,NLAM,NION,
	1               " "," "," "," ")
	  END IF
!
	  DEALLOCATE (CHI_LAM)
	  DEALLOCATE (ETA_LAM)
	  DEALLOCATE (LAM_VEC)
	  DEALLOCATE (CHI_TOT_LAM)
	  DEALLOCATE (ETA_TOT_LAM)
! 
!
! Graphical options.
!
	ELSE IF(XOPT .EQ. 'GR') THEN
	  CALL GRAMON_PGPLOT(XAXIS,YAXIS,NAME,' ')
	  XAXIS=XAXSAV
!
	ELSE IF(XOPT .EQ. 'GRNL') THEN
	  CALL GRAMON_PGPLOT(' ',' ',' ',' ')
	  XAXIS=XAXSAV
!
	ELSE IF(XOPT .EQ. 'GRED')THEN
	  TA(1:ND)=XV(1:ND)
!
	  I=LOG10(ED(1));J=LOG10(ED(ND))
	  IF(I .LT. ED(1))I=I+1
	  NXED=(J-I)*2+1
	  NXED=MIN(NXED,NSC)
	  SCED(1)=I
	  DO I=2,NXED
	    SCED(I)=SCED(1)+(I-1)*0.5D0
	  END DO
	  DO I=1,ND
	    TB(I)=DLOG10(ED(I))
	  END DO
	  CALL MON_INTERP(XED,NXED,IONE,SCED,NXED,TA,ND,TB,ND)
	  TOPLABEL='Log(N\de\u)'
	  CALL GRAMON_PGPLOT(XAXIS,YAXIS,NAME,'TOPLAB')
	  XAXIS=XAXSAV
!
	ELSE IF(XOPT .EQ. 'GROSS')THEN
	  DO I=1,ND
	    CHIROSS(I)=CLUMP_FAC(I)*ROSS_MEAN(I)
	  END DO
	  CALL TORSCL_V2(TAUROSS,CHIROSS,R,TB,TC,ND,METHOD,TYPE_ATM,L_FALSE)
	  WRITE(T_OUT,*)'Rossland optical depth is : ',TAUROSS(ND)
!
	  TB(1:ND)=LOG10(TAUROSS(1:ND))
	  TA(1:ND)=XV(1:ND)
!
	  I=TB(1); J=TB(ND)
	  IF(I .LT. TB(1))I=I+1
	  IF(J .GT. TB(ND))J=J-1
	  SCED(1)=I
	  NXED= (J-I)*2+1
	  NXED=MIN(NXED,NSC)
	  DO I=2,NXED
	    SCED(I)=SCED(I-1)+0.5D0
	  END DO
	  CALL MON_INTERP(XED,NXED,IONE,SCED,NXED,TA,ND,TB,ND)
	  TOPLABEL='Log(\gt\dRoss\u)'
	  CALL GRAMON_PGPLOT(XAXIS,YAXIS,NAME,'TOPLAB')
	  XAXIS=XAXSAV
!
! 
	ELSE IF(XOPT .EQ. 'SM')THEN
	  FLAG=.FALSE.
	  CALL USR_OPTION(T1,'VAL','3.0','Fractional change allowed.')
!
	  DEFAULT=WR_STRING(ND/2)
	  CALL USR_OPTION(K,'DEPTH',DEFAULT,'Smooth from DEPTH to 1')
	  CALL USR_HIDDEN(TYPE,'TYPE','POP','Smooth in POP (def) or DC')
!
	  FOUND=.FALSE.
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES)THEN
	      CALL SMOOTH_V2(ATM(ID)%XzV_F,ATM(ID)%XzVLTE_F,T1,
	1         ATM(ID)%NXzV_F,TYPE,K,
	1         ND,XSPEC,ION_ID(ID),FLAG)
	      CALL SMOOTH_V2(ATM(ID)%DXzV_F,POPION,T1,1,TYPE,K,
	1         ND,XSPEC,ION_ID(ID),FLAG)
	      FOUND=.TRUE.
	    END IF
	  END DO
!
	  IF(.NOT. FOUND)THEN
	    WRITE(T_OUT,*)' Invalid population type or species unavailable.'
	    GOTO 1
	  END IF
	  IF(.NOT. FLAG)THEN
	    WRITE(T_OUT,*)'Warning - invalid SM extension.'
	  END IF
!
! 
!
! Allow the OCCUPATION PROBAILITIES to be plotted.
!
	ELSE IF(XOPT .EQ. 'WLD')THEN
	  CALL USR_OPTION(LEV,ITEN,IONE,'LEVS','1,2,3,4,5,6,7,8,9,10',
	1      'Levels to plot (10 max)')
	  J=1
	  FOUND=.FALSE.
	  DO WHILE(LEV(J) .GT. 0 .AND. LEV(J) .LE. 200)
	    DO ID=1,NUM_IONS
	      IF(XSPEC .EQ. ION_ID(ID))THEN
	        DO I=1,ND
	          YV(I)=ATM(ID)%W_XzV_F(LEV(J),I)
	        END DO
	        FOUND=.TRUE.
	      END IF
	    END DO
	    IF(.NOT. FOUND)THEN
	      WRITE(T_OUT,*)'Invalid species descritor ', X
	    END IF
	    J=J+1
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='WLD'
	  END DO
!
! 
!
! ***********************************************************************
! ***********************************************************************
!
! DC  --- Plot of Log(departure coefficient).
!           IF(LIN_DC is TRUE, a linear plot is preseneted)
!
! POP --- LOG plot of the population of an individual level.
!
! RAT --- Plot of the ionization fraction of a species. Fraction can be
!           relative to the ion density, or the total species density.
!
! ***********************************************************************
! ***********************************************************************
!
	ELSE IF(XOPT .EQ.'DC' .OR. XOPT .EQ. 'POP'
	1             .OR. XOPT .EQ. 'RAT' .OR. XOPT .EQ. 'TX') THEN
	  DO I=1,10
	    LEV(I)=0
	  END DO
	  IF(XOPT(1:2) .EQ. 'DC')CALL USR_HIDDEN(LIN_DC,'LIN','F','Linear D.C. plots?')
	  CALL USR_HIDDEN(SPEC_FRAC,'SPEC_FRAC','F','Species fraction?')
	  IF(XOPT .EQ. 'RAT')THEN
	   HAM=.TRUE.
	  ELSE
	   HAM=.FALSE.
	  END IF
	  DEFAULT=WR_STRING(HAM)
	  CALL USR_HIDDEN(HAM,'SUM',DEFAULT,'Sum all ion levels?')
	  IF(HAM)THEN
	    LEV(1)=1000
	    LEV(2)=0
	  ELSE
	    CALL USR_OPTION(LEV,ITEN,IONE,'LEVS','1,2,3,4,5,6,7,8,9,10',
	1	 'Levels to plot (10 max)')
	  END IF
!
	  FOUND=.FALSE.
	  DO ID=1,NUM_IONS
	    IF(XSPEC .EQ. UC(ION_ID(ID)) .AND. ATM(ID)%XzV_PRES)THEN
	      DO I=1,10
	        IF(LEV(I) .EQ. 0 .OR. LEV(I) .GT. 20000)GOTO 1
	          FLAG=.FALSE.
	          CALL SET_DC_OR_POP_OR_TX_V2(YV,LEV(I),ATM(ID)%XzV_F,ATM(ID)%LOG_XzVLTE_F,
	1            ATM(ID)%EDGEXzV_F,ATM(ID)%NXzV_F,ND,T,X,UC(ION_ID(ID)),FLAG)
	          IF(.NOT. ATM(ID+1)%XzV_PRES)THEN
	            T1=0.0D0
	            CALL SET_DC_OR_POP_OR_TX(YV,LEV(I),ATM(ID)%DXzV_F,ATM(ID)%DXzV_F,
	1                   T1,IONE,ND,T,X,'D'//TRIM(UC(ION_ID(ID))),FLAG)
	          END IF
!
	        IF(XOPT .EQ. 'DC')THEN
	          YAXIS='Log(b)'            
	          IF(LIN_DC)THEN		!Saves altering SETDC_OR_POP
	            DO J=1,ND
	              YV(J)=10.0**(YV(J))
	             END DO
	             YAXIS='Log(b)'
	          END IF                 
	        ELSE IF(XOPT .EQ. 'TX')THEN
	          YAXIS='T\dX\u(10\u4\ \dK)'            
	        ELSE IF(XOPT .EQ. 'RAT')THEN
	          IF(SPEC_FRAC)THEN
	          ELSE
	            DO J=1,ND
	              YV(J)=YV(J)-DLOG10(POP_ATOM(J))
	            END DO
	            YAXIS='Log(N\dx\u/N\di\u)'
	          END IF
	        ELSE
	          YAXIS='Log(N)'
	          IF(SPEC_FRAC)THEN
	            ISPEC=SPECIES_LNK(ID)
	            DO J=1,ND
	              YV(J)=YV(J)-LOG10(POPDUM(J,ISPEC))
	            END DO
	            YAXIS='Log(N/Nsp)'
	          END IF
	        END IF
!
	        CALL DP_CURVE(ND,XV,YV)
	      END DO
	    END IF
	  END DO
!
! Plots departure coeeficient for each level at each depth.
!
	ELSE IF(XOPT .EQ. 'DCS')THEN
!
	  WRITE(6,*)' '
	  WRITE(6,*)BLUE_PEN//'Use B option in PLT_SPEC for plotting'//DEF_PEN
	  WRITE(6,*)' '
	  FOUND=.FALSE.
	  DO ID=1,NUM_IONS
	    IF( (XSPEC .EQ. UC(ION_ID(ID)) .OR. XSPEC .EQ. 'ALL') .AND. ATM(ID)%XzV_PRES)THEN
	      K=0
	      DO L=1,ND
	        DO J=1,ATM(ID)%NXzV_F
	          K=K+1
	          ZV(K)=XV(L)
	          YV(K)=LOG10(ATM(ID)%XzV_F(J,L)/ATM(ID)%XzVLTE_F(J,L))
	        END DO
	      END DO
	      CALL DP_CURVE(K,ZV,YV)
	      FOUND=.TRUE.
	      YAXIS='Log b'
	    END IF
	  END DO
!
	ELSE IF(XOPT .EQ.'TPOP' .OR. XOPT .EQ. 'LPOP' .OR. XOPT .EQ. 'LDC')THEN
!
	  CALL USR_OPTION(K,'DEPTH','20','Depth at which quantities are to be plotted.')
	  FOUND=.FALSE.
	  XAXIS='Level'
	  DO ID=1,NUM_IONS
	    IF(XSPEC .EQ. UC(ION_ID(ID)) .AND. ATM(ID)%XzV_PRES)THEN
	      L=ATM(ID)%NXzV_F
	      IF(XOPT .EQ. 'LDC')THEN
	        DO I=1,ATM(ID)%NXzV_F
	          YV(I)=LOG10( ATM(ID)%XzV_F(I,K) / ATM(ID)%XzVLTE_F(I,K) )
	          ZV(I)=I
	        END DO
	        YAXIS='log b'
	      ELSE 
	        DO I=1,ATM(ID)%NXzV_F
	          YV(I)=LOG10(ATM(ID)%XzV_F(I,K))
	          ZV(I)=I
	        END DO
	        YAXIS='log N'
	      END IF
	      FOUND=.TRUE.
	    END IF
	  END DO
	  IF(XOPT .EQ. 'TPOP' .AND. FOUND)THEN
	    WRITE(6,*)'Wavelength set to 5000 Angstroms'
	    FL=ANG_TO_HZ/5000.0D0
	    T1=LOG10( 3.0D-10*OPLIN*R(K)/V(K)/(1.0D0+SIGMA(K))/FL )
	    DO I=1,L
	      YV(I)=YV(I)+T1
	    END DO
	    YAXIS='Tau/f'
	  END IF
	  IF(FOUND)THEN
	    CALL DP_CURVE(L,ZV,YV)
	  ELSE
	    WRITE(6,*)'Species not found'
	  END IF
!
! 
!
! Compute the ionization ratio for the total population of one ionization
! state to the total population of the next ionization state.
!
	ELSE IF(XOPT .EQ.'ION') THEN
!
	  FOUND=.FALSE.
	  DO ID=1,NUM_IONS
	    IF(XSPEC .EQ. ION_ID(ID))THEN
	      IF(ATM(ID+1)%XzV_PRES)THEN
	        CALL COMPION(ATM(ID)%XzV_F,ATM(ID)%NXzV_F,
	1                ATM(ID+1)%XzV_F,ATM(ID+1)%NXzV_F,YV,ND)
	      ELSE
	        CALL COMPION(ATM(ID)%XzV_F,ATM(ID)%NXzV_F,ATM(ID)%DXzV_F,IONE,YV,ND)
	      END IF
	      FOUND=.TRUE.
	    END IF
	  END DO
	  IF(.NOT. FOUND)THEN
	    WRITE(T_OUT,*)' Invalid species'
	    GOTO 1
	  END IF
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='Log(X\u(n-1)+\d/X\un+\d)'
! 
! ^L
!
! Computes the column density for 2 successive ionization stages, and for
! the atom. If XV is the column density for species XzV, the program returns
! XV/X and XSIX/X where X is the species column density. Inserted for comparison
! with SETI results.
!
	ELSE IF(XOPT .EQ. 'QF') THEN
	  DO ISPEC=1,NSPEC
	    IF(XSPEC .EQ. SPECIES(ISPEC))THEN
	      ZV(ND)=0.0D0
	      YAXIS='q'
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	        IF(ATM(ID)%XzV_PRES)THEN
	          DO I=1,ND
	            T2=0.0D0
	            DO J=1,ATM(ID)%NXzV_F
	              T2=T2+ATM(ID)%XzV_F(J,I)
	            END DO
	            TA(I)=POPDUM(I,ISPEC)
	            TB(I)=T2
	            TC(I)=ATM(ID)%DXzV_F(I)
	          END DO
	          CALL TORSCL(DTAU,TA,R,Z,XM,ND,METHOD,TYPE_ATM)
	          TA(1:ND)=DTAU(1:ND)
	          CALL TORSCL(DTAU,TB,R,Z,XM,ND,METHOD,TYPE_ATM)
	          TB(1:ND)=DTAU(1:ND)
	          CALL TORSCL(DTAU,TC,R,Z,XM,ND,METHOD,TYPE_ATM)
	          TC(1:ND)=DTAU(1:ND)
!
	          DO I=1,ND
	            YV(I)=LOG10(TB(I)/TA(I))
	            ZV(I)=LOG10(TC(I)/TA(I))
	          END DO
	          CALL DP_CURVE(ND,XV,YV)
	        END IF
	      END DO
	      IF(ZV(ND) .NE. 0.0D0)CALL DP_CURVE(ND,XV,ZV)
	      YAXIS='CD(XV)/CD(X); CD(XSIX)/CD(X)'
	    END IF
	  END DO
! 
!
! Compute the ionization ratio for the total population of one ionization
! state to the total population of the next ionization state.
!
	ELSE IF(XOPT .EQ. 'IF') THEN
	  CALL USR_OPTION(SPEC_FRAC,'SPEC_FRAC','T','Species fraction?')
	  DO ISPEC=1,NSPEC
	    IF(XSPEC .EQ. SPECIES(ISPEC))THEN
	      IF(SPEC_FRAC)THEN
	        DO I=1,ND
	          TA(I)=POPDUM(I,ISPEC)
	        END DO
	        YAXIS='Log '//TRIM(SPECIES_ABR(ISPEC))//
	1                       '\un+\d/N('//TRIM(SPECIES_ABR(ISPEC))//')'
	      ELSE
	        DO I=1,ND
	          TA(I)=POP_ATOM(I)
	        END DO
	        YAXIS='Log '//TRIM(SPECIES_ABR(ISPEC))//'\un+\d/N(total)'
	      END IF
	      CNT=0
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	        IF(ATM(ID)%XzV_PRES)THEN
	          DO I=1,ND
	            T1=0.0D0			!Using T1 avoids underflow
	            DO J=1,ATM(ID)%NXzV_F
	              T1=T1+ATM(ID)%XzV_F(J,I)
	            END DO
	            YV(I)=LOG10(T1/TA(I))
	            ZV(I)=LOG10(ATM(ID)%DXzV_F(I)/TA(I))
	          END DO
	          CALL DP_CURVE(ND,XV,YV)
	          CNT=CNT+1
	          WRITE(6,'(A,I2,A,A)')' Curve ',CNT,' is due to: ',TRIM(ION_ID(ID))
	        END IF
	      END DO
	      CALL DP_CURVE(ND,XV,ZV)
	      CNT=CNT+1
	      ID=SPECIES_END_ID(ISPEC)
	      WRITE(6,'(A,I2,A,A)')' Curve ',CNT,' is due to: ',TRIM(ION_ID(ID))
	    END IF
	  END DO
! 
	ELSE IF(XOPT .EQ. 'MODSUM')THEN
!
! This option was installed to revise the incorrect vaules of R(Tau) and
! Tau output by CMFGEN to MOD_SUM when clumping is present. Problem potentially
! applies to all MOD_SUM files with a CMFGEN program date earlier than
! 01-Sep-199.
!
! Compute the Radius and Velocity at Tau=10, and at Tau=2/3.
!
	  TCHI(1:ND)=ROSS_MEAN(1:ND)*CLUMP_FAC(1:ND)
	  CALL DERIVCHI(dCHIdR,TCHI,R,ND,METHOD)
          CALL NORDTAU(DTAU,TCHI,R,R,dCHIdR,ND)
	  TA(1:ND)=0.0D0 ; DO I=2,ND ; TA(I) = TA(I-1)+DTAU(I-1) ; END DO
	  TB(1)=2.0D0/3.0D0  ; TB(2)=10.0D0
	  CALL MON_INTERP(TC,ITWO,IONE,TB,ITWO,R,ND,TA,ND)
	  CALL MON_INTERP(AV,ITWO,IONE,TB,ITWO,V,ND,TA,ND)
!
	  OPEN(UNIT=LU_IN,FILE='MOD_SUM',STATUS='OLD',ACTION='READ')
	  OPEN(UNIT=LU_OUT,FILE='REV_MOD_SUM',STATUS='UNKNOWN')
	   READ(LU_IN,'(A)')STRING
	   DO WHILE( INDEX(STRING,'Tau=') .EQ. 0)
	     WRITE(LU_OUT,'(A)')TRIM(STRING)
	     READ(LU_IN,'(A)')STRING
	   END DO
!
! Output summary of Teff, R, and V at RSTAR, Tau=10, and TAU=2/3.
!
	   NEXT_LOC=1  ;   STRING=' '
	   CALL WR_VAL_INFO(STRING,NEXT_LOC,'Tau',TA(ND))
	   T1=1.0D+10*R(ND)/RAD_SUN(); CALL WR_VAL_INFO(STRING,NEXT_LOC,'R*/Rsun',T1)
	   T1=TEFF_SUN()*(LUM/T1**2)**0.25
	   CALL WR_VAL_INFO(STRING,NEXT_LOC,'T*  ',T1)
	   CALL WR_VAL_INFO(STRING,NEXT_LOC,'V(km/s)',V(ND))
	   WRITE(LU_OUT,'(A)')TRIM(STRING)
!
	   NEXT_LOC=1  ;   STRING=' '
	   CALL WR_VAL_INFO(STRING,NEXT_LOC,'Tau',TB(2))		!10.0D0
	   T1=1.0D+10*TC(2)/RAD_SUN()
	   CALL WR_VAL_INFO(STRING,NEXT_LOC,'R /Rsun',T1)
	   T1=TEFF_SUN()*(LUM/T1**2)**0.25
	   CALL WR_VAL_INFO(STRING,NEXT_LOC,'Teff',T1)
	   CALL WR_VAL_INFO(STRING,NEXT_LOC,'V(km/s)',AV(2))
	   WRITE(LU_OUT,'(A)')TRIM(STRING)
!
	   NEXT_LOC=1  ;   STRING=' '
	   CALL WR_VAL_INFO(STRING,NEXT_LOC,'Tau',TB(1))		!0.67D0
	   T1=1.0D+10*TC(1)/RAD_SUN()
	   CALL WR_VAL_INFO(STRING,NEXT_LOC,'R /Rsun',T1)
	   T1=TEFF_SUN()*(LUM/T1**2)**0.25
	   CALL WR_VAL_INFO(STRING,NEXT_LOC,'Teff',T1)
	   CALL WR_VAL_INFO(STRING,NEXT_LOC,'V(km/s)',AV(1))
	   WRITE(LU_OUT,'(A)')TRIM(STRING)
!
	   READ(LU_IN,'(A)')STRING		!Tau=10 string
	   READ(LU_IN,'(A)')STRING		!Tau=1 string
!
	   IOS=0
	   READ(LU_IN,'(A)',IOSTAT=IOS)STRING
	   DO WHILE(IOS .EQ. 0)
	     WRITE(LU_OUT,'(A)')TRIM(STRING)
	     READ(LU_IN,'(A)',IOSTAT=IOS)STRING
	   END DO
	   WRITE(LU_OUT,'(A)')'R(Tau) values have been revised.'
	   CLOSE(LU_OUT); CLOSE(LU_IN)
!
	ELSE IF(XOPT .EQ. 'MEANOPAC')THEN
!
! Output Opacities, Optical depth scales and optical depth increments to file 
! in an identical format to MEANOPAC created by CMFGEN.
!
	  TCHI(1:ND)=ROSS_MEAN(1:ND)*CLUMP_FAC(1:ND)
	  CALL DERIVCHI(dCHIdR,TCHI,R,ND,METHOD)
          CALL NORDTAU(TA,TCHI,R,R,dCHIdR,ND)
	  TCHI(1:ND)=FLUX_MEAN(1:ND)*CLUMP_FAC(1:ND)
	  CALL DERIVCHI(dCHIdR,TCHI,R,ND,METHOD)
          CALL NORDTAU(TB,TCHI,R,R,dCHIdR,ND)
	  CALL ESOPAC(ESEC,ED,ND)		!Elec. scattering opacity.
	  TCHI(1:ND)=ESEC(1:ND)*CLUMP_FAC(1:ND)
	  CALL DERIVCHI(dCHIdR,TCHI,R,ND,METHOD)
          CALL NORDTAU(DTAU,TCHI,R,R,dCHIdR,ND)
	  TA(ND)=0.0D0; TB(ND)=0.0D0; DTAU(ND)=0.0D0
!
	  CALL GEN_ASCI_OPEN(LU_OUT,'MEANOPAC','UNKNOWN',' ',' ',IZERO,IOS)
	    WRITE(LU_OUT,
	1    '( ''     R        I   Tau(Ross)   /\Tau   Rat(Ross)'//
	1    '  Chi(Ross)  Chi(ross)  Chi(Flux)   Chi(es) '//
	1    '  Tau(Flux)  Tau(es)  Rat(Flux)  Rat(es)'' )' )
	    IF(R(1) .GE. 1.0D+05)THEN
	      FMT='( 1X,1P,E10.4,2X,I3,1X,1P,E9.3,2(2X,E8.2),1X,'//
	1          '4(2X,E9.3),2(2X,E8.2),2(2X,E8.2) )'
	    ELSE
	      FMT='( 1X,F10.4,2X,I3,1X,1P,E9.3,2(2X,E8.2),1X,'//
	1          '4(2X,E9.3),2(2X,E8.2),2(2X,E8.2) )'
	    END IF
	    DO I=1,ND
	      IF(I .EQ. 1)THEN
	        T1=0.0D0		!Rosseland optical depth scale
	        T2=0.0D0		!Flux optical depth scale
	        T3=0.0D0		!Elec. scattering optical depth scale.
	        TC(1:3)=0.0D0
	      ELSE
	        T1=T1+TA(I-1)
	        T2=T2+TB(I-1)
	        T3=T3+DTAU(I-1)
	        TC(1)=TA(I)/TA(I-1)
	        TC(2)=TB(I)/TB(I-1)
	        TC(3)=DTAU(I)/DTAU(I-1)
	      END IF
	      WRITE(LU_OUT,FMT)R(I),I,T1,TA(I),TC(1),
	1        ROSS_MEAN(I),ROSS_MEAN(I),FLUX_MEAN(I),ESEC(I),
	1        T2,T3,TC(1),TC(2)
	    END DO
	    WRITE(LU_OUT,'(//,A,/,A,/,A)')
	1       'NB: Mean opacities do not include effect of clumping.',
	1       'NB: Optical depth scale includes effect of clumping.',
	1       'INT_dBdT(opacity) set to ROSS_MEAN.'
	  CLOSE(UNIT=LU_OUT)
!
! Option to determine the velocity and radius of a star at a give
! Rosseland optical depth.
!
	ELSE IF(XOPT .EQ. 'CHKV')THEN
	  CALL USR_HIDDEN(T1,'TAU','100.0',' ')
!
! Determine radius where TAUROSS=T1
!
	  IF(.NOT. ROSS)THEN
	    WRITE(T_OUT,*)'Error - Roesseland optical depth not computed'
	  ELSE
	    IF(TAUROSS(ND) .LT. T1)THEN
	      T2=LOG(T1/TAUROSS(ND-2)) / LOG(TAUROSS(ND-1)/TAUROSS(ND-2))
	      NEW_RSTAR=EXP( LOG(R(ND-2)) - T2*LOG(R(ND-2)/R(ND)) )
	      T2=LOG(R(ND-2)/NEW_RSTAR) / LOG(R(ND-2)/R(ND))
	      NEW_VSTAR=EXP( LOG(V(ND-2)) + T2*LOG( V(ND)/V(ND-2) ) )
	    ELSE
	      I=ND-1
	      DO WHILE(TAUROSS(I) .GT. T1)
	        I=I-1
	      END DO
	      T2=LOG(T1/TAUROSS(I)) / LOG(TAUROSS(I+1)/TAUROSS(I))
	      NEW_RSTAR=EXP( LOG(R(I)) - T2*LOG(R(I)/R(I+1)) )
	      T2=LOG(R(I)/NEW_RSTAR) / LOG(R(I)/R(I+1))
	      NEW_VSTAR=EXP(  LOG(V(I)) + T2*LOG( V(I+1)/V(I) )  )
	    END IF
	    WRITE(T_OUT,730)T1,NEW_RSTAR,NEW_VSTAR
730	    FORMAT(1X,'Estimate of R* where Rosseland Optical depth is',
	1     F8.2,' is ',F10.4,/,
	1     1X,'Velocity at R* is ',1PE14.5,'km/s')
	  END IF
!
	ELSE IF(XOPT .EQ. 'COLR')THEN
	  IF(DPTH_INDX .LT. 1 .OR. DPTH_INDX .GT. ND)DPTH_INDX=ND/2
	  DEFAULT=WR_STRING(DPTH_INDX)
	  CALL USR_OPTION(DPTH_INDX,'Depth',DEFAULT,'Input depth to check collision quantities')
	  STRING=' '
	  WRITE(6,*)' '
	  WRITE(6,*)'Options are:'
	  WRITE(6,*)'   NR: Output net rates;       DR: Downward rates'
	  WRITE(6,*)'   CR: Output cooling terms;   CL: Cooling lines'
	  WRITE(6,*)' '
	  DEFAULT=COL_OPT
	  CALL USR_OPTION(COL_OPT,'TYPE',DEFAULT,'NR, DR, CR, CL')
	  IF(UC(COL_OPT(1:2)) .EQ. 'NR')THEN
	    STRING='NET_RATES'
	  ELSE IF(UC(COL_OPT(1:2)) .EQ. 'CR')THEN
	    STRING='COOL_RATES'
	  ELSE IF(UC(COL_OPT(1:2)) .EQ. 'CL')THEN
	    STRING='COOL_LINES'
	  ELSE
	    STRING=' '
	  END IF
	  I=DPTH_INDX
	  TMP_ED=ED(I)
	  T1=T(I)
	  DO ID=1,NUM_IONS
	    IF(XSPEC .EQ. UC(ION_ID(ID)))THEN
	      CALL SUBCOL_MULTI_V4(
	1         OMEGA_F,dln_OMEGA_F_dlnT,
	1         OMEGA_S,dln_OMEGA_S_dlnT,
	1         ATM(ID)%XzV_F(1,I),ATM(ID)%XzVLTE_F(1,I),UNIT_VEC,ATM(ID)%NXzV_F,
	1         ATM(ID)%XzV_F(1,I),ATM(ID)%XzVLTE_F(1,I),ATM(ID)%W_XzV_F(1,I),
	1         ATM(ID)%EDGEXzV_F,ATM(ID)%AXzV_F,ATM(ID)%GXzV_F,
	1         ATM(ID)%XzVLEVNAME_F,ATM(ID)%NXzV_F,ATM(ID)%ZXzV,
	1         ID,TRIM(ION_ID(ID))//'_COL_DATA',OMEGA_GEN_V3,
	1         ATM(ID)%F_TO_S_XzV,TEMP,T1,TMP_ED,IONE)
	      IF(STRING .EQ. 'COOL_LINES')THEN
	  	CALL USR_OPTION(FLAG,'NEW_FILE',FIRST_COLR,'New file?')
	        FIRST_COLR='F'
	        J=10
	        CALL WR_COL_LINES(OMEGA_S,ATM(ID)%AXzV_F,ATM(ID)%XzV_F(1,I),ATM(ID)%XzVLTE_F(1,I),
	1               ATM(ID)%W_XzV_F(1,I),ATM(ID)%GXzV_F,
	1               ATM(ID)%EDGEXzV_F,ATM(ID)%XzVLEVNAME_F,ATM(ID)%NXzV_F,
	1               TRIM(ION_ID(ID)),R(I),V(I),SIGMA(I),I,J,LU_COL,FLAG)
	      ELSE
	        CALL WR_COL_RATES(OMEGA_S,ATM(ID)%XzV_F(1,I),ATM(ID)%XzVLTE_F(1,I),
	1               ATM(ID)%EDGEXzV_F,ATM(ID)%XzVLEVNAME_F,ATM(ID)%NXzV_F,
	1               TRIM(XSPEC)//'R',LU_COL,STRING)
	      END IF
	    END IF
	  END DO
!
!	ELSE IF(XOPT .EQ. 'NONT')THEN
!          DO ID=1,NUM_IONS
!            IF(XSPEC .EQ. UC(ION_ID(ID)))THEN
!	       CALL DISP_BETHE_APPROX_V5(Q,NL,NUP,XKT,dXKT_ON_XKT,NKT,ID,DPTH_INDX)
!	    END IF
!	  END DO
!
	ELSE IF(XOPT .EQ. 'COL' .OR. XOPT .EQ. 'CRIT')THEN
	  TMP_ED=1.0D0
	  CALL USR_OPTION(T1,'T','1.0','Input T')
	  DO ID=1,NUM_IONS
	    IF(XSPEC .EQ. UC(ION_ID(ID)))THEN
	      CALL SUBCOL_MULTI_V4(
	1         OMEGA_F,dln_OMEGA_F_dlnT,
	1         OMEGA_S,dln_OMEGA_S_dlnT,
	1         UNIT_VEC,UNIT_VEC,ZERO_VEC,ATM(ID)%NXzV_F,
	1         UNIT_VEC,UNIT_VEC,UNIT_VEC,
	1         ATM(ID)%EDGEXzV_F,ATM(ID)%AXzV_F,ATM(ID)%GXzV_F,
	1         ATM(ID)%XzVLEVNAME_F,ATM(ID)%NXzV_F,ATM(ID)%ZXzV,
	1         ID,TRIM(ION_ID(ID))//'_COL_DATA',OMEGA_GEN_V3,
	1         ATM(ID)%F_TO_S_XzV,TEMP,T1,TMP_ED,IONE)
	          IF(XOPT .EQ. 'COL')THEN
	            CALL WR_COL(OMEGA_F,ATM(ID)%XzVLEVNAME_F,ATM(ID)%NXzV_F,XSPEC,LU_COL,' ')
	          ELSE
	            CALL WR_CRIT(OMEGA_F,T1,HDKT,ATM(ID)%AXzV_F,ATM(ID)%EDGEXzV_F,ATM(ID)%GXzV_F,
	1                    ATM(ID)%XzVLEVNAME_F,ATM(ID)%NXzV_F,XSPEC,LU_COL)
	         END IF
	    END IF
	  END DO
!
! 
!
	ELSE IF(XOPT .EQ. 'RNT')THEN
	  DEFAULT='NON_THERM_SPEC_INFO'
	  CALL USR_OPTION(FILENAME,'FILE',DEFAULT,'File name with fractional energies')
	  CALL RD_NON_THERM_SPEC(YV,ZV,WV,ND,FILENAME,LU_IN,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error unable to open file. IOS=',IOS
	    GOTO 1
	  END IF
!
	  WRITE(T_OUT,'(A)')' '
	  WRITE(T_OUT,'(A)')'Curve ordering is Felec, Fion, Fexec '
	  CALL DP_CURVE(ND,XV,YV)
	  CALL DP_CURVE(ND,XV,ZV)
	  CALL DP_CURVE(ND,XV,WV)
	  YAXIS='Channel fraction'
!
	ELSE IF(XOPT .EQ. 'PHOT')THEN
	  CALL USR_OPTION(TEMP,'Nu','0.0','Input frequency (10^15)Hz')
	  CALL USR_OPTION(PHOT_ID,'Nu','1','Photoionization route')
	  FLAG=.FALSE.				!Don't return edge value.
	  IF(TEMP .EQ. 0.0)FLAG=.TRUE.
	  I=1
!
	  K=0
	  OMEGA_F(:,:)=0.0D0
	  DO ID=1,NUM_IONS
	    IF(XSPEC .EQ. UC(ION_ID(ID)))THEN
	      K=ATM(ID)%NXzV_F
	      CALL SUB_PHOT_GEN(ID,OMEGA_F,TEMP,ATM(ID)%EDGEXzV_F,
	1            ATM(ID)% NXzV_F,PHOT_ID,FLAG)
	      EXIT
	    END IF
	 END DO
	 IF(K .EQ. 0)THEN
	    WRITE(T_OUT,*)'Species not recognized'
	    GOTO 1
	  END IF
!
	  WRITE(T_OUT,*)'Cross sections (Mbarns) for frequency ',TEMP
	  DO I=1,K
	    WRITE(T_OUT,'(1X,I4,3X,1P,E12.5)')I,OMEGA_F(I,1)/1.0D-08
	  END DO
	  WRITE(LU_CROSS,*)'Cross sections (Mbarns) for frequency ',TEMP
	  DO I=1,K
	    WRITE(LU_CROSS,'(1X,I4,3X,1P,E12.5)')I,OMEGA_F(I,1)/1.0D-08
	  END DO
! 
!  
	ELSE IF(XOPT .EQ. 'PLTPHOT')THEN
	  CALL PLTPHOT_SUB(XSPEC,XV,YV,WV,ZV,N_PLT_MAX,OMEGA_F,N_MAX,XAXIS,YAXIS,XRAYS,VSM_DIE_KMS,ND)
!
	ELSE IF(XOPT .EQ. 'PLTARN')THEN
	  CALL PLT_ARN(XSPEC,ND,XV,YV,N_PLT_MAX)
	  XAXIS='eV'; YAXIS='\gs(cm\u-2\d)'	
!
	ELSE IF(XOPT .EQ. 'RDDIE')THEN
	  CALL USR_OPTION(DIE_REG,'REG','F','Include dielectronic lines')
	  CALL USR_OPTION(DIE_WI,'WI','F','Include WI dielectronic lines')
	  CALL USR_OPTION(VSM_DIE_KMS,'VSM','3000.0D0',
	1                                   'Smoothing width in km/s')
	  DO ID=1,NUM_IONS
	    IF(XSPEC .EQ. UC(ION_ID(ID)))THEN
	      IF(DIE_REG .OR. DIE_WI)THEN
	        CALL RD_PHOT_DIE_V1(ID,ATM(ID)%EDGEXzV_F,ATM(ID)%XzVLEVNAME_F,
	1             ATM(ID)%NXzV_F,ATM(ID)%GIONXzV_F,
	1             VSM_DIE_KMS,DIE_REG,DIE_WI,
	1             ION_ID(ID),LU_IN,LU_LOG,'DIE'//TRIM(ION_ID(ID)))
	      END IF
	    END IF
	  END DO
!  
c 
!
! Section to compute the recombination rates for various species.
! Default is for ionizatio to the ground state. The photionzation route
! is a hidden parameters. Both GION and EXC_EN must be specified. These
! are input to accomodate TERM spitting.
!
	ELSE IF(XOPT .EQ. 'NRR')THEN
	  CALL USR_OPTION(T1,'T','1.0','Input T (in 10^4 K)')
	  CALL USR_HIDDEN(PHOT_ID,'PHOT_ID','1',
	1      'Photoionization route')
	  EXC_EN=0.0D0
	  IF(PHOT_ID .NE. 1)THEN
	    CALL USR_OPTION(EXC_EN,'EXC_EN',' ',
	1     'Excitaiton Energy (cm^-1) of final state')
	  END IF
	  EXC_EN=1.0D-15*C_CMS*EXC_EN
	  CALL USR_OPTION(TMP_GION,'GION',' ',
	1      'G for ION (No def)')
	  J=0		!Set in case no identification.
	  I=INDEX(X,' ')
!
	  DO ID=1,NUM_IONS
	    I=ID
	    IF(XSPEC .EQ. UC(ION_ID(ID)) .AND. ATM(ID)%XzV_PRES)THEN
	      CALL RECOM_CHK_V2(TA,ATM(ID)%EDGEXzV_F,ATM(ID)%GXzV_F,TMP_GION,
	1            ATM(ID)% NXzV_F,EXC_EN,PHOT_ID,SUB_PHOT_GEN,I,T1)
	      STRING=ION_ID(ID)
              J=ATM(ID)%NXzV_F
	      EXIT
	    END IF
	  END DO
!
	  IF(J .NE. 0)THEN
	    WRITE(LU_REC,*)' '
	    WRITE(LU_REC,*)'Temperatre (10^4 K) is',T1
	    WRITE(T_OUT,*)'Recombination rates for ',TRIM(STRING)
	    WRITE(LU_REC,*)'Recombination rates for ',TRIM(STRING)
	    WRITE(T_OUT,'(1X,A,1X,I3)')'Photoionization route number is',PHOT_ID
	    WRITE(LU_REC,'(1X,A,1X,I3)')'Photoionization route number is',PHOT_ID
	    WRITE(T_OUT,*)' '
	    WRITE(LU_REC,*)' '
	    T1=0.0D0
	    DO I=1,J
	      WRITE(T_OUT, '(I4,ES12.4,2X,A)')I,TA(I),TRIM(ATM(ID)%XzVLEVNAME_F(I))
	      WRITE(LU_REC,'(I4,ES12.4,2X,A)')I,TA(I),TRIM(ATM(ID)%XzVLEVNAME_F(I))
	      T1=T1+TA(I)
	    END DO
	    WRITE(T_OUT,'(1X,A,1X,1PE12.4)')'Total recombination rate is:',T1
	    WRITE(T_OUT,'(1X,A,1X,1PE10.2)')'Value for GION used was:',TMP_GION
	    WRITE(LU_REC,'(1X,A,1X,1PE12.4)')'Total recombination rate is:',T1
	    WRITE(LU_REC,'(1X,A,1X,1PE10.2)')'Value for GION used was:',TMP_GION
	  END IF
!
! 
!
	ELSE IF(XOPT .EQ. 'GNTFF')THEN
	    CALL USR_OPTION(CHI,3,3,'GNT',' ','LAM(um) , T & Z')
	    CHI(1)=1.0D-04*ANG_TO_HZ/CHI(1)
	    CHI(4)=GFF(CHI(1),CHI(2),CHI(3))
	    WRITE(T_OUT,2000)CHI(4)
2000	    FORMAT(3X,'The Free-free gaunt factor is',1X,1PE10.3)
!
	ELSE IF(XOPT .EQ. 'GNTBF')THEN
	  CALL USR_OPTION(CHI,3,3,'GNT',' ','Lam(um) ,Level & Z')
	    CHI(1)=1.0D-04*ANG_TO_HZ/CHI(1)
	    I=(CHI(2)+0.000002)
	    CHI(4)=GBF(CHI(1),I,CHI(3))
	    WRITE(T_OUT,2001)CHI(4)
2001	    FORMAT(3X,'The Bound-free gaunt factor is',1X,1PE10.3)
!
	ELSE IF(XOPT .EQ. 'LAM')THEN
!
! NB. FL and FREQ have been equivalenced.
!
	  DEFAULT=WR_STRING(FREQ)
	  CALL USR_HIDDEN(T1,'FREQ',DEFAULT,' ')
	  T2=LAMVACAIR(T1)		!Wavelength(Angstroms)
	  WRITE(T_OUT,'(1X,''Lambda(in air)='',1PE14.6)')T2
	  WRITE(T_OUT,'(1X,''Lambda(in vac)='',1PE14.6)')ANG_TO_HZ/T1
	ELSE IF(XOPT .EQ. 'TIT')THEN
	  CALL USR_OPTION(NAME,'Title',' ','Title for all graphs')
	ELSE IF(XOPT .EQ. 'METHOD')THEN
	  CALL USR_OPTION(METHOD,'LOGLOG',' ','Tau option (LOGLOG,LOGMON, ZERO')
!
! 
!
! SETREC sets TA to DCI if DCI is non-zero, otherwise TA is set to C2(1, ).
! Thus can compute C2-CI recombination rate even if CI is not present.
!
	ELSE IF(XOPT .EQ. 'RR')THEN
	   FOUND=.FALSE.
	   DO ID=1,NUM_IONS
	     IF(XSPEC .EQ. ION_ID(ID))THEN
	       CALL SETREC(TA,ATM(ID)%DXzV,ATM(ID+1)%XzV_F,ATM(ID+1)%NXzV_F,ND)
	       FOUND=.TRUE.
	     END IF
	   END DO
	   IF(.NOT. FOUND)THEN
	     WRITE(T_OUT,*)'Invalid Recombination Request'
	     GOTO 3
	   END IF
!
! Compute the recombination rate. Two recombination rates will be
! computed. The first will ignore the temperature variation, whilst the
! second will normalize the number of recombinations to 10^4K by assuming
! that the recombintion rate goes as T^{-0.8}. For d=1kpc, alpha=10^{-12}
! the recombination rate will have the units ergs/cm^2/s . (need to
! divide bu the transition wavelength in mum)
!
	   CALL USR_HIDDEN(ELEC,'TVAR','T','Include T variation?')
	   CALL USR_HIDDEN(T2,'EXP','0.8',
	1	'Exponent for T variation?')
!
	   T1=-26.68059
!        
	   DO I=1,ND
	     TA(I)=(ED(I)/1.0D+10)*TA(I)*( R(I)**3 )
	     TB(I)=TA(I)*( T(I)**(-T2) )
	   END DO
	   IF(ELEC)THEN
!
! Need to use TC to compute YV as can get floating point underflow
! at outer boundary.
!
	     TC(1)=TB(1)
	     DO I=2,ND
	       TC(I)=TC(I-1)+(TB(I)+TB(I-1))*
	1            ( LOG(R(I-1))-LOG(R(I)) )*0.5
	       YV(I-1)=LOG10(TC(I-1))+T1
	     END DO
	     YV(ND)=LOG10(TC(ND))+T1
	     CALL DP_CURVE(ND,XV,YV)
	      YAXIS='Log(c.Ne.X\uN+\d.T\u0.8\d)'
	   ELSE
	     TC(1)=TA(1)
	     DO I=2,ND
	       TC(I)=TC(I-1)+(TA(I)+TA(I-1))*
	1            ( LOG(R(I-1))-LOG(R(I)) )*0.5
	       YV(I-1)=LOG10(TC(I-1))+T1
	     END DO
	     YV(ND)=LOG10(TC(ND))+T1
	     CALL DP_CURVE(ND,XV,YV)
	     YAXIS='Log(c.Ne.X\uN+\d)'
	   END IF
!
! 
!
! Plots the Rayleigh scattering opacity as a function of depth.
!
	ELSE IF(XOPT .EQ. 'RAY')THEN
	  CALL USR_OPTION(FREQ,'LAM','0.0',FREQ_INPUT)
	  IF(FREQ .EQ. 0)THEN
	  ELSE IF(KEV_INPUT)THEN
	      FREQ=FREQ*KEV_TO_HZ
	  ELSE IF(ANG_INPUT)THEN
	      FREQ=ANG_TO_HZ/FREQ
	  END IF
          IF(ATM(1)%XzV_PRES)THEN
	    CHI_RAY(1:ND)=0.0D0
            CALL RAYLEIGH_SCAT(CHI_RAY,ATM(1)%XzV_F,ATM(1)%AXzV_F,ATM(1)%EDGEXZV_F,
	1             ATM(1)%NXzV_F,FREQ,ND)
	    CHI_RAY(1:ND)=1.0D-10*CHI_RAY(1:ND)*CLUMP_FAC(1:ND)
          END IF
	  CALL DP_CURVE(ND,XV,CHI_RAY)
	  YAXIS='\gx(cm\u-1\d)'
!
! Plots the Rayleigh scattering cross-section as a function of opacity.
! In units of the Thompson-cross section.
!
	ELSE IF(XOPT .EQ. 'RAYL')THEN
	  CALL USR_OPTION(LAM_ST,'LAM_ST','100.0',FREQ_INPUT)
	  CALL USR_OPTION(LAM_EN,'LAM_EN','10000.0',FREQ_INPUT)
	  CALL USR_OPTION(K,'NPNTS=','1000','Number of points')
	  K=MIN(K,N_PLT_MAX)
	  T2=EXP(LOG(LAM_EN/LAM_ST)/(K-1))
	  T1=LAM_ST
	  I=0
	  DO WHILE(I .LT. K)
	    I=I+1
	    FREQ=ANG_TO_HZ/T1
	    WV(I)=T1
            IF(ATM(1)%XzV_PRES)THEN
	      YV(I)=0.0D0
              CALL RAYLEIGH_SCAT(YV(I),ATM(1)%XzV_F,ATM(1)%AXzV_F,ATM(1)%EDGEXZV_F,
	1             ATM(1)%NXzV_F,FREQ,IONE)
	    ELSE
	      GOTO 1
            END IF
	    T1=T1*T2
	  END DO
	  YV(1:K)=YV(1:K)/ATM(1)%XzV_F(1,1)/6.65D-15
	  CALL DP_CURVE(K,WV,YV)
	  XAXIS='\gl(\A)'
	  YAXIS='\gs/\gs\dT\u'
!
	ELSE IF(XOPT .EQ. 'ETA')THEN
	  ELEC=.FALSE.
	  CALL USR_OPTION(ELEC,'RCUBE','F','Multiply ETA by R^3?')
	  IF(ELEC)THEN
	    DO I=1,ND
	      YV(I)=DLOG10(R(I)*R(I)*R(I)*ETA(I)+1.0D-250)+20.0D0
	    END DO
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='Log(r\u3\d\ge(ergs/s/Hz)'
	  ELSE
	    DO I=1,ND
	      YV(I)=DLOG10(ETA(I)+1.0D-250)-10.0D0
	    END DO
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='Log(\ge(ergs/cm\u3\d/s/Hz)'
	  END IF
!
	ELSE IF(XOPT .EQ. 'DIFFT')THEN
	  ELEC=.FALSE.
	  CALL USR_OPTION(ELEC,'DB','T','Diffusion time to outer boundary?')
	  DO I=1,ND
	    CHIROSS(I)=CLUMP_FAC(I)*ROSS_MEAN(I)
	  END DO
	  CALL TORSCL(TAUROSS,CHIROSS,R,TB,TC,ND,METHOD,TYPE_ATM)
	  WRITE(T_OUT,*)'Rossland optical depth is : ',TAUROSS(ND)
	  DO I=1,ND-1
	    YV(I)=(1.0D0+TAUROSS(I+1)-TAUROSS(I))*(R(I)-R(I+1))*1.0D+05/C_KMS/24.0D0/3600.0D0
	  END DO
	  J=ND-1
	  IF(ELEC)THEN
	    ZV(1)=0.0D0
	    DO I=2,ND
	      ZV(I)=ZV(I-1)+YV(I-1) 
	    END DO
	    CALL DP_CURVE(J,XV,ZV)
	    YAXIS='Diffusion time to outer boundary (days)'
	  ELSE
	    WRITE(6,*)'Computed diffusion time between grid points'
	    CALL DP_CURVE(J,XV,YV)
	    YAXIS='Diffusion time (days)'
	  END IF
!
	ELSE IF(XOPT .EQ.'OP' .OR.
	1       XOPT .EQ. 'KAPPA' .OR.
	1       XOPT .EQ. 'TAUC' .OR.
	1       XOPT .EQ. 'DTAUC') THEN
!
	  IF(.NOT. ELEC)THEN
	    DO I=1,ND
	      CHI(I)=CHI(I)-ESEC(I)
	    END DO
	  END IF
!
	  IF(XOPT .EQ. 'TAUC')THEN
	    CALL TORSCL(TA,CHI,R,TB,TC,ND,METHOD,TYPE_ATM)
	    DO I=1,ND
	      YV(I)=DLOG10(TA(I))
	    END DO
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='Log(\gt)'
	  ELSE IF(XOPT .EQ. 'DTAUC')THEN
	    CALL TORSCL(TA,CHI,R,TB,TC,ND,METHOD,TYPE_ATM)
	    DO I=1,ND-1
	      YV(I)=DLOG10(TA(I+1)-TA(I))
	    END DO
	    I=ND-1
	    CALL DP_CURVE(I,XV,YV)
	    YAXIS='Log(\gD\gt)'
	  ELSE IF(XOPT .EQ. 'KAPPA')THEN
	     YV(1:ND)=1.0D-10*CHI(1:ND)/MASS_DENSITY(1:ND)/CLUMP_FAC(1:ND)
	     CALL DP_CURVE(ND,XV,YV)
	     YAXIS='\gk(cm\u3 \d/g)'
	  ELSE
!
! We subtract 10 to put CHI in units of cm-1
!
	    WRITE(6,*)'Warning: The plotted CHI has units of cm^-1 - not program units'
	    DO I=1,ND
	      YV(I)=DLOG10(CHI(I))-10.0
	    END DO
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='Log \gx(cm\u-1\d)'
	  END IF
!
	ELSE IF(XOPT .EQ. 'TWOOP')THEN
	  CALL USR_OPTION(LAM_ST,'LAMST','911.0D0',FREQ_INPUT)
	  CALL USR_OPTION(LAM_EN,'LAMEN','9000.0D0',FREQ_INPUT)
	  IF(KEV_INPUT)THEN
	    LAM_ST=LAM_ST*KEV_TO_HZ
	    LAM_EN=LAM_EN*KEV_TO_HZ
	  ELSE IF(ANG_INPUT)THEN
	    LAM_ST=ANG_TO_HZ/LAM_ST
	    LAM_EN=ANG_TO_HZ/LAM_EN
	  END IF
	  CALL USR_OPTION(NFREQ,'NPTS','-100','+ve lin. spacing, -ve log')
!
	  IF(NFREQ .LT. 0)THEN
	    LINX=.FALSE.
	    NFREQ=-NFREQ
	    T1=LOG(LAM_EN/LAM_ST)/(NFREQ-1)
	    DO I=1,NFREQ
	      XNU(I)=EXP( LOG(LAM_ST)+(I-1)*T1 )
	    END DO
	  ELSE
	    LINX=.TRUE.
	    T1=(LAM_EN-LAM_ST)/(NFREQ-1)
	    DO I=1,NFREQ
	      XNU(I)=LAM_ST+(I-1)*T1
	    END DO
	  END IF
!
	  DO I=1,NUM_IONS
	    IF(XSPEC .EQ. UC(ION_ID(I)))THEN
	      ID=I
	      EXIT
	    END IF
	  END DO
	  DEFAULT='LTE'
	  CALL USR_OPTION(TWO_PHOT_OPTION,'OPT','LTE','Two photon option: LTE, NOSTIM, RAD, and OLD_DEF')
!
	  XAXSAV=XAXIS
	  CALL USR_OPTION(ELEC,'OPAC','F','Plot opacity instead of emissivity?')
	  DO ML=1,NFREQ
	    ETA(1:ND)=0.0D0; CHI(1:ND)=0.0D0
	    CALL TWO_PHOT_OPAC_DISP_V3(ETA,CHI,ATM(ID)%XzV_F,T,XNU(ML),TWO_PHOT_OPTION,ID,ND,ATM(ID)%NXzV_F)
	    IF(ELEC)THEN
	      YV(ML)=CHI(1)
	      WV(ML)=CHI(ND/2)
	      ZV(ML)=CHI(ND)
	    ELSE
	      YV(ML)=ETA(1)
	      WV(ML)=ETA(ND/2)
	      ZV(ML)=ETA(ND)
	    END IF
	  END DO
	  IF(ELEC)THEN
	    YAXIS='Opacity'
	  ELSE
	    YAXIS='Emissivity'
	  END IF
!
	  CALL USR_OPTION(ELEC,'PHOT','F','Plot photon emission rate instead of emissivity?')
	  IF(ELEC)THEN
	    DO ML=1,NFREQ
	      T1=4.0D0*PI/6.626D-27/1.0D+15/XNU(ML)/1.0D+10
	      YV(ML)=YV(ML)*T1/ATM(ID)%XzV_F(2,1)
	      WV(ML)=WV(ML)*T1/ATM(ID)%XzV_F(2,ND/2)
	      ZV(ML)=ZV(ML)*T1/ATM(ID)%XzV_F(2,ND)
	    END DO
	    YAXIS='# of photons/s/Hz'
	  END IF
!
	  XAXIS='\gv(10\u15\d Hz)'
	  CALL DP_CURVE(NFREQ,XNU,YV)
	  CALL DP_CURVE(NFREQ,XNU,WV)
	  CALL DP_CURVE(NFREQ,XNU,ZV)
!
! Estimate the effective absorbative optical depth scale for Gamma-rays.
!
	ELSE IF(XOPT .EQ. 'TAUGAM' .OR. XOPT .EQ. 'GAMABS')THEN
	  TA(1:ND)=0.5D0			!Number of electrons per baryon
	  DO ISPEC=1,NSPEC
	    IF('HYD' .EQ. SPECIES(ISPEC))THEN
	      TA(1:ND)=0.5D0*(1.0D0+POPDUM(1:ND,ISPEC)/POP_ATOM(1:ND))
	    END IF
	  END DO
	  XM(1:ND)=0.06D0*TA(1:ND)*MASS_DENSITY(1:ND)*1.0D+10
	  WRITE(6,*)XM(ND)*R(ND),XM(ND)*(R(ND-1)-R(ND))
!
	  IF(XOPT .EQ. 'GAMABS')THEN
	    CALL GAM_ABS(R,V,XM,ND)
	    YAXIS='Rad. Dec. Energy (ergs/cm\u2\d/sec)'
!
	  ELSE
	    CALL TORSCL(TA,XM,R,TB,TC,ND,METHOD,TYPE_ATM)
	    CALL USR_OPTION(ELEC,'dTAU','F','Plot dTAU)')
	    IF(ELEC)THEN
	      DO I=1,ND-1
	        YV(I)=LOG10(TA(I+1)/TA(I))
	      END DO
	      CALL DP_CURVE(ND-1,XV,YV)
	      YAXIS='dLog \gt\d\g\u'
	    ELSE
	      DO I=1,ND
	        YV(I)=DLOG10(TA(I))
	      END DO
	      CALL DP_CURVE(ND,XV,YV)
	      YAXIS='Log(\gt\d\g\u)'
	    END IF
	  END IF
!
! These options allow you to plot tau at a particular R (TAUR), or
! alternatively, R at a given value of Tau (RTAU).
!
! We can also write out the continuum opacities as a function of
! wavelength.
!
! To be read TAU_at_R and R_at_TAU respectively.
!
	ELSE IF(XOPT .EQ. 'RTAU'  .OR. XOPT .EQ. 'VTAU' .OR. XOPT .EQ. 'TAUR' .OR. 
	1       XOPT .EQ. 'EDTAU' .OR. 
	1       XOPT .EQ. 'KAPR'  .OR. XOPT .EQ. 'CHIR' .OR. XOPT .EQ. 'ALBEDO' .OR.
	1       XOPT .EQ. 'ETAR'  .OR. XOPT .EQ. 'WROPAC')THEN
	  IF(XRAYS)WRITE(T_OUT,*)'Xray opacities (i.e. K shell) are included'
	  IF(.NOT. XRAYS)WRITE(T_OUT,*)'Xray opacities (i.e. K shell) are NOT included'
	  CALL USR_OPTION(LAM_ST,'LAMST',' ',FREQ_INPUT)
	  CALL USR_OPTION(LAM_EN,'LAMEN',' ',FREQ_INPUT)
	  IF(KEV_INPUT)THEN
	    LAM_ST=LAM_ST*KEV_TO_HZ
	    LAM_EN=LAM_EN*KEV_TO_HZ
	  ELSE IF(ANG_INPUT)THEN
	    LAM_ST=ANG_TO_HZ/LAM_ST
	    LAM_EN=ANG_TO_HZ/LAM_EN
	  END IF
	  CALL USR_OPTION(NFREQ,'NPTS','10','+ve lin. spacing, -ve log')
!
	  IF(NFREQ .LT. 0)THEN
	    LINX=.FALSE.
	    NFREQ=-NFREQ
	    T1=LOG(LAM_EN/LAM_ST)/(NFREQ-1)
	    DO I=1,NFREQ
	      Q(I)=EXP( LOG(LAM_ST)+(I-1)*T1 )
	      ZV(I)=Q(I) 			!Use ZV .: don't corrupt XV
	    END DO
	  ELSE
	    LINX=.TRUE.
	    T1=(LAM_EN-LAM_ST)/(NFREQ-1)
	    DO I=1,NFREQ
	      Q(I)=LAM_ST+(I-1)*T1
	      ZV(I)=Q(I) 			!Use ZV .: don't corrupt XV
	    END DO
	  END IF
	  CALL USR_OPTION(ELEC,'ELEC','T','Include elec?')
!
	  CALL USR_HIDDEN(LINX,'LINX','F','Linear X spacing?')
	  IF(KEV_INPUT)THEN
	    DO I=1,NFREQ
	      ZV(I)=ZV(I)/KEV_TO_HZ
	    END DO
	    XAXIS='E(keV)'
	  ELSE IF(ANG_INPUT)THEN
	    DO I=1,NFREQ
	      ZV(I)=ANG_TO_HZ/ZV(I)
	    END DO
	    XAXIS='\gl(\A)'
	  ELSE
	    XAXIS='\gn(10\u15\d Hz)'
	  END IF
	  IF(.NOT. LINX)THEN
	    XAXIS='Log '//XAXIS
	    DO I=1,NFREQ
	      ZV(I)=LOG10(ZV(I))
	    END DO
	  END IF
	  CALL USR_HIDDEN(LINY,'LINY','F','Linear Y Axis')
!
	  IF(XOPT .EQ. 'WROPAC')THEN
	    WRITE(6,*)' '
	    WRITE(6,*)'Writing contiinuum opacities to CONT_OPAC '
	    WRITE(6,*)' '
	    OPEN(LU_OUT,FILE='CONT_OPAC',STATUS='UNKNOWN')
	      WRITE(LU_OUT,'(I5,T30,A)')ND,'!Number of depth points'
	      WRITE(LU_OUT,'(I5,T30,A)')NFREQ,'!Number of frequency points'
	      WRITE(LU_OUT,'(A)')' '
	      WRITE(LU_OUT,'(A)')' Radius grid'
	      WRITE(LU_OUT,'(6ES14.6)')R(1:ND)
	      WRITE(LU_OUT,'(A)')' '
	      WRITE(LU_OUT,'(A)')' Velocity(km/s)'
	      WRITE(LU_OUT,'(6ES14.6)')V(1:ND)
	      WRITE(LU_OUT,'(A)')' '
	      WRITE(LU_OUT,'(A)')' Density (gm/cm^3)'
	      WRITE(LU_OUT,'(6ES14.6)')MASS_DENSITY(1:ND)
	      WRITE(LU_OUT,'(A)')' '
	      WRITE(LU_OUT,'(A)')' Clumping factor'
	      WRITE(LU_OUT,'(6ES14.6)')CLUMP_FAC(1:ND)
	      WRITE(LU_OUT,'(A)')' '
	      WRITE(LU_OUT,'(A)')' Kappa Table (cm^2/gm)'
	  ELSE IF(XOPT .EQ. 'TAUR' .OR. XOPT .EQ. 'CHIR' .OR. XOPT .EQ. 'KAPR' .OR. XOPT .EQ. 'ETAR')THEN
	    CALL USR_OPTION(RVAL,'RAD',' ','Radius in R* (-ve for depth index)')
	    IF(RVAL .GT. 0)THEN
	      RVAL=RVAL*R(ND)
	      IF(RVAL .GT. R(1))RVAL=R(1)
	      IF(RVAL .LT. R(ND))RVAL=R(ND)
	      R_INDX=ND
	      DO WHILE(RVAL .GT. R(R_INDX))
	        R_INDX=R_INDX-1
	      END DO
	      R_INDX=MIN(R_INDX,ND-1)
	    ELSE
	      R_INDX=NINT(ABS(RVAL))
	      R_INDX=MAX(1,R_INDX)
	      R_INDX=MIN(R_INDX,ND)
	      RVAL=R(R_INDX)
	      R_INDX=MIN(R_INDX,ND-1)
	    END IF
	    WRITE(6,*)'   RVAL is',RVAL
	    WRITE(6,*)'RVAL/R* is',RVAL/R(ND)
	    YAXIS='Log(\gt[R])'
	    IF(LINY)YAXIS='\gt[R]'
	    IF(XOPT .EQ. 'CHIR')YAXIS='\gx (cm\u-1\d)'
	    IF(XOPT .EQ. 'KAPR')YAXIS='\gk (cm\u2\d/g)'
	    IF(XOPT .EQ. 'ETAR')YAXIS='\ge (ergs cm\u3 \ds\u-1\d)'
	  ELSE IF(XOPT .EQ. 'VTAU')THEN
	    CALL USR_OPTION(TAU_VAL,'TAU',' ','Tau value for which V is to be determined')
	    YAXIS='V(km/s)'
	  ELSE IF(XOPT .EQ. 'EDTAU')THEN
	    CALL USR_OPTION(TAU_VAL,'TAU',' ','Tau value for which Ne is to be determined')
	    YAXIS='Ne(/cm\u3\d)'
	  ELSE
	    CALL USR_OPTION(TAU_VAL,'TAU',' ','Tau value for which R is to be determined')
	    YAXIS='Log(R[\gt]/R\d*\u)'
	    IF(LINY)YAXIS='R[\gt]/R\d*\u'
	  END IF
!
	  CALL USR_HIDDEN(IN_R_SUN,'IN_R_SUN','F','Y In R_SUN?')
	  IF(IN_R_SUN)THEN
	    IF(LINY)YAXIS='R/R\dO\u'
	    IF(.NOT. LINY)YAXIS='Log(R/R\dO\u)'
	  END IF
!
	  DO ML=1,NFREQ
!
! Compute continuum opacity and emissivity at the line frequency.
!
	    FL=Q(ML)
	    INCLUDE 'OPACITIES.INC'
!	    INCLUDE 'HOME:[jdh.cmf.carb.test]OPACITIES.INC'
!
	    IF(.NOT. ELEC)THEN
	      DO I=1,ND
	        CHI(I)=CHI(I)-ESEC(I)
	      END DO
	    END IF
!
	    IF(XOPT .EQ. 'WROPAC')THEN
	      WRITE(LU_OUT,'(/,ES16.6)')ANG_TO_HZ/FL
	      WRITE(LU_OUT,'(6ES14.4)')1.0D-10*CHI(1:ND)/MASS_DENSITY(1:ND)
	    ELSE IF(XOPT .EQ. 'ETAR')THEN
	      T1=(R(R_INDX)-RVAL)/(R(R_INDX)-R(R_INDX+1))
	      T2=ETA(R_INDX+1)*CLUMP_FAC(R_INDX+1)
	      T3=ETA(R_INDX)*CLUMP_FAC(R_INDX)
	      YV(ML)=1.0D-10*( T1*T2 + (1.0-T1)*T3 )
	      IF(.NOT. LINY)YV(ML)=LOG10(YV(ML))
	    ELSE IF(XOPT .EQ. 'CHIR')THEN
	      T1=(R(R_INDX)-RVAL)/(R(R_INDX)-R(R_INDX+1))
	      T2=CHI(R_INDX+1)*CLUMP_FAC(R_INDX+1)
	      T3=CHI(R_INDX)*CLUMP_FAC(R_INDX)
	      YV(ML)=1.0D-10*( T1*T2 + (1.0-T1)*T3 )
	      IF(.NOT. LINY)YV(ML)=LOG10(YV(ML))
	    ELSE IF(XOPT .EQ. 'KAPR')THEN
	      T1=(R(R_INDX)-RVAL)/(R(R_INDX)-R(R_INDX+1))
	      T2=CHI(R_INDX+1)/MASS_DENSITY(R_INDX+1)
	      T3=CHI(R_INDX)/MASS_DENSITY(R_INDX)
	      YV(ML)=1.0D-10*( T1*T2 + (1.0-T1)*T3 )
	      IF(.NOT. LINY)YV(ML)=LOG10(YV(ML))
	    ELSE 
!
! Adjust opacities for the effect of clumping.
!
	      DO I=1,ND
	        ETA(I)=ETA(I)*CLUMP_FAC(I)
	        CHI(I)=CHI(I)*CLUMP_FAC(I)
	        ESEC(I)=ESEC(I)*CLUMP_FAC(I)
	        CHI_RAY(I)=CHI_RAY(I)*CLUMP_FAC(I)
	        CHI_SCAT(I)=CHI_SCAT(I)*CLUMP_FAC(I)
	      END DO
!
	      CALL TORSCL(TA,CHI,R,TB,TC,ND,METHOD,TYPE_ATM)
	      IF(XOPT .EQ. 'ALBEDO')THEN
	        I=1
	        DO WHILE(TAU_VAL .GT. TA(I) .AND. I .LT. ND)
	          I=I+1
	        END DO
	        IF(TAU_VAL .GT. TA(ND-1))THEN
	          YV(ML)=ESEC(ND)/CHI(ND)
	        ELSE IF(TAU_VAL .LE. TA(1))THEN
	          YV(ML)=ESEC(1)/CHI(1)
	        ELSE
	          T2=(TA(I)-TAU_VAL)/(TA(I)-TA(I-1))
	          YV(ML)=( (1.0-T2)*ESEC(I)+T2*ESEC(I-1) )/( (1.0-T2)*CHI(I)+T2*CHI(I-1) )
	        END IF
	        YAXIS='Albedo'
	      ELSE
	        IF(XOPT .EQ. 'TAUR')THEN
	          T2=(R(R_INDX)-RVAL)/(R(R_INDX)-R(R_INDX+1))
	          YV(ML)=T2*TA(R_INDX+1) + (1.0-T2)*TA(R_INDX)
	          IF(.NOT. LINY)YV(ML)=LOG10(YV(ML))
	        ELSE IF(XOPT .EQ. 'VTAU' .OR. XOPT .EQ. 'EDTAU')THEN
	          TB(1:ND)=V(1:ND)
	          IF(XOPT .EQ. 'EDTAU')TB(1:ND)=ED(1:ND)
	          I=1
	          DO WHILE(TAU_VAL .GT. TA(I) .AND. I .LT. ND)
	            I=I+1
	          END DO
	          IF(TAU_VAL .GT. TA(ND-1))THEN
	            YV(ML)=TB(ND)
	          ELSE IF(TAU_VAL .LE. TA(1))THEN
	            YV(ML)=TB(1)*TA(1)/TAU_VAL
	          ELSE
	            T2=(TA(I)-TAU_VAL)/(TA(I)-TA(I-1))
                    YV(ML)=( (1.0-T2)*TB(I)+T2*TB(I-1) )
	          END IF
	        ELSE
	          I=1
	          DO WHILE(TAU_VAL .GT. TA(I) .AND. I .LT. ND)
	            I=I+1
	          END DO
	          IF(TAU_VAL .GT. TA(ND-1))THEN
	            YV(ML)=1.0
	          ELSE IF(TAU_VAL .LE. TA(1))THEN
	            YV(ML)=R(1)*TA(1)/TAU_VAL/R(ND)
	          ELSE
	            T2=(TA(I)-TAU_VAL)/(TA(I)-TA(I-1))
                    YV(ML)=( (1.0-T2)*R(I)+T2*R(I-1) )/R(ND)
	          END IF
	          IF(IN_R_SUN)YV(ML)=YV(ML)*R(ND)/6.96
	          IF(.NOT. LINY)YV(ML)=LOG10(YV(ML))
	        END IF
	      END IF
	    END IF
	  END DO
	  IF(XOPT .NE. 'WROPAC')CALL DP_CURVE(NFREQ,ZV,YV)
!
	ELSE IF(XOPT .EQ. 'INTERP')THEN
!
! Routine interpolates the last varible plotted. Variable can only
! be plotted against R. Crude, but effective.
! Routine assumes YV is in the correct form for the interpolation.
!
	  DO I=1,ND
	    TA(I)=YV(I)
	  END DO
!
! Do the interpolation - Based on INTERPTHREE.
!
	  DO I=1,NDX
	    TB(I)=0.0D0
	    DO J=0,3
	      TB(I)=TB(I)+COEF(J,I)*TA(J+INDX(I))
	    END DO
	    YV(I)=TB(I)
	  END DO
!
	  DO I=1,NDX
	    XNU(I)=DLOG10( REXT(I)/REXT(NDX) )
	  END DO
	  CALL DP_CURVE(NDX,XNU,YV)
!
! Option computes a new radius grid, equally spaced in Log(Tau).
!
	ELSE IF(XOPT .EQ. 'NEWRG')THEN
!
! Compute the optical depth scale, which is stored, in LOG form, in TA.
!
	  IF(.NOT. ELEC)THEN
	    DO I=1,ND
	      CHI(I)=CHI(I)-ESEC(I)
	    END DO
	  END IF
	  CALL TORSCL(TA,CHI,R,TB,TC,ND,METHOD,TYPE_ATM)
	  DO I=1,ND
	    TA(I)=DLOG10(TA(I))
	  END DO
!
	  DEFAULT=WR_STRING(ND)
	  CALL USR_OPTION(ND_TMP,'ND=',DEFAULT,'Number of radal depths')
	  IF(ND_TMP .GT. NP_MAX)THEN
	    WRITE(T_OUT,*)'Error in MAINGEN'
	    WRITE(T_OUT,*)'Value for ND is too large'
	    WRITE(T_OUT,*)'Maximum value for ND is',ND
	    GOTO 1
	  END IF
!
! Compute the new optical depth scale, equally spaced in Log(TAU), with
! a finer spacing near both boundaies.
!
	  DTAU(1)=(TA(ND)-TA(1))/(ND_TMP-5)
	  TB(1)=TA(1)
	  TB(2)=TA(1)+DTAU(1)/10
	  TB(3)=TA(1)+DTAU(1)/2
	  DO I=4,ND_TMP-3
	    TB(I)=TA(1)+(I-3)*DTAU(1)
	  END DO
	  TB(ND_TMP)=TA(ND)
	  TB(ND_TMP-1)=TA(ND)-DTAU(1)/10
	  TB(ND_TMP-2)=TA(ND)-DTAU(1)/2
!
	  CALL MON_INTERP(TC,ND_TMP,IONE,TB,ND,R,ND,TA,ND)
	  CALL MON_INTERP(XM,ND_TMP,IONE,TB,ND,V,ND,TA,ND)
	  CALL MON_INTERP(RJ,ND_TMP,IONE,TB,ND,SIGMA,ND,TA,ND)
!
	  WRITE(T_OUT,*)'Plotting Tau, V, Sigma versus Log(r/RP).'
	  WRITE(T_OUT,*)'Plotting Tau, V, Sigma versus depth index.'
	  WRITE(T_OUT,*)'Choose plot in plot package.'

	  T1=TC(ND)
	  DO I=1,ND
	    XNU(I)=LOG10(TC(I)/T1)
	  END DO
	  CALL DP_CURVE(NDX,XNU,TC)
	  CALL DP_CURVE(NDX,XNU,XM)
	  CALL DP_CURVE(NDX,XNU,RJ)
!
	  DO I=1,ND
	    XNU(I)=I
	  END DO
	  CALL DP_CURVE(NDX,XNU,TC)
	  CALL DP_CURVE(NDX,XNU,XM)
	  CALL DP_CURVE(NDX,XNU,RJ)
!
	  CALL GEN_ASCI_OPEN(27,'NEW_R_GRID','UNKNOWN',' ',' ',IZERO,IOS)
	    WRITE(27,*)ND_TMP
	    WRITE(27,'(3X,1P,3E17.8)')(TC(I),XM(I),RJ(I),I=1,ND_TMP)
	  CLOSE(UNIT=27)
! 
!
!
	ELSE IF(XOPT .EQ. 'EW')THEN
!
	  DO I=1,ND
	    SOURCE(I)=ZETA(I)+THETA(I)*RJ(I)
	  END DO
!
! TA is used for the line flux. Integral of TA dlog(r) is
! the line EW.
!
	  FORCE_MULT(1:ND)=0.0D0
	  CALL SOBEW_GRAD_V2(SOURCE,CHI,CHI_SCAT,CHIL,ETAL,
	1             V,SIGMA,R,P,FORCE_MULT,LUM,
	1             JQW,HQW,TA,T1,S1,
	1             FREQ,'HOLLOW',DBB,IC,THICK,.FALSE.,NC,NP,ND,METHOD)
!
          OPEN(UNIT=18,FILE='DSOB_FORCE_MULT',STATUS='UNKNOWN')
	    WRITE(18,'(3X,A1,10X,A1,15X,A1,13X,A1)')'I','R','V','M'
            DO I=1,ND
              WRITE(18,'(1X,I3,3X,3ES14.6)')I, R(I),V(I),FORCE_MULT(I)
            END DO
          CLOSE(UNIT=18)
!
! S1 is the continuum flux in Jy for an object at 1kpc.
! T1 is the line equivalent width in Angstroms.
! T3 is the line flux in ergs/cm^2/s
!
	  T2=LAMVACAIR(FREQ)		!Wavelength(Angstroms)
	  T3=T1*S1*1.0D-23*FREQ*1.0D+15/T2
	  WRITE(LU_NET,40008)T1,S1,T2,T3
	  WRITE(T_OUT,'(A)')RED_PEN
	  WRITE(T_OUT,40008)T1,S1,T2,T3
	  WRITE(T_OUT,'(A)')DEF_PEN
!
	  CALL USR_OPTION(ELEC,'PLOT','F','PLot Line Origin?')
!
! We can now use the current X-axis for the X-axis. Zeta(x)
! is defined similarly to before: Zeta(x) dx gives the
! emission about x in an interval dx.
!
! We use JNU for COEF which must be dimensioned JNU(ND,*) [with * > 3]
! NB: T1 contains the line EW.
!
	  IF(ELEC)THEN
	    WRITE(T_OUT,*)' '
	    WRITE(T_OUT,*)BLUE_PEN//'Area under curve is normalized to unity'//DEF_PEN
	    WRITE(T_OUT,*)' '
	    DO I=1,ND-2
	      IF( (XV(I+1)-XV(I))*(XV(I+2)-XV(I+1)) .LE. 0.0D0)THEN
	        WRITE(T_OUT,*)RED_PEN//'Error -- can only plot line origin against a mononotonix X-axis'//DEF_PEN
	        GOTO 1
	      END IF
	    END DO 
	    CALL MON_INT_FUNS_V2(JNU,XV,R,ND)
	    DO I=1,ND
	      YV(I)=TA(I)/R(I)/JNU(I,3)/T1
	    END DO
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='\gz'
	  ELSE
	    CALL DP_CURVE(ND,XV,FORCE_MULT)
	    YAXIS='Force Multiplier'
	  END IF
!
	ELSE IF(XOPT .EQ. 'EP' .OR. XOPT .EQ. 'BETA')THEN
!
! Used to indicate where the line emission is coming from.
! It should be plotted against log(R) . Area under curve is
! emission. Takes continuous opacity (not .e.s.) into account.
! This section requires tha the line and continuous opacity have
! been previously computed.
!
	  CALL USR_HIDDEN(ELEC,'ELEC','T',
	1      'Ignore Electron Scattering?')
	  IF(ELEC)THEN
	    DO I=1,ND
	      CHI(I)=CHI(I)-ESEC(I)
	    END DO
	  END IF
	  CALL USR_HIDDEN(ELEC,'NEW','T','Continuum included?')
	  CALL USR_HIDDEN(LINV,'LINV','F','Linear R Scale?')
!
! We are using ETA for BETA, the escape probability.
!
	  IF(ELEC)THEN
	    CALL BETANEW(CHI,CHIL,SIGMA,ETA,R,Z,V,FREQ,ND)
	  ELSE
!	    CALL TORSCL(TA,CHI,R,TB,TC,ND,METHOD,TYPE_ATM)
 	    CALL WRBETA(CHIL,SIGMA,ETA,R,V,FREQ,ND)
	  END IF
!
	  IF(XOPT .EQ. 'EP')THEN
	    YAXIS='\gc'
	    IF(LINV)THEN
	      DO I=1,ND
	        ZV(I)=R(I)/R(ND)
	        YV(I)=ETA(I)*R(I)*R(I)*ETAL(I)
	      END DO
	      XAXIS='r/R\d*\u'
	    ELSE
	      DO I=1,ND
	        ZV(I)=DLOG10(R(I)/R(ND))
!	        YV(I)=ETA(I)*R(I)*R(I)*R(I)*ETAL(I)*DEXP(-TA(I))
	        YV(I)=ETA(I)*R(I)*R(I)*R(I)*ETAL(I)
	      END DO
	      XAXIS='Log(r/R\d*\u)'
	    END if
	    T1=0.0
	    DO I=1,ND-1
	      T1=T1+0.5*(YV(I)+YV(I+1))*(ZV(I)-ZV(I+1))
	    END DO
!
! Normalize to unit area.
!
	    CALL USR_HIDDEN(SCALE,'SCL','T','Your scalinge for EP?')
	    IF(.NOT. SCALE)THEN
	      WRITE(T_OUT,*)'SCALE PARAM=',T1
	      CALL USR_OPTION(T1,'T1',' ','Scale parameter')
	    END IF
	    DO I=1,ND
	      YV(I)=YV(I)/T1
	    END DO
!
	  ELSE
	    Yaxis='\gb'
	    XAXIS='Log(r/R\d*\u)'
	    DO I=1,ND
	      ZV(I)=DLOG10(R(I)/R(ND))
	      YV(I)=ETA(I)
	    END DO
	  END IF
	  CALL DP_CURVE(ND,XV,YV)
!
! 
!
! Require CHIL to have been computed in setup.
!
	ELSE IF(XOPT .EQ. 'TAUL')THEN
	  CALL USR_HIDDEN(ELEC,'STAT','F','Stationary opactical depth?')
	  CALL USR_HIDDEN(RADIAL,'RADS','F','Radial Sobolev optical depth?')
	  IF(ELEC)THEN
	    CALL TORSCL(TA,CHIL,R,TB,TC,ND,METHOD,TYPE_ATM)
!
! Assumes V_D=10kms.
!
	    T1=DLOG10(1.6914D-11/FREQ)
	    DO I=1,ND
	      IF(TA(I) .GT. 0)THEN
	        YV(I)=T1+DLOG10(TA(I))
	      ELSE IF(TA(I) .LT. 0)THEN
	        YV(I)=T1+DLOG10(-TA(I))-20.0D0
	      ELSE
	        YV(I)=-30.0
	      END IF
	    END DO
	    YAXIS='Log(\gt\dstat\u)'
	  ELSE
	    DO I=1,ND
	      YV(I)=CHIL(I)*R(I)*2.998E-10/FREQ/V(I)
	      IF(RADIAL)YV(I)=YV(I)/(1.0D0+SIGMA(I))
	      IF(YV(I) .GT. 0)THEN
	        YV(I)=LOG10(YV(I))
	      ELSE
	        YV(I)=-20.0
	      END IF
	    END DO
	    YAXIS='Log(\gt\dSob\u)'
	  END IF
	  CALL DP_CURVE(ND,XV,YV)
!
! Require CHIL to have been computed in setup.
!
	ELSE IF(XOPT .EQ. 'CHIL')THEN
!
! Assumes V_D=10kms.
!
	  CALL USR_OPTION(ELEC,'KAPPA','T','Plot kappa insted of chi?')
	  WRITE(6,*)'A Doppler velocity of 10 km/s is assumed'
	  WRITE(6,*)'Warning: CHIL now has units of cm^-1 - not program units'
	  IF(ELEC)THEN
	    T1=1.6914D-21/FREQ
	    YV(1:ND)=T1*CHIL(1:ND)/MASS_DENSITY(1:ND)/CLUMP_FAC(1:ND)
	    YAXIS='\gk(cm\u3 \d/g)'
	    CALL DP_CURVE(ND,XV,YV)
	  ELSE
	    VALID_VALUE=.TRUE.
	    ZV(1:ND)=0.0D0; YV(1:ND)=0.0D0
	    T1=DLOG10(1.6914D-11/FREQ)-10.0D0
	    DO I=1,ND
	      IF(CHIL(I) .GT. 0)THEN
	        YV(I)=T1+DLOG10(CHIL(I))
	      ELSE IF(CHIL(I) .LT. 0)THEN
	        ZV(I)=T1+DLOG10(-CHIL(I))
	        VALID_VALUE=.FALSE.
	      ELSE
	        YV(I)=-30.0
	      END IF
	    END DO
	    YAXIS='Log \gx\dL\u(cm\u-1\d)'
	    CALL DP_CURVE(ND,XV,YV)
	    IF(.NOT. VALID_VALUE)THEN
	      WRITE(6,*)'Warning -- negative line opacities'
	      WRITE(6,*)'Sending to data plots to DURVE - use MARK to illustrate'
	      CALL DP_CURVE(ND,XV,ZV)
	    END IF
	  END IF
! 
!
! Write out line/continuum opacities and emissivities for use with
! the polarization codes, or profile codes.
!
	ELSE IF(XOPT .EQ. 'WRC' .OR. XOPT .EQ. 'WRL')THEN
	  CONT_INT=TWOHCSQ*(FREQ**3)/( EXP(HDKT*FREQ/T(ND))-1.0D0 )
	  CALL USR_OPTION(NEW_FORMAT,'NEW_FORM','F',
	1        'Output multiple lines in new format file?')
	  IF(NEW_FORMAT)THEN
	    CALL USR_OPTION(NEW_FILE,'NEW_FILE','F','Open new file?')
	  END IF
!
	  IF(NEW_FORMAT .AND. NEW_FILE)THEN
	    OPEN(UNIT=25,FILE='LINEDATA',STATUS='NEW')
   	    WRITE(25,'(1X,A,T25,A,T40,A)')'12-May-1998',
	1            '[Date]','Revised format date'
   	    WRITE(25,'(1X,A,T25,A,T40,A)')TRIM(NAME),
	1             '[MOD_ID]','Model'
 	    WRITE(25,'(1X,A,T25,A,T40,A)')'TRUE',
	1            '[DIF]','Diffusion approximation'
	    WRITE(25,'(1X,1PE15.8,T25,A,T40,A)')CONT_INT,
	1            '[IC]','Schuster intensity'
	    WRITE(25,*)' '
!
	    WRITE(25,'(1X,A,T25,A,T40,A)')'Continuum:',
	1               '[TR_ID]','Transition identification'
	    WRITE(25,'(1X,1PE15.8,T25,A,T40,A)')FREQ,
	1            '[FREQ]','Frequency (10^15 Hz)'
	    WRITE(25,'(1X,1PE12.5,T25,A,T40,A)')LAMVACAIR(FREQ),
	1            '[LAM]','Wavelength (Ang)'
!
	    WRITE(25,*)' '
	    WRITE(25,'(1X,A,/,(1P,1X,9E14.6))')'R',(R(I),I=1,ND)
	    WRITE(25,'(1X,A,/,(1P,1X,9E14.6))')'T',(T(I),I=1,ND)
	    WRITE(25,'(1X,A,/,(1P,1X,9E14.6))')'SIGMA',(SIGMA(I),I=1,ND)
	    WRITE(25,'(1X,A,/,(1P,1X,9E14.6))')'V',(V(I),I=1,ND)
	    WRITE(25,'(1X,A,/,(1P,1X,9E14.6))')'ETA',(ETA(I),I=1,ND)
	    WRITE(25,'(1X,A,/,(1P,1X,9E14.6))')'CHI_TH',
	1                              ( (CHI(I)-ESEC(I)) ,I=1,ND)
	    WRITE(25,'(1X,A,/,(1P,1X,9E14.6))')'ESEC',(ESEC(I),I=1,ND)
	  END IF
!
	  IF(NEW_FORMAT .AND. XOPT .EQ. 'WRL')THEN
	    IF(.NOT. NEW_FILE)THEN
	      OPEN(UNIT=25,FILE='LINEDATA',STATUS='OLD',IOSTAT=IOS)
	      IF(IOS .NE. 0)THEN
	         WRITE(T_OUT,*)'Error - unable to open old LINEDATA file'
	      END IF
	    END IF
	    WRITE(25,*)' '
	    WRITE(25,'(1X,A,A,I3,A,I3,A,T25,A,T40,A)')
	1            TRIM(XSPEC),'(',LEV(2),'-',LEV(1),')',
	1            '[TR_ID]','Transition identification'
	    WRITE(25,'(1X,1PE15.8,T25,A,T40,A)')FREQ,
	1            '[FREQ]','Frequency (10^15 Hz)'
	    WRITE(25,'(1X,1PE12.5,T25,A,T40,A)')LAMVACAIR(FREQ),
	1            '[LAM]','Wavelength (Ang)'
	    WRITE(25,'(1X,F6.2,T25,A,T40,A)')AMASS,'[AMASS]','Atomic mass'
	    WRITE(25,*)' '
	    WRITE(25,'(1X,A,/,(1P,1X,9E14.6))')'ETAL',(ETAL(I),I=1,ND)
	    WRITE(25,'(1X,A,/,(1P,1X,9E14.6))')'CHIL',(CHIL(I),I=1,ND)
	  END IF

	  IF(.NOT. NEW_FORMAT)THEN
	    OPEN(UNIT=25,FILE='LINEDATA',STATUS='NEW',IOSTAT=IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(T_OUT,*)'Error opening LINEDATA file: probably because file exists'
	      DEFAULT='LINEDATA_1'
	      CALL USR_OPTION(FILENAME,'FILE',DEFAULT,
	1                  'File name (existing file will be overwrittem')
	      OPEN(UNIT=25,FILE=FILENAME,STATUS='UNKNOWN')
	    END IF
   	    WRITE(25,'(1X,A,T25,A,T40,A)')'09-Mar-1994',
	1            '[Date]','Revised format date'
   	    WRITE(25,'(1X,A,T25,A,T40,A)')TRIM(NAME),
	1            '[MOD_ID]','Model'
	    IF(XOPT .EQ. 'WRL')THEN
	        WRITE(25,'(1X,A,A,I3,A,I3,A,T25,A,T40,A)')
	1            TRIM(XSPEC),'(',LEV(2),'-',LEV(1),')',
	1            '[TR_ID]','Transition identification'
	      ELSE
	        WRITE(25,'(1X,A,T25,A,T40,A)')'Continuum:',
	1               '[TR_ID]','Transition identification'
	      END IF
 	      WRITE(25,'(1X,A,T25,A,T40,A)')'TRUE',
	1            '[DIF]','Diffusion approximation'
	      WRITE(25,'(1X,1PE15.8,T25,A,T40,A)')FREQ,
	1            '[FREQ]','Frequency (10^15 Hz)'
	      WRITE(25,'(1X,1PE12.5,T25,A,T40,A)')LAMVACAIR(FREQ),
	1            '[LAM]','Wavelength (Ang)'
	      WRITE(25,'(1X,F6.2,T25,A,T40,A)')AMASS,
	1            '[AMASS]','Atomic mass'
	      WRITE(25,'(1X,1PE15.8,T25,A,T40,A)')CONT_INT,
	1            '[IC]','Schuster intensity'
	      WRITE(25,'(1X,I3,T25,A,T40,A)')ND,
	1            '[ND]','Number of depth points'
	      WRITE(25,*)' '
	      WRITE(25,
	1       '(1X,T9,A,T23,A,T35,A,T51,A,T64,A,T76,A,T91,A,T105,A,T119,A)')
	1            'R','T','SIGMA','V','ETA','CHI_TH','ESEC',
	1            'ETAL','CHIL'
	      IF(XOPT .EQ. 'WRC')THEN
	        DO I=1,ND
	          CHIL(I)=1.0D-10
	          ETAL(I)=1.0D-10
	        END DO
	      END IF
	      DO I=1,ND
 	        WRITE(25,555)R(I),T(I),SIGMA(I),V(I),
 	1                 ETA(I),(CHI(I)-ESEC(I)),
	1                 ESEC(I),ETAL(I),CHIL(I)
	      END DO
	    CLOSE(UNIT=25)
555	    FORMAT(1X,1P9E14.6)
	  END IF
!
! Option to output nformation concerning bound-bound transitions.
! Output file has same format as that generated by CMFGEN. Only a section
! of wavelength space needs to be output.
!
	ELSE IF(XOPT .EQ. 'WRTRANS')THEN
!
	  DEFAULT=WR_STRING(LAM_ST)
	  CALL USR_OPTION(LAM_ST,'LAMST',DEFAULT,FREQ_INPUT)
	  LAM_EN=LAM_ST*1.005D0
	  DEFAULT=WR_STRING(LAM_EN)
	  CALL USR_OPTION(LAM_EN,'LAMEN',DEFAULT,FREQ_INPUT)
	  IF(KEV_INPUT)THEN
	    NU_ST=LAM_ST*KEV_TO_HZ
	    NU_EN=LAM_EN*KEV_TO_HZ
	  ELSE IF(ANG_INPUT)THEN
	    NU_ST=ANG_TO_HZ/LAM_ST
	    NU_EN=ANG_TO_HZ/LAM_EN
	  ELSE
	    NU_ST=LAM_ST
	    NU_EN=LAM_EN
	  END IF
!
	  NL=GET_INDX_DP(NU_ST,VEC_FREQ,N_LINE_FREQ)
	  NUP=GET_INDX_DP(NU_EN,VEC_FREQ,N_LINE_FREQ)
!	  
	  I=160	!Record length - allow for long names
	  CALL GEN_ASCI_OPEN(LU_OUT,'TRANS_INFO','UNKNOWN',' ','WRITE',I,IOS)
	  WRITE(LU_OUT,*)
	1     '     I  NL_F  NUP_F      Nu',
	1     '       Lam(A)      gf      /\V(km/s)    Transition' 
	  IF(NUP-NL .LT. 50)WRITE(T_OUT,*)
	1     '     I  NL_F  NUP_F      Nu',
	1     '       Lam(A)      gf      /\V(km/s)    Transition' 
	  DO ML=NL,NUP
	    T1=LAMVACAIR(VEC_FREQ(ML))
	    T2=2.998D+05*(VEC_FREQ(MAX(1,ML-1))-VEC_FREQ(ML))/VEC_FREQ(ML)
	    IF(T2 .GT. 2.998E+05)T2=2.998E+05
	    ID=VEC_ION_INDX(ML)
	    T3=VEC_OSCIL(ML)*ATM(ID)%GXzV_F(VEC_MNL_F(ML))	!gf
	    IF(T1 .LT. 1.0E+04)THEN
	      WRITE(LU_OUT,
	1      '(1X,I6,2I6,F10.6,2X,F10.3,ES10.2,1X,F10.2,4X,A)')
	1         ML,VEC_MNL_F(ML),VEC_MNUP_F(ML),
	1         VEC_FREQ(ML),T1,T3,T2,TRIM(VEC_TRANS_NAME(ML))
	      IF(NUP-NL .LT. 50)WRITE(T_OUT,
	1      '(1X,I6,2I6,F10.6,2X,F10.3,ES10.2,1X,F10.2,4X,A)')
	1         ML,VEC_MNL_F(ML),VEC_MNUP_F(ML),
	1         VEC_FREQ(ML),T1,T3,T2,TRIM(VEC_TRANS_NAME(ML))
	    ELSE             
	      WRITE(LU_OUT,
	1      '(1X,I6,2(1X,I6),F10.6,1X,ES11.4,ES10.2,1X,F10.2,4X,A)')
	1         ML,VEC_MNL_F(ML),VEC_MNUP_F(ML),
	1         VEC_FREQ(ML),T1,T3,T2,TRIM(VEC_TRANS_NAME(ML))
	      IF(NUP-NL .LT. 50)WRITE(T_OUT,
	1      '(1X,I6,2(1X,I6),F10.6,1X,ES11.4,ES10.2,1X,F10.2,4X,A)')
	1         ML,VEC_MNL_F(ML),VEC_MNUP_F(ML),
	1         VEC_FREQ(ML),T1,T3,T2,TRIM(VEC_TRANS_NAME(ML))
	    END IF
	  END DO
	CLOSE(UNIT=LU_OUT)

	ELSE IF(XOPT .EQ. 'WRRTK')THEN
	  DO I=1,ND
	    WRITE(25,'(I3,ES15.5,3ES14.4)')I,R(I)*1.0D+10,T(I)*1.0D+04,
	1          MASS_DENSITY(I),1.0D-10*ROSS_MEAN(I)/MASS_DENSITY(I)
	  END DO
	ELSE IF(XOPT .EQ. 'WRRVSIG')THEN
	  CALL USR_OPTION(ELEC,'DCF','F','Include the density and clumping factor in file?')
	  CALL GEN_ASCI_OPEN(LU_OUT,'NEW_RVSIG','UNKNOWN',' ','WRITE',I,IOS)
	  WRITE(LU_OUT,'(A)')'!'
	  IF(ELEC)THEN
	    WRITE(LU_OUT,'(A,7X,A,9X,10X,A,12X,A,2(7X,A),3X,A)')'!','R','V(km/s)','Sigma',
	1                   '    Density','Clump. Fac.','Depth'
	  ELSE
	    WRITE(LU_OUT,'(A,7X,A,9X,10X,A,12X,A,3X,A)')'!','R','V(km/s)','Sigma','Depth'
	  END IF
	  WRITE(LU_OUT,'(A)')'!'
	  WRITE(LU_OUT,'(A)')' '
	  WRITE(LU_OUT,'(I4,20X,A)')ND,'!Number of depth points`'
	  WRITE(LU_OUT,'(A)')' '
	  IF(ELEC)THEN
	    DO I=1,ND
	      WRITE(LU_OUT,'(F18.8,ES17.7,F17.7,2ES18.8,4X,I4)')R(I),V(I),SIGMA(I),MASS_DENSITY(I),CLUMP_FAC(I),I
	    END DO
	  ELSE
	    DO I=1,ND
	      WRITE(LU_OUT,'(F18.8,ES17.7,F17.7,4X,I4)')R(I),V(I),SIGMA(I),I
	    END DO
	  END IF
	  CLOSE(LU_OUT)
	  WRITE(T_OUT,'(A,ES17.10)')'     R(ND)=',R(ND)
	  WRITE(T_OUT,'(A,ES17.10)')'RMAX/R(ND)=',R(1)/R(ND)
	  WRITE(T_OUT,'(A)')'RVSIG data output to NEW_RVSIG'
	ELSE IF(XOPT .EQ. 'WRPLOT')THEN
	  WRITE(T_OUT,*)' OPTION `WRPLOT` NOT AVAILABLE'
	  WRITE(T_OUT,*)' Use WP option in plot package'
	ELSE IF(XOPT .EQ. 'RDPLOT')THEN
	  WRITE(T_OUT,*)' OPTION `RDPLOT` NOT AVAILABLE'
	  WRITE(T_OUT,*)' Use RP option in plot package'
!
! 
!
! Output suumary of commands, or a brief command summary.
! The files containing the descriptions need to be linked via DEFINE
! statements (VMS), or through soft links (UNIX).
!
	ELSE IF(XOPT .EQ. 'LI' .OR. XOPT .EQ. 'LIST' .OR.
	1       XOPT .EQ. 'HE' .OR. XOPT .EQ. 'HELP')THEN
	  IF(XOPT(1:2) .EQ. 'LI')THEN
	    CALL GEN_ASCI_OPEN(LU_IN,'MAINGEN_OPT_DESC','OLD',' ','READ',
	1                  IZERO,IOS)
	  ELSE 
	    CALL GEN_ASCI_OPEN(LU_IN,'MAINGEN_OPTIONS','OLD',' ','READ',
	1                  IZERO,IOS)
	  END IF
!
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error opening description file'
	    GOTO 1
	  END IF
	  READ(LU_IN,*)I,K			!For page formating (I=22,K=12)
	  DO WHILE(1.EQ. 1)
	    DO J=1,I
	      READ(LU_IN,'(A)',END=700)STRING
	      L=LEN_TRIM(STRING)
	      IF(L .EQ. 0)THEN
	        WRITE(T_OUT,'(1X)')
	      ELSE
	        WRITE(T_OUT,'(1X,A)')STRING(1:L)
	      END IF
	    END DO
	    READ(T_IN,'(A)')STRING
	    IF(STRING(1:1) .EQ. 'E' .OR. STRING(1:1) .EQ. 'e' .OR.
	1       STRING(1:1) .EQ. 'Q' .OR.
	1       STRING(1:1) .EQ. 'q')GOTO 700		!Exit from listing.
	    I=K
	  END DO
700	  CONTINUE
	  CLOSE(UNIT=LU_IN)
!
	ELSE IF(XOPT .EQ. 'EX' .OR. XOPT .EQ. 'EXIT') THEN
	  STOP
	ELSE IF(X(1:4) .EQ. 'BOX=') THEN
	  CALL WR_BOX_FILE(MAIN_OPT_STR)
	ELSE
	  PRINT*,'OPTION REQUESTED DOES NOT EXIST'
	END IF
!
 1	  GO TO 3
!
! Formats used for writing out headers, net rates tec.
!
40001	  FORMAT(/,1X,A70)
40002	  FORMAT(1X,'NL=',I3,5X,'NUP=',I3)
40003	  FORMAT(3X,1P5E16.5)
!
	  END
