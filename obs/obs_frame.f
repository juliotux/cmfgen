C
C Calling routine to compute the Observer's fluxes. The fluxes are computed
C using an OBSERVER'S FRAME formulation. J, ETA, and CHI, computed in the
C comoving frame, must be supplied.
C
	PROGRAM OBS_FRAME
	IMPLICIT NONE
C
C Altered 06-Jan-1999 --- Now call OBS_FRAME_SUB_V2.
C                         TAU_MAX, ES_DTAU, and INT_METHOD included.
C Altered 03-Jan-1999 --- MAX_DEL_V_RES_ZONE was being incorrectly computed.
C Altered 11-Feb-1998 --- RJ_VEC was being set outside valid range.
C Created 9-Dec-1998 :
C
	INTEGER ND,NP,NC

	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
	INTEGER, PARAMETER :: IFOUR=4
	INTEGER, PARAMETER :: ISIX=6
	INTEGER, PARAMETER :: ITEN=10
	INTEGER, PARAMETER :: T_OUT=6		!Terminal IO
C
	INTEGER, PARAMETER :: LUIN=8
	INTEGER, PARAMETER :: LUMOD=9
	INTEGER, PARAMETER :: LU_FLUX=20
C
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
C
	REAL*8, ALLOCATABLE :: R(:)
	REAL*8, ALLOCATABLE :: V(:)
	REAL*8, ALLOCATABLE :: SIGMA(:)
	REAL*8, ALLOCATABLE :: T(:)
	REAL*8, ALLOCATABLE :: ED(:)
	REAL*8, ALLOCATABLE :: ESEC(:)
C
C These vectors are neeed so that we can use the standard RVTJ read routine.
C Only CLUMP_FAC is utilized. 
C
	REAL*8, ALLOCATABLE :: ROSS_MEAN(:)
	REAL*8, ALLOCATABLE :: FLUX_MEAN(:)
	REAL*8, ALLOCATABLE :: POP_ATOM(:)
	REAL*8, ALLOCATABLE :: MASS_DENSITY(:)
	REAL*8, ALLOCATABLE :: CLUMP_FAC(:)
	REAL*8, ALLOCATABLE :: POPION(:)
C
	REAL*8, ALLOCATABLE :: VTURB(:)
	REAL*8, ALLOCATABLE :: MAX_DEL_V_RES_ZONE(:)
C
C P is the impact parameter. PQW is R^2 p dp
C
	REAL*8, ALLOCATABLE :: P(:)
	REAL*8, ALLOCATABLE :: MU(:)
	REAL*8, ALLOCATABLE :: HQW_AT_RMAX(:)
C
	REAL*8 RMDOT
	REAL*8 RLUM
	REAL*8 ABUND_HYD
	CHARACTER*20 FILE_DATE
	CHARACTER*20 TIME
	CHARACTER*10 NAME_CONVENTION
C
C 
	INTEGER NLF
	REAL*8, ALLOCATABLE :: FREQ_CMF(:)
	REAL*8, ALLOCATABLE :: CHI_CMF(:,:)
	REAL*8, ALLOCATABLE :: ETA_CMF(:,:)
	REAL*8, ALLOCATABLE :: NEW_RJ_CMF(:,:)
C
	REAL*8, ALLOCATABLE :: FREQ_RJ(:)
	REAL*8, ALLOCATABLE :: RJ_CMF(:,:)
	REAL*8, ALLOCATABLE :: RJ_VEC(:)
	REAL*8, ALLOCATABLE :: NEW_RJ_VEC(:)
C
	INTEGER NOS
	REAL*8, ALLOCATABLE :: OBS_FREQ(:)
	REAL*8, ALLOCATABLE :: OBS_FLUX(:)
C 
C
	EXTERNAL ERROR_LU,SPEED_OF_LIGHT
	INTEGER ERROR_LU,LUER
	REAL*8 SPEED_OF_LIGHT,C_KMS
C
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
C
C
	REAL*8 TAU_MAX
	REAL*8 ES_DTAU
	REAL*8 FRAC_DOP
	REAL*8 DEL_V_OBS
	REAL*8 VTURB_MIN
	REAL*8 VTURB_MAX
	REAL*8 MIN_FREQ
	REAL*8 MAX_FREQ
	REAL*8 MIN_WAVE
	REAL*8 MAX_WAVE
	REAL*8 T1
	CHARACTER*10 INT_METHOD
C
	INTEGER NCF
	INTEGER NCF_RJ
	INTEGER NOBS
	INTEGER ND_RD
	INTEGER INIT_REC
	INTEGER IOS
	INTEGER I,K,LS,ML
	INTEGER CONT_REC
	INTEGER RECL
C
	LOGICAL INTERP_RJ_NEC
	LOGICAL CONVOLVE_J
	LOGICAL DO_CLUMP
C 
C
C Set constants.
C
	CHIBF=2.815E-06
	CHIFF=3.69E-29
	HDKT=4.7994145
	TWOHCSQ=0.0147452575
	OPLIN=2.6540081E+08
	EMLIN=5.27296E-03
	OPLIN=2.6540081E+08
	EMLIN=5.27296E-03
C
	CONT_REC=3
	LUER=ERROR_LU()
	C_KMS=1.0D-05*SPEED_OF_LIGHT()
C
	CALL GEN_ASCI_OPEN(LUIN,'OBS_INP','OLD',' ','READ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening OBS_INP in CMFGEN, IOS=',IOS
	  STOP
	END IF
	CALL RD_DBLE(MIN_WAVE,'MIN_WAVE',LUIN,LUMOD,
	1       'Minimum observers wavelenghth (Ang)')
	CALL RD_DBLE(MAX_WAVE,'MAX_WAVE',LUIN,LUMOD,
	1       'Maximum observers wavelenghth (Ang)')
	CALL RD_DBLE(DEL_V_OBS,'DEL_V_OBS',LUIN,LUMOD,
	1       'Spacing (km/s) for observers frame frequencies ')
C
C NB: VTURB_MIN and VTURB_MAX do NOT directly effect the spectrum. They
C only effect the accuracy of its computation.
C
	CALL RD_DBLE(VTURB_MIN,'VTURB_MIN',LUIN,LUMOD,
	1       'Minimum turbulent velocity (km/s)')
	CALL RD_DBLE(VTURB_MAX,'VTURB_MAX',LUIN,LUMOD,
	1       'Maximum turbulent velocity (units of Vinf)')
	CALL RD_DBLE(FRAC_DOP,'FRAC_DOP',LUIN,LUMOD,
	1       'Fractional spacing (Doppler widths) across resonance one')
C
	CALL RD_LOG(CONVOLVE_J,'CONV_ES_J',LUIN,LUMOD,
	1       'Convolve J with e.s. redistribution function?')
C
	CALL RD_DBLE(TAU_MAX,'TAU_MAX',LUIN,LUMOD,
	1       'Maximum TAU along ray at which integration ceases.')
	CALL RD_DBLE(ES_DTAU,'ES_DTAU',LUIN,LUMOD,
	1       'Maximum DTAU on Electron scattering optical depth scale.')
	CALL RD_NCHAR(INT_METHOD,'INT_METH',ITEN,LUIN,LUMOD,
	1      'Intgeration method [STAU or ETAZ]')
C
	CLOSE(LUIN)
C
C Read in atmosphere parameters.
C
	CALL RD_RVTJ_PARAMS_V2(RMDOT,RLUM,ABUND_HYD,TIME,NAME_CONVENTION,
	1                             ND,NC,NP,'RVTJ',LUIN)
	ALLOCATE (R(ND))
	ALLOCATE (V(ND))
	ALLOCATE (SIGMA(ND))
	ALLOCATE (T(ND))
	ALLOCATE (ED(ND))
	ALLOCATE (ROSS_MEAN(ND))
	ALLOCATE (FLUX_MEAN(ND))
	ALLOCATE (POP_ATOM(ND))
	ALLOCATE (MASS_DENSITY(ND))
	ALLOCATE (POPION(ND))
	ALLOCATE (CLUMP_FAC(ND))
	CALL RD_RVTJ_VEC(R,V,SIGMA,ED,T,ROSS_MEAN,FLUX_MEAN,
	1       POP_ATOM,POPION,MASS_DENSITY,CLUMP_FAC,ND,LUIN)
	CLOSE(LUIN)
C
	CALL READ_DIRECT_INFO_V3(I,RECL,FILE_DATE,'ETA_DATA',LUIN,IOS)
	IF(IOS .NE. 0)STOP
	OPEN(UNIT=LUIN,FILE='ETA_DATA',STATUS='OLD',
	1      ACCESS='DIRECT',FORM='UNFORMATTED',ACTION='READ')
	READ(LUIN,REC=CONT_REC)INIT_REC,NCF,ND_RD
	IF(ND_RD .NE. ND)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error in OBS_FRAME'
	  WRITE(I,*)'ND and ND_RD must agree'
	  WRITE(I,*)'ND=',ND,'ND_RD=',ND_RD
	  STOP
	END IF
	ALLOCATE (FREQ_CMF(NCF))
	ALLOCATE (ETA_CMF(ND,NCF))
	DO ML=1,NCF
	  READ(LUIN,REC=INIT_REC+ML-1)(ETA_CMF(K,ML),K=1,ND),FREQ_CMF(ML)
	END DO
C
	CALL READ_DIRECT_INFO_V3(I,RECL,FILE_DATE,'CHI_DATA',LUIN,IOS)
	IF(IOS .NE. 0)STOP
	OPEN(UNIT=LUIN,FILE='CHI_DATA',STATUS='OLD',
	1      ACCESS='DIRECT',FORM='UNFORMATTED',ACTION='READ')
	READ(LUIN,REC=CONT_REC)INIT_REC,NCF,ND_RD
	IF(ND_RD .NE. ND)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error in OBS_FRAME'
	  WRITE(I,*)'ND and ND_RD must agree'
	  WRITE(I,*)'ND=',ND,'ND_RD=',ND_RD
	  STOP
	END IF
	ALLOCATE (CHI_CMF(ND,NCF))
	DO ML=1,NCF
	  READ(LUIN,REC=INIT_REC+ML-1)(CHI_CMF(K,ML),K=1,ND),T1
	  IF(T1 .NE. FREQ_CMF(ML))THEN
	    I=ERROR_LU()
	    WRITE(I,*)'ERROR in OBS_FRAME'
	    WRITE(I,*)'ETA and CHI frequencies are not identical'
	    STOP
	  END IF
	END DO
C 
C
C ***************************************************************************
C
C Read in mean intensity as a function of depth and frequency. If requested,
C J will be convolved with the electron scattering redistribution function.
C If the frequency grid for J differs from that for ETA and CHI, it is
C interpolated onto the ETA (and CHI) frequency grid.
C
C ***************************************************************************
C
	CALL READ_DIRECT_INFO(I,RECL,'RJ_DATA',LUIN)
	OPEN(UNIT=LUIN,FILE='RJ_DATA',STATUS='OLD',
	1      ACCESS='DIRECT',FORM='UNFORMATTED',ACTION='READ')
	READ(LUIN,REC=CONT_REC)INIT_REC,NCF_RJ,ND_RD
	IF(ND_RD .NE. ND)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error in OBS_FRAME'
	  WRITE(I,*)'ND and ND_RD must agree'
	  WRITE(I,*)'ND=',ND,'ND_RD=',ND_RD
	  STOP
	END IF
	ALLOCATE (RJ_CMF(ND,NCF_RJ))
	ALLOCATE (FREQ_RJ(NCF_RJ))
	DO ML=1,NCF_RJ
	  READ(LUIN,REC=INIT_REC+ML-1)(RJ_CMF(K,ML),K=1,ND),FREQ_RJ(ML)
	END DO
C
C If requested, convolve J with the electron-scattering redistribution
C funtion. K is used for LUIN, and LUOUIT but is not accessed.
C
	IF(CONVOLVE_J)THEN
	  I=ND*NCF_RJ
	  CALL COMP_J_CONV_V2(RJ_CMF,I,FREQ_RJ,T,ND,NCF_RJ,
	1         K,'J PASSED VIA CALL',CONT_REC,L_FALSE,L_FALSE,
	1         K,'RETURN J VIA CALL')
	END IF
C
	ALLOCATE (ESEC(ND))
	ESEC(1:ND)=6.65D-15*ED(1:ND)
C
C Detemine if RJ_CMF is on same grid as ETA and CHI. If not put on same grid
C via interpolation. We then combine RJ with ETA.
C
	IF(NCF .EQ. NCF_RJ)THEN
	  INTERP_RJ_NEC=.FALSE.
	  DO ML=1,NCF
	    IF(FREQ_CMF(ML) .NE. FREQ_RJ(ML))INTERP_RJ_NEC=.TRUE.
	  END DO
	ELSE
	  INTERP_RJ_NEC=.TRUE.
	END IF
	IF(INTERP_RJ_NEC)THEN
	  ALLOCATE (RJ_VEC(NCF_RJ))
	  ALLOCATE (NEW_RJ_VEC(NCF))
	  ALLOCATE (NEW_RJ_CMF(ND,NCF))
	  DO K=1,ND
	    RJ_VEC(1:NCF_RJ)=RJ_CMF(K,1:NCF_RJ)
	    CALL MON_INTERP(NEW_RJ_VEC,NCF,IONE,FREQ_CMF,NCF,
	1                             RJ_VEC,NCF_RJ,FREQ_RJ,NCF_RJ)
	    NEW_RJ_CMF(K,1:NCF)=NEW_RJ_VEC(1:NCF)
	  END DO
	  DO ML=1,NCF
	    DO I=1,ND
	      ETA_CMF(I,ML)=ETA_CMF(I,ML)+NEW_RJ_CMF(I,ML)*ESEC(I)
	    END DO
	  END DO
	  DEALLOCATE (RJ_CMF)
	  DEALLOCATE (NEW_RJ_CMF)
	  DEALLOCATE (RJ_VEC)
	  DEALLOCATE (NEW_RJ_VEC)
	ELSE
	  DO ML=1,NCF
	    DO I=1,ND
	      ETA_CMF(I,ML)=ETA_CMF(I,ML)+RJ_CMF(I,ML)*ESEC(I)
	    END DO
	  END DO
	  DEALLOCATE (RJ_CMF)
	END IF
C
C 
C ***************************************************************************
C ***************************************************************************
C
C Determine the turbulent velocity as a funtion of depth. This is only
C used to define the points inserted along each ray via MAX_DEL_V_RES_ZONE.
C
	ALLOCATE (VTURB(ND))
	ALLOCATE (MAX_DEL_V_RES_ZONE(ND))
C
	DO I=1,ND
	  VTURB(I)=VTURB_MIN+(VTURB_MAX-VTURB_MIN)*V(I)/V(1)
	END DO
C
	DO I=1,ND
	  MAX_DEL_V_RES_ZONE(I)=VTURB(I)*FRAC_DOP
	END DO
C
C If CLUMPING is important, need to correct emissivities and opacities.
C NB: ETA and CHI read from the files do NOT contain the clumping factor.
C
	DO_CLUMP=.FALSE.
	DO I=1,ND
	  IF( ABS(CLUMP_FAC(I)-1.0D0) .GT. 1.0D-06)DO_CLUMP=.TRUE.
	END DO
C
	IF(DO_CLUMP)THEN
	  DO ML=1,NCF
	    ETA_CMF(:,ML)=ETA_CMF(:,ML)*CLUMP_FAC(:)
	    CHI_CMF(:,ML)=CHI_CMF(:,ML)*CLUMP_FAC(:)
	  END DO
	END IF
C
C 
C ***************************************************************************
C 
C Define the impact parameters and appropriate quarature weights for obtaing
C the observed flux. At present the non-core rays should be defined by the 
C radius grid.
C
C ***************************************************************************
C
	NC=5+ND/4
	NP=ND+NC
	ALLOCATE (P(NP))
	ALLOCATE (MU(NP))
	ALLOCATE (HQW_AT_RMAX(NP))
C
C Compute impact parameter valiues.
C
	CALL IMPAR(P,R,R(ND),NC,ND,NP)
	T1=R(ND)/(NC-1)
	DO I=2,NC-1
	  P(I)=P(I-1)+T1
	END DO
	P(NC)=P(NC+1)-0.001*(P(NC+1)-P(NC))
C
C By definition, p * dp equals R**2 * mu * dmu. Integration over mu is
C more stable, and is to be preferred. To get better accuracy with the
C integration, NORDULUND weights will be used (Changed 11-Dec-1986).
C
	T1=R(1)*R(1)
	DO LS=1,NP
	  MU(LS)=SQRT(T1-P(LS)*P(LS))/R(1)
	END DO
	CALL HWEIGHT(MU,HQW_AT_RMAX,NP)
C
C ***************************************************************************
C ***************************************************************************
C
C Define the observer's frame frequency grid. This grid is in units of
C 10^15 Hz and must be ordered from highest to lowest (consistent with
C the comoving frame frequencies).
C
C ***************************************************************************
C ***************************************************************************
C
C 0.01D0=(1.0D+3/1.0D-10)/1.0D+15
C
	MAX_FREQ=0.01D0*C_KMS/MIN_WAVE
	MIN_FREQ=0.01D0*C_KMS/MAX_WAVE
	NOS=DLOG(MAX_FREQ/MIN_FREQ)/DLOG(1.0D0+DEL_V_OBS/C_KMS)-1
C
	ALLOCATE (OBS_FLUX(NOS))
	ALLOCATE (OBS_FREQ(NOS))
C
	OBS_FREQ(1)=MAX_FREQ
	DO I=2,NOS
	  OBS_FREQ(I)=MAX_FREQ/(1.0+DEL_V_OBS/C_KMS)**(I-1)
	END DO
C
	WRITE(LUER,*)' '
	WRITE(LUER,'(A,I3)')' Number of depth points is: ',ND
	WRITE(LUER,'(A,I6)')' Number of CMF frequency points is: ',NCF
	WRITE(LUER,'(A,I5)')' Number of observer''s frequency points is: ',NOS
	WRITE(LUER,*)' '
C
C ***************************************************************************
C ***************************************************************************
C
C Compute the observer's frame fluxes. The fluxes are returned in Janskies.
C
	CALL OBS_FRAME_SUB_V2(ETA_CMF,CHI_CMF,FREQ_CMF,
	1            R,V,T,ED,ND,NCF,
	1            P,HQW_AT_RMAX,NC,NP,
	1            OBS_FREQ,OBS_FLUX,NOS,
	1            MAX_DEL_V_RES_ZONE,
	1            TAU_MAX,ES_DTAU,INT_METHOD)
	CALL TUNE(3,' ')
C
	CALL GEN_ASCI_OPEN(LU_FLUX,'OBSFRAME','UNKNOWN',' ',' ',IZERO,IOS)
	  CALL WRITV_V2(OBS_FREQ,NOS,ISIX,'Continuum Frequencies',LU_FLUX)
	  CALL WRITV_V2(OBS_FLUX,NOS,IFOUR,
	1                            'Observed intensity (Janskys)',LU_FLUX)
	CLOSE(UNIT=LU_FLUX)
C
	STOP
	END
