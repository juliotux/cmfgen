!
! Routine to compute line opacity and emissivity. Routine can be
! called in dTdR section, T section, or in the continuum calculation section.
! After the call storage locations for line quantities (but not variation)
! have been allocated. 
! Returned:
!         EINA etc 
!         CHIL_MAT, ETAL_MAT
!         LINE_PROF_SIM, LINE_QW_SIM
!	  LINE_OPAC_CON, LINE_EMIS_CON
!         L_STAR_RATIO, U_STAR_RATION
!
        SUBROUTINE SET_LINE_OPAC(POPS,NU,FREQ_INDX,LAST_LINE,N_LINE_FREQ,
	1               LST_DEPTH_ONLY,LUER,ND,NT,NCF,MAX_SIM)
	USE MOD_CMFGEN
	USE OPAC_MOD
	USE CONTROL_VARIABLE_MOD
	USE LINE_VEC_MOD
	USE LINE_MOD
        IMPLICIT NONE
!
! Altered 20-Feb-2006 : Minor bug fix --- incorrect acces VAR_IN_USE_CNT when BA not being computed
! Created 21-Dec-2004
!
	INTEGER ND
	INTEGER NT
	INTEGER NCF
	INTEGER MAX_SIM
	INTEGER FIRST_LINE 
	INTEGER LAST_LINE
	INTEGER N_LINE_FREQ
	INTEGER LUER
!
	REAL*8 NU(NCF)
	REAL*8 POPS(NT,ND)
	LOGICAL LST_DEPTH_ONLY
!
	COMMON/LINE/ OPLIN,EMLIN
	REAL*8 OPLIN,EMLIN
!
	REAL*8 T1,T2
	REAL*8 NU_DOP
	REAL*8 FL
!
	INTEGER I,J,K,L
	INTEGER ID
	INTEGER FREQ_INDX
	INTEGER D_ST
!
	INTEGER NL,NUP
	INTEGER MNL,MNL_F
	INTEGER MNUP,MNUP_F
!
! Inilization section:
!
	FL=NU(FREQ_INDX)
!
! If we are computing dTdR, we only need to evaluate the opacity at the LASt depth.
!
	D_ST=1
	IF(LST_DEPTH_ONLY)D_ST=ND
!
! We use NEW_LINE_STORAGE to inidcate which storage locations are to be initialized.
!
	NEW_LINE_STORAGE(:)=.FALSE.
!
! Ensure none of the storage locations for the line, for J, or for the 
! variation of J with CHIL etc are being pointed at.
!
! WEAK_LINE is only accessed in the CONTINUUM section, and is reset. By setting it
! TRUE, we avoid accessing VAR_IN_SE and VAR_LEV_ID.
!
	IF(FREQ_INDX .EQ. 1)THEN
          LINE_STORAGE_USED(1:MAX_SIM)=.FALSE.
          LOW_POINTER(1:MAX_SIM)=0
          UP_POINTER(1:MAX_SIM)=0
	  WEAK_LINE(1:MAX_SIM)=.TRUE.
!
          VAR_IN_USE_CNT(:)=0
          VAR_LEV_ID(:)=0
          IMP_TRANS_VEC(:)=.FALSE.
!
          NUM_OF_WEAK_LINES=0.0D0
	END IF
!
! 
!
! Description of main vectors:
!
! LINES_THIS_FREQ --- Logical vector [NCF] indicating whether this frequency
!                        is part of the resonance zone (i.e. Doppler profile) of 
!                        one (or more) lines.
!
! LINE_ST_INDX_IN_NU --- Integer vector [N_LINES] which specifies the starting
!                          frequency index for this lines resonance zone.
!
! LINE_END_INDX_IN_NU --- Integer vector [N_LINES] which specifies the final
!                          frequency index for this lines resonance zone.
!
! FIRST_LINE   ---- Integer specifying the index of the highest frequency
!                         line which we are taking into account in the
!                         transfer.
!
! LAST_LINE  ---- Integer specifying the index of the lowest frequency
!                         line which we are taking into account in the
!                         transfer.
!                                        
! LINE_LOC   ---- Integer array. Used to locate location of a particular line
!                         in the SIM vectors/arrays.
!
! SIM_LINE_POINTER --- Integer array --- locates the line corresponding to
!                         the indicated storage location in the SIM vectors/
!                         arrays.
!
! Check whether we have to treat another line. We use a DO WHILE, rather
! than an IF statement, to handle lines which begin at the same (upper)
! frequency.
!
	DO WHILE( LAST_LINE .LT. N_LINE_FREQ .AND.
	1                FREQ_INDX .EQ. LINE_ST_INDX_IN_NU(LAST_LINE+1) )
!
! Have another line --- need to find its storage location.
!
	  I=1
	  DO WHILE(LINE_STORAGE_USED(I))
	    I=I+1
	    IF(I .GT. MAX_SIM)THEN
	      FIRST_LINE=N_LINE_FREQ
	      DO SIM_INDX=1,MAX_SIM		!Not 0 as used!
	        FIRST_LINE=MIN(FIRST_LINE,SIM_LINE_POINTER(SIM_INDX)) 
	      END DO
	      IF( FREQ_INDX .GT. LINE_END_INDX_IN_NU(FIRST_LINE))THEN
!
! Free up storage location for line.
!
	        I=LINE_LOC(FIRST_LINE)
	        LINE_STORAGE_USED(I)=.FALSE.
	        SIM_LINE_POINTER(I)=0
!
! Free up storage location for variation with respect to lower and upper
! levels.
!
	        NL=LOW_POINTER(I)
	        IF(NL .NE. 0 .AND. .NOT. WEAK_LINE(I))THEN
	          VAR_IN_USE_CNT(NL)=VAR_IN_USE_CNT(NL)-1
	          IF(VAR_IN_USE_CNT(NL) .EQ. 0)VAR_LEV_ID(NL)=0
	        END IF
	        LOW_POINTER(I)=0
!
	        NUP=UP_POINTER(I)
	        IF(NUP .NE. 0 .AND. .NOT. WEAK_LINE(I))THEN
	          VAR_IN_USE_CNT(NUP)=VAR_IN_USE_CNT(NUP)-1
	          IF(VAR_IN_USE_CNT(NUP) .EQ. 0)VAR_LEV_ID(NUP)=0
	        END IF
	        UP_POINTER(I)=0
	      ELSE
	        WRITE(LUER,*)'Too many lines have overlapping '//
	1                      'resonance zones'
	        WRITE(LUER,*)'Current frequency is:',FL
	        STOP
	      END IF
	    END IF
	  END DO
	  SIM_INDX=I
	  LAST_LINE=LAST_LINE+1
	  LINE_STORAGE_USED(SIM_INDX)=.TRUE.
	  LINE_LOC(LAST_LINE)=SIM_INDX
	  SIM_LINE_POINTER(SIM_INDX)=LAST_LINE
	  NEW_LINE_STORAGE(SIM_INDX)=.TRUE.
!
! Have located a storage location. Now must compute all relevant quantities
! necessary to include this line in the transfer calculations.
!
	  SIM_NL(SIM_INDX)=VEC_NL(LAST_LINE)
	  SIM_NUP(SIM_INDX)=VEC_NUP(LAST_LINE)
	  NL=SIM_NL(SIM_INDX)
	  NUP=SIM_NUP(SIM_INDX)
!
	  EINA(SIM_INDX)=VEC_EINA(LAST_LINE)
	  OSCIL(SIM_INDX)=VEC_OSCIL(LAST_LINE)
	  FL_SIM(SIM_INDX)=VEC_FREQ(LAST_LINE)
!
	  TRANS_NAME_SIM(SIM_INDX)=TRIM(VEC_SPEC(LAST_LINE))
!
! This is a temporary measure. We currently set AMASS to AMASS_DOP for all
! species.
!
!	  AMASS_SIM(SIM_INDX)=AMASS_ALL(NL)
	  AMASS_SIM(SIM_INDX)=AMASS_DOP
!
! 
!
! Compute U_STAR_RATIO and L_STAR_RATIO which are used to switch from
! the opacity/emissivity computed with a FULL_ATOM to an equivalent form
! but written in terms of the SUPER-LEVELS. 
!
! L refers to the lower level of the transition.
! U refers to the upper level of the transition.
!
! At present we must treat each species separately (for those with both FULL
! and SUPER_LEVEL model atoms).
!
! MNL_F (MNUP_F) denotes the lower (upper) level in the full atom.
! MNL (MNUP) denotes the lower (upper) level in the super level model atom.
!
	MNL_F=VEC_MNL_F(LAST_LINE)
	MNUP_F=VEC_MNUP_F(LAST_LINE)
	DO K=D_ST,ND
	  L_STAR_RATIO(K,SIM_INDX)=1.0D0
	  U_STAR_RATIO(K,SIM_INDX)=1.0D0
	END DO
!
! T1 is used to represent b(level)/b(super level). If no interpolation of
! the b values in a super level has been performed, this ratio will be unity .
! This ratio is NOT treated in the linearization.
!
	  DO ID=1,NUM_IONS-1
	    IF(VEC_SPEC(LAST_LINE) .EQ.ION_ID(ID))THEN
	      MNL=ATM(ID)%F_TO_S_XzV(MNL_F)
	      MNUP=ATM(ID)%F_TO_S_XzV(MNUP_F)
	      DO K=D_ST,ND
	        IF(ATM(ID)%XzVLTE_F(MNL_F,K) .NE. 0.0D0 .AND. ATM(ID)%XzVLTE(MNL,K) .NE. 0.0D0)THEN
	          T1=(ATM(ID)%XzV_F(MNL_F,K)/ATM(ID)%XzVLTE_F(MNL_F,K)) /
	1             (ATM(ID)%XzV(MNL,K)/ATM(ID)%XzVLTE(MNL,K))
	          L_STAR_RATIO(K,SIM_INDX)=T1*ATM(ID)%W_XzV_F(MNUP_F,K)*
	1            ATM(ID)%XzVLTE_F(MNL_F,K)/ATM(ID)%XzVLTE(MNL,K)/
	1            ATM(ID)%W_XzV_F(MNL_F,K)
	          T2=(ATM(ID)%XzV_F(MNUP_F,K)/ATM(ID)%XzVLTE_F(MNUP_F,K)) /
	1               (ATM(ID)%XzV(MNUP,K)/ATM(ID)%XzVLTE(MNUP,K))
	          U_STAR_RATIO(K,SIM_INDX)=T2*ATM(ID)%XzVLTE_F(MNUP_F,K)/
	1              ATM(ID)%XzVLTE(MNUP,K)
	        END IF
	      END DO
	      GLDGU(SIM_INDX)=ATM(ID)%GXzV_F(MNL_F)/ATM(ID)%GXzV_F(MNUP_F)
	      TRANS_NAME_SIM(SIM_INDX)=TRIM(TRANS_NAME_SIM(SIM_INDX))//
	1     '('//TRIM(ATM(ID)%XzVLEVNAME_F(MNUP_F))//'-'//
	1         TRIM(ATM(ID)%XzVLEVNAME_F(MNL_F))//')'
	      EXIT
	    END IF
	  END DO
! 
!
! Compute line opacity and emissivity for this line.
!
	  T1=OSCIL(SIM_INDX)*OPLIN
	  T2=FL_SIM(SIM_INDX)*EINA(SIM_INDX)*EMLIN
	  NL=SIM_NL(SIM_INDX)
	  NUP=SIM_NUP(SIM_INDX)
	  DO I=D_ST,ND
	    CHIL_MAT(I,SIM_INDX)=T1*(L_STAR_RATIO(I,SIM_INDX)*POPS(NL,I)-
	1            GLDGU(SIM_INDX)*U_STAR_RATIO(I,SIM_INDX)*POPS(NUP,I))
	    ETAL_MAT(I,SIM_INDX)=T2*POPS(NUP,I)*U_STAR_RATIO(I,SIM_INDX)
	    IF(CHIL_MAT(I,SIM_INDX) .EQ. 0)THEN
	      CHIL_MAT(I,SIM_INDX)=0.01*T1*POPS(NL,I)*L_STAR_RATIO(I,SIM_INDX)
	      WRITE(LUER,*)'Zero line opacity in CMFGEN_SUB'
	      WRITE(LUER,*)'This needs to be fixed'
	      J=LEN_TRIM(TRANS_NAME_SIM(SIM_INDX))
	      WRITE(LUER,'(1X,A)')TRANS_NAME_SIM(SIM_INDX)(1:J)
	    END IF
	  END DO
!
	  LINE_OPAC_CON(SIM_INDX)=T1
	  LINE_EMIS_CON(SIM_INDX)=T2
!
	END DO	!Checking whether a  new line is being added.
!
! 
!
! Check whether current frequency is a resonance frequency for each line.
!
	DO SIM_INDX=1,MAX_SIM
	  RESONANCE_ZONE(SIM_INDX)=.FALSE.
	  END_RES_ZONE(SIM_INDX)=.FALSE.
	  IF(LINE_STORAGE_USED(SIM_INDX))THEN
	    L=SIM_LINE_POINTER(SIM_INDX)
	    IF( FREQ_INDX .GE. LINE_ST_INDX_IN_NU(L) .AND.
	1          FREQ_INDX .LT. LINE_END_INDX_IN_NU(L))THEN
	      RESONANCE_ZONE(SIM_INDX)=.TRUE.
	    ELSE IF(FREQ_INDX .EQ. LINE_END_INDX_IN_NU(L))THEN
 	      RESONANCE_ZONE(SIM_INDX)=.TRUE.
	      END_RES_ZONE(SIM_INDX)=.TRUE.
	    END IF
	  END IF
	END DO
!
! Compute Doppler profile. At present this is assumed, for simplicity, to be
! depth independent.
!
	T1=1.0D-15/1.77245385095516D0		!1.0D-15/SQRT(PI)
	DO SIM_INDX=1,MAX_SIM
	  IF(RESONANCE_ZONE(SIM_INDX))THEN
	    NU_DOP=FL_SIM(SIM_INDX)*12.85*SQRT( TDOP/AMASS_SIM(SIM_INDX) +
	1                        (VTURB/12.85)**2 )/2.998D+05
	    LINE_PROF_SIM(SIM_INDX)=EXP( -( (FL-FL_SIM(SIM_INDX))/
	1              NU_DOP )**2 )*T1/NU_DOP
	  ELSE
	    LINE_PROF_SIM(SIM_INDX)=0.0D0
	  END IF
	END DO                                    
!
! Compute the LINE quadrature weights. Defined so that JBAR= SUM[LINE_QW*J]
!
	IF(FREQ_INDX .EQ. 1)THEN
	  T1=(NU(1)-NU(2))*0.5D+15
	ELSE IF(FREQ_INDX .EQ. NCF)THEN
	  T1=(NU(NCF-1)-NU(NCF))*0.5D+15
	ELSE
	  T1=(NU(FREQ_INDX-1)-NU(FREQ_INDX+1))*0.5D+15
	END IF
	DO SIM_INDX=1,MAX_SIM
	  LINE_QW_SIM(SIM_INDX)=LINE_PROF_SIM(SIM_INDX)*T1
	END DO
!
! Ensure that LAST_LINE points to the next LINE that is going to be handled 
! in the BLANKETING portion of the code.
!
	DO WHILE(LAST_LINE .LT. N_LINE_FREQ.AND.
	1            VEC_TRANS_TYPE(LAST_LINE+1)(1:4) .NE. 'BLAN')
	       LAST_LINE=LAST_LINE+1
	END DO
!	   
	RETURN
	END
