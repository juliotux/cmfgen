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
! Incorporated 2-Jan-2014: Changes for depth depndent line profiles.
! Altered 05-Apr-2011 : L_STAR_RATIO and U_STAR_RATIO now computed using XzVLTE_F_ON_S (29-Nov-2010).
!                         Done to facilitate use of lower temperaturs.
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
	REAL*8 TA(ND),TB(ND),TC(ND),ED_MOD(ND)		!Work vectors
!
! Constants for opacity etc. These are set in CMFGEN.
!
        COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
        COMMON/LINE/ OPLIN,EMLIN
        REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
        REAL*8 OPLIN,EMLIN
!
	REAL*8 T1,T2,T3
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
	      T1=1.0D+08
	      DO SIM_INDX=1,MAX_SIM             !Not 0 as used!
	        IF(LINE_END_INDX_IN_NU(SIM_LINE_POINTER(SIM_INDX)) .LT. T1)THEN
	          FIRST_LINE=SIM_LINE_POINTER(SIM_INDX)
	          T1=LINE_END_INDX_IN_NU(SIM_LINE_POINTER(SIM_INDX))
                END IF
              END DO
!	      FIRST_LINE=N_LINE_FREQ
!	      DO SIM_INDX=1,MAX_SIM		!Not 0 as used!
!	        FIRST_LINE=MIN(FIRST_LINE,SIM_LINE_POINTER(SIM_INDX)) 
!	      END DO
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
	  IF(FIX_DOP)THEN
	    AMASS_SIM(SIM_INDX)=AMASS_DOP
	  ELSE
	    AMASS_SIM(SIM_INDX)=AT_MASS(SPECIES_LNK(VEC_ID(LAST_LINE)))
	  END IF
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
	        T1=(ATM(ID)%XzV_F(MNL_F,K)/ATM(ID)%XzV(MNL,K))/ATM(ID)%XzVLTE_F_ON_S(MNL_F,K)
	        L_STAR_RATIO(K,SIM_INDX)=T1*(ATM(ID)%W_XzV_F(MNUP_F,K)/ATM(ID)%W_XzV_F(MNL_F,K))*ATM(ID)%XzVLTE_F_ON_S(MNL_F,K)
	        T2=(ATM(ID)%XzV_F(MNUP_F,K)/ATM(ID)%XzV(MNUP,K))/ATM(ID)%XzVLTE_F_ON_S(MNUP_F,K)
	        U_STAR_RATIO(K,SIM_INDX)=T2*ATM(ID)%XzVLTE_F_ON_S(MNUP_F,K)
	        dL_RAT_dT(K,SIM_INDX)=L_STAR_RATIO(K,SIM_INDX)*
	1         (-1.5D0-HDKT*ATM(ID)%EDGEXzV_F(MNL_F)/T(K)-ATM(ID)%dlnXzVLTE_dlnT(MNL,K))/T(K)
	        dU_RAT_dT(K,SIM_INDX)=U_STAR_RATIO(K,SIM_INDX)*
	1         (-1.5D0-HDKT*ATM(ID)%EDGEXzV_F(MNUP_F)/T(K)-ATM(ID)%dlnXzVLTE_dlnT(MNUP,K))/T(K)
	      END DO
	      GLDGU(SIM_INDX)=ATM(ID)%GXzV_F(MNL_F)/ATM(ID)%GXzV_F(MNUP_F)
	      TRANS_NAME_SIM(SIM_INDX)=TRIM(TRANS_NAME_SIM(SIM_INDX))//
	1     '('//TRIM(ATM(ID)%XzVLEVNAME_F(MNUP_F))//'-'//
	1         TRIM(ATM(ID)%XzVLEVNAME_F(MNL_F))//')'
	      EXIT
	    END IF
	  END DO
!
	  IF(.NOT. INCLUDE_dSLdT)THEN
	    DO K=D_ST,ND
	      dU_RAT_dT(K,SIM_INDX)=0.0D0
	      dL_RAT_dT(K,SIM_INDX)=0.0D0
	    END DO
	  END IF
! 
!
! If desired we can scale the line opacity/emissivities of transitions
! among SL's. This allows us to obtain full consistency in the
! upward/downwrd rates, and the cooling/heating rates. For the
! cooling/heting rates, this is similar to SCL_LIN_COOL_RATES option.
! However this option provides consistent with the radiative equilibrium
! term appearing in the transfer equation. Implemented for SN which can
! be totally dominated by scattering.
!
! This section should be consistent with corrections in CMFGEN_SUB
! (in SCL_LINE_COO_RATES sections).
!
!	  IF(SCL_SL_LINE_OPAC)THEN
!	    T3=(AVE_ENERGY(NL)-AVE_ENERGY(NUP))/FL_SIM(SIM_INDX)
!	    IF(ABS(T3-1.0D0) .GT. SCL_LINE_HT_FAC)T3=1.0D0
!	    DO I=D_ST,ND
!	      IF(POP_ATOM(I) .LE. SCL_LINE_DENSITY_LIMIT)THEN
!	        L_STAR_RATIO(I,SIM_INDX)=T3*L_STAR_RATIO(I,SIM_INDX)
!	        U_STAR_RATIO(I,SIM_INDX)=T3*U_STAR_RATIO(I,SIM_INDX)
!	      END IF
!	    END DO
!	  END IF
!
! Compute line opacity and emissivity for this line.
!
	  T3=1.0D0
	  IF(SCL_SL_LINE_OPAC)THEN
	    T3=(AVE_ENERGY(NL)-AVE_ENERGY(NUP))/FL_SIM(SIM_INDX)
	    IF(ABS(T3-1.0D0) .GT. SCL_LINE_HT_FAC)T3=1.0D0
	  END IF
	  T1=OSCIL(SIM_INDX)*OPLIN*T3
	  T2=FL_SIM(SIM_INDX)*EINA(SIM_INDX)*EMLIN*T3
	  LINE_OPAC_CON(SIM_INDX)=T1
	  LINE_EMIS_CON(SIM_INDX)=T2
!
	  NL=SIM_NL(SIM_INDX)
	  NUP=SIM_NUP(SIM_INDX)
	  DO I=D_ST,ND
	    CHIL_MAT(I,SIM_INDX)=T1*(L_STAR_RATIO(I,SIM_INDX)*POPS(NL,I)-
	1            GLDGU(SIM_INDX)*U_STAR_RATIO(I,SIM_INDX)*POPS(NUP,I))
	    ETAL_MAT(I,SIM_INDX)=T2*POPS(NUP,I)*U_STAR_RATIO(I,SIM_INDX)
	    IF(CHIL_MAT(I,SIM_INDX) .EQ. 0)THEN
	      CHIL_MAT(I,SIM_INDX)=0.01D0*T1*POPS(NL,I)*L_STAR_RATIO(I,SIM_INDX)
	      WRITE(LUER,*)'Zero line opacity in CMFGEN_SUB'
	      WRITE(LUER,*)'This needs to be fixed'
	      J=LEN_TRIM(TRANS_NAME_SIM(SIM_INDX))
	      WRITE(LUER,'(1X,A)')TRANS_NAME_SIM(SIM_INDX)(1:J)
	    END IF
	  END DO
!
	END DO	!Checking whether a new line is being added.
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
! Compute intrinsic line profile. This can be depth independent (simplest option) or
! depth dependent. Note: AMASS_SIM has been set to AMASS_DOP for FIX_DOP.

	IF(FIX_DOP .OR. GLOBAL_LINE_PROF .EQ. 'DOP_SPEC')THEN
	  T1=1.0D-15/1.77245385095516D0		!1.0D-15/SQRT(PI)
	  DO SIM_INDX=1,MAX_SIM
	    IF(RESONANCE_ZONE(SIM_INDX))THEN
	      NU_DOP=FL_SIM(SIM_INDX)*12.85D0*SQRT( TDOP/AMASS_SIM(SIM_INDX) + (VTURB/12.85D0)**2 )/2.998D+05
	      LINE_PROF_SIM(ND,SIM_INDX)=EXP( -( (FL-FL_SIM(SIM_INDX))/NU_DOP )**2 )*T1/NU_DOP
	    ELSE
	      LINE_PROF_SIM(ND,SIM_INDX)=0.0D0
	    END IF
	    DO I=D_ST,ND-1
	      LINE_PROF_SIM(I,SIM_INDX)=LINE_PROF_SIM(ND,SIM_INDX)
	    END DO
	  END DO                      
	ELSE
!
! Because of storage issues, need to compute all ND profiles. Thus there is
! currently no LST_DEPTH option.
!
! We also have a temporary limit on ED to prvent the Stark profile from becoming
! too large. This may need to change.
!
	  TB(1:ND)=0.0D0; TC(1:ND)=0.0D0
	  ED_MOD(1:ND)=ED(1:ND)
	  DO I=1,ND
	    ED_MOD(I)=MIN(15.0D0,ED_MOD(I))
	  END DO
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES .AND. ION_ID(ID) .EQ. 'HI')TB(1:ND)=ATM(ID)%DxzV(1:ND)
	    IF(ATM(ID)%XzV_PRES .AND. ION_ID(ID) .EQ. 'HeI')TC(1:ND)=ATM(ID)%DxzV(1:ND)
	  END DO
	  DO SIM_INDX=1,MAX_SIM
	    IF(RESONANCE_ZONE(SIM_INDX))THEN
	      J=SIM_LINE_POINTER(SIM_INDX); I=FREQ_INDX
	      ID=VEC_ID(J); T1=ATM(ID)%ZXzV+1; T3=0.0D0
	      CALL SET_PROF_V4(TA,NU,I,
	1               LINE_ST_INDX_IN_NU(J),LINE_END_INDX_IN_NU(J),
	1               ED,TB,TC,T,VTURB_VEC,ND,
	1               PROF_TYPE(J),PROF_LIST_LOCATION(J),
	1               VEC_FREQ(J),VEC_MNL_F(J),VEC_MNUP_F(J),
	1               AMASS_SIM(SIM_INDX),T1,VEC_ARAD(J),T3,
	1               TDOP,AMASS_DOP,VTURB,
	1               END_RES_ZONE(SIM_INDX),L_TRUE,7)
	      LINE_PROF_SIM(1:ND,SIM_INDX)=TA(1:ND)
	      IF(VEC_SPEC(J)(1:1) .EQ. 'H')THEN
	        WRITE(135,'(A,T10,2ES14.5,2I4,3E12.4)')PROF_TYPE(J),VEC_FREQ(J),
	1            3.0D+05*(NU(I)/VEC_FREQ(J)-1.0D0),
	1            VEC_MNL_F(J),VEC_MNUP_F(J),TA(1),TA(40),TA(ND)
	      END IF
	    ELSE
	      LINE_PROF_SIM(1:ND,SIM_INDX)=0.0D0
	    END IF
	  END DO
	END IF
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
!
	IF(LST_DEPTH_ONLY)THEN
	  DO SIM_INDX=1,MAX_SIM
	    LINE_QW_SIM(ND,SIM_INDX)=LINE_PROF_SIM(ND,SIM_INDX)*T1
	  END DO
	ELSE IF(FIX_DOP .OR. GLOBAL_LINE_PROF .EQ. 'DOP_SPEC')THEN
!$OMP PARALLEL DO PRIVATE(SIM_INDX,I)
	  DO SIM_INDX=1,MAX_SIM
	    LINE_QW_SIM(1,SIM_INDX)=LINE_PROF_SIM(1,SIM_INDX)*T1
	    DO I=2,ND
	      LINE_QW_SIM(I,SIM_INDX)=LINE_QW_SIM(1,SIM_INDX)
	    END DO
	  END DO
	ELSE
!$OMP PARALLEL DO PRIVATE(SIM_INDX,I)
	  DO SIM_INDX=1,MAX_SIM
	    DO I=1,ND
	      LINE_QW_SIM(I,SIM_INDX)=LINE_PROF_SIM(I,SIM_INDX)*T1
	    END DO
	  END DO
	END IF
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
