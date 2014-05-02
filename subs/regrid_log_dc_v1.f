!
! Program reads in population estimates from a file and uses
! linear interpolation in the log plane to lay these estimates
! on the "new" radius grid. 
!
	SUBROUTINE REGRID_LOG_DC_V1(DHEN,R,ED,T,DI,CLUMP_FAC,EDGE,F_TO_S,INT_SEQ,
	1               POPATOM,N,ND,LU_IN,INTERP_OPTION,FILE_NAME)
	IMPLICIT NONE
!
! Altered 18-Jan-2014 - Fixed minor bug -- needed ABS when checkikng equality of R(1) and OLD_R(1).
! Altered 16-Oct-2012 - Added BA_TX_CONV and NO_TX_CONV.
! Altered: 18-Nov-2011: Big fix. Not all DPOP was being set and this could cause a crash
!                            when taking LOGS (levels unused).
! Created: 18-Dec-2010: This routine replace REGRIDWSC_V3, REGRIDB_ON_NE, REGRID_TX_R.
!                       Interplation options are 'R', 'ED', 'SPH_TAU', and 'RTX'.
!                       Routine handled INPUT files with departure coefficients, or
!                       with LOG(DCs). 
!
	INTEGER N,ND
	INTEGER LU_IN
	REAL*8 DHEN(N,ND)
	REAL*8 R(ND)
	REAL*8 T(ND)
	REAL*8 ED(ND)
	REAL*8 DI(ND)
	REAL*8 EDGE(N)
	REAL*8 POPATOM(ND)
	REAL*8 CLUMP_FAC(ND)
	INTEGER F_TO_S(N)
	INTEGER INT_SEQ(N)
	CHARACTER(LEN=*) INTERP_OPTION
!
	REAL*8 TAU(ND)
	REAL*8 NEW_ED(ND)
	REAL*8 NEW_X(ND)
	REAL*8 NEW_T(ND)
	REAL*8 LOG_TEN
!
	REAL*8, ALLOCATABLE :: DPOP(:,:)
	REAL*8, ALLOCATABLE :: OLD_CLUMP_FAC(:)
	REAL*8, ALLOCATABLE :: OLD_DI(:)
	REAL*8, ALLOCATABLE :: OLD_ED(:)
	REAL*8, ALLOCATABLE :: OLD_R(:)
	REAL*8, ALLOCATABLE :: OLD_T(:)
	REAL*8, ALLOCATABLE :: OLD_TAU(:)
	REAL*8, ALLOCATABLE :: TA(:)
	REAL*8, ALLOCATABLE :: TB(:)
!
	REAL*8, ALLOCATABLE :: OLD_X(:)
!
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
	INTEGER ERROR_LU,LUER,LUWARN,WARNING_LU
	EXTERNAL ERROR_LU,WARNING_LU
!
! Local Variables.
!
	REAL*8 T_EXCITE,FX,DELTA_T
	REAL*8 T1,T2
	REAL*8 RTAU1,RTAU1_OLD
!
	INTEGER I,J,K
	INTEGER NZ,NOLD,NDOLD
	INTEGER NX,NX_ST,NX_END
	INTEGER COUNT,IOS
	LOGICAL, SAVE :: FIRST=.TRUE.
	LOGICAL BAD_TX_CONV
	LOGICAL NO_TX_CONV(ND)
	LOGICAL TAKE_LOGS
	LOGICAL CHECK_DC
	LOGICAL CLUMP_PRES
	CHARACTER(LEN=80) STRING
	CHARACTER(LEN=*) FILE_NAME
!
	LUER=ERROR_LU()
	LUWARN=WARNING_LU()
!
! Read in values from previous model.
!
	OPEN(UNIT=LU_IN,STATUS='OLD',FILE=FILE_NAME,IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening ',TRIM(FILE_NAME),' in REGRID_LOG_DC_V1'
	  WRITE(LUER,*)'IOS=',IOS
	  STOP
	END IF
!
! Check whether the file has a record containing 'Format date'. Its presence
! effects the way we read the file. CLUMP_FAC is presently not uses in
! this routine, as we regrid in R.
!
! We now assume base 10 for the input -- there was a brief period with base e.
!
	I=0
	STRING=' '
	DO WHILE(INDEX(STRING,'!Format date') .EQ. 0 .AND. I .LE. 10)
	  I=I+1
	  READ(LU_IN,'(A)')STRING
	END DO
	CHECK_DC=.FALSE.
	CLUMP_PRES=.TRUE.
	IF( INDEX(STRING,'!Format date') .EQ. 0)THEN
	  REWIND(LU_IN)
	  CHECK_DC=.TRUE.
	  CLUMP_PRES=.FALSE.
	END IF
	IF( INDEX(STRING,'!Format date') .NE. 0 .AND. INDEX(STRING,'18-Nov-2010') .NE. 0)THEN
	  TAKE_LOGS=.FALSE.
	  LOG_TEN=1.0D0
	ELSE IF( INDEX(STRING,'!Format date') .NE. 0 .AND. INDEX(STRING,'10-Dec-2010') .NE. 0)THEN
	  TAKE_LOGS=.FALSE.
	  LOG_TEN=LOG(10.0D0)
	ELSE
	  TAKE_LOGS=.TRUE.
	END IF
!
! T1=RP_OLD
!
	READ(LU_IN,*)T1,T2,NOLD,NDOLD
!
! Allocate required memory.
!
	I=MAX(NDOLD,NOLD,N,ND)
	ALLOCATE (TA(I))
	ALLOCATE (TB(I))
!
	ALLOCATE (DPOP(NOLD,NDOLD))
	ALLOCATE (OLD_R(NDOLD))
	ALLOCATE (OLD_T(NDOLD))
	ALLOCATE (OLD_ED(NDOLD))
	ALLOCATE (OLD_DI(NDOLD))
	ALLOCATE (OLD_X(NDOLD))
	ALLOCATE (OLD_TAU(NDOLD))
	ALLOCATE (OLD_CLUMP_FAC(NDOLD)); OLD_CLUMP_FAC=1.0D0
!
! NZ defines the number of atomic levels which can be found by
! direct interpolation
!
	NZ=N
	IF(N .GT. NOLD)NZ=NOLD
	DO I=1,NDOLD
	  IF(CLUMP_PRES)THEN
	    READ(LU_IN,*)OLD_R(I),OLD_DI(I),OLD_ED(I),OLD_T(I),T1,T2,OLD_CLUMP_FAC(I)
	  ELSE
	    READ(LU_IN,*)OLD_R(I),OLD_DI(I),OLD_ED(I),OLD_T(I)
	  END IF
	  READ(LU_IN,*)(TA(J),J=1,NOLD)
	  DO J=1,NOLD                         !We use NOLD for when we take Logs.
	    DPOP(J,I)=TA(J)
	  END DO
	END DO
	CLOSE(LU_IN)
!
! Decide if DPOP refers to b, b-1 or log b. Convert all to log b (base e).
!
	IF( DABS( TA(NOLD) ) .LT. 0.2D0 .AND. CHECK_DC)DPOP=DPOP+1.0D0
	IF(TAKE_LOGS)THEN
	  DPOP=DLOG(DPOP)
	ELSE IF(LOG_TEN .NE. 1.0D0)THEN
	  DPOP=DPOP*LOG_TEN
	END IF
!
	NX_ST=1
	NX_END=ND
	IF(INTERP_OPTION .EQ. 'R')THEN
	  IF(DABS(OLD_R(NDOLD)/R(ND)-1.0D0) .GT. 0.0001D0)THEN
	    IF(FIRST)THEN
	      WRITE(LUWARN,*)'Warning - core radius not identical in REGRID_LOG_DC_V1'
	      WRITE(LUER,*)'Warning - core radius not identical in REGRID_LOG_DC_V1'
	      WRITE(LUER,*)'Rescaling to make Rcore identical --- ',TRIM(FILE_NAME)
	      WRITE(LUER,*)'Additional warnings will be output to WARNINGS'
	      FIRST=.FALSE.
	    ELSE
	      WRITE(LUWARN,*)'Rescaling to make Rcore identical --- ',TRIM(FILE_NAME)
	    END IF
	    DO I=1,NDOLD
	      OLD_R(I)=R(ND)*( OLD_R(I)/OLD_R(NDOLD) )
	    END DO
	    OLD_R(NDOLD)=R(ND)
	  ELSE
	    OLD_R(NDOLD)=R(ND)
	  END IF
	  IF( ABS(1.0D0-OLD_R(1)/R(1)) .LE. 1.0D-10 )OLD_R(1)=R(1)
	  IF(OLD_R(2) .GE. OLD_R(1))THEN
	    WRITE(LUER,*)'Reset OLD_R(1) in REGRID_T_ED but now OLD_R(2) .GE. OLD_R(1))'
	    STOP
	  END IF
!
	  OLD_X=LOG(OLD_R)
	  NEW_X=LOG(R)
	  DO WHILE(R(NX_ST) .GT. OLD_R(1))
	    NX_ST=NX_ST+1
	  END DO
	  NX=ND-NX_ST+1
	ELSE IF(INTERP_OPTION .EQ. 'ED')THEN
	  OLD_X=LOG(OLD_ED*OLD_CLUMP_FAC)
	  NEW_X=LOG(ED)
          DO WHILE (NEW_X(NX_ST) .LT. OLD_X(1))
            NX_ST=NX_ST+1
          END DO
          NX_END=ND
          DO WHILE (NEW_X(NX_END) .GT. OLD_X(NDOLD))
            NX_END=NX_END-1
          END DO
          NX=NX_END-NX_ST+1
	ELSE IF(INTERP_OPTION .EQ. 'SPH_TAU')THEN
!
! Determine radius at which optical depth is unity. We use NEW_ED as a
! temporary vector for ED*CLUMP_FAC.
!
          NEW_ED(1:ND)=ED(1:ND)*CLUMP_FAC(1:ND)
          TAU(1)=6.65D-15*NEW_ED(1)*R(1)
          K=1
	  DO I=2,ND
            TAU(I)=TAU(I-1)+6.65D-15*(NEW_ED(I-1)+NEW_ED(I))*(R(I-1)-R(I))*0.5D0
            IF(TAU(I) .LE. 1.0D0)K=I
          END DO
	  IF(K .EQ. 1 .OR. K .EQ. ND)THEN
	    LUER=ERROR_LU()
            WRITE(LUER,*)'Error computing RTAU1 in REGRID_LOG_DC_V1'
          END IF
	  T1=(1.0D0-TAU(K))/(TAU(K+1)-TAU(K))
          RTAU1=T1*R(K+1)+(1.0D0-T1)*R(K)
!
! Compute the spherical optical depth scale. Assume atmosphere has constant
! V at outer boundary.
!
          DO I=1,ND
            NEW_ED(I)=NEW_ED(I)*RTAU1*RTAU1/R(I)/R(I)
          END DO
          TAU(1)=6.65D-15*NEW_ED(1)*R(1)/3.0D0
          DO I=2,ND
            TAU(I)=TAU(I-1)+6.65D-15*(NEW_ED(I-1)+NEW_ED(I))*(R(I-1)-R(I))*0.5D0
          END DO
!
! Determine radius at which optical depth is unity in old model. 
!
	  OLD_ED=OLD_ED*OLD_CLUMP_FAC
          OLD_TAU(1)=6.65D-15*OLD_ED(1)*OLD_R(1)
          K=1
	  DO I=2,NDOLD
            OLD_TAU(I)=OLD_TAU(I-1)+6.65D-15*(OLD_ED(I-1)+OLD_ED(I))*(OLD_R(I-1)-OLD_R(I))*0.5D0
            IF(OLD_TAU(I) .LE. 1.0D0)K=I
          END DO
	  IF(K .EQ. 1 .OR. K .EQ. ND)THEN
	    LUER=ERROR_LU()
            WRITE(LUER,*)'Error computing RTAU1_OLD in REGRID_LOG_DC_V1'
          END IF
          T1=(1.0D0-OLD_TAU(K))/(OLD_TAU(K+1)-OLD_TAU(K))
          RTAU1_OLD=T1*OLD_R(K+1)+(1.0D0-T1)*OLD_R(K)
!
! Compute the spherical optical depth scale.
!
          DO I=1,NDOLD
            OLD_ED(I)=OLD_ED(I)*(RTAU1_OLD/OLD_R(I))**2
          END DO
          OLD_TAU(1)=6.65D-15*OLD_ED(1)*OLD_R(1)/3.0D0
          DO I=2,NDOLD
            OLD_TAU(I)=OLD_TAU(I-1)+6.65D-15*(OLD_ED(I-1)+OLD_ED(I))*(OLD_R(I-1)-OLD_R(I))*0.5D0
          END DO
!
! Detrmine whether new mesh extends beyond oldmesh.
! NX is used to define the region over which we may use a linear
! interpolation. Beyond the old radius mesh we set the departure coefficients
! constant or equal to 1.
!
	  NX_ST=1
	  DO WHILE (TAU(NX_ST) .LT. OLD_TAU(1))
	    NX_ST=NX_ST+1
	  END DO
	  NX_END=ND
	  DO WHILE (TAU(NX_END) .GT. OLD_TAU(NDOLD))
	    NX_END=NX_END-1
	  END DO
	  NX=NX_END-NX_ST+1
!
	  OLD_X=LOG(OLD_TAU)
	  NEW_X=LOG(TAU)
	END IF
!
! Interpolate the ion density.
!
	DO I=1,NDOLD
	  OLD_DI(I)=LOG(OLD_DI(I))
	END DO
	CALL LINPOP(NEW_X(NX_ST),DI(NX_ST),NX,OLD_X,OLD_DI,NDOLD)
	DO I=NX_ST,NX_END
	  DI(I)=EXP(DI(I))
	END DO
!
! For locations outside grid, we assume ion density is fixed relative to the atom density.
!
        DO I=1,NX_ST-1
          DI(I)=POPATOM(I)*DI(NX_ST)/POPATOM(NX_ST)
        END DO
	DO I=NX_END+1,ND
	  DI(I)=POPATOM(I)*DI(NX_END)/POPATOM(NX_END)
	END DO
!
! Compute LOG(excitation temperatures), if INTERP_OPTION is RTX.
!
	IF(INTERP_OPTION .EQ. 'RTX')THEN
	  DO I=1,NZ
	    T_EXCITE=OLD_T(NDOLD)
	    DO J=NDOLD,1,-1
	      DELTA_T=100
	      DO WHILE(ABS(DELTA_T/T_EXCITE) .GT. 1.0D-08)
	        FX=DPOP(I,J)*EXP(HDKT*EDGE(I)*(1.0D0/OLD_T(J)-1.0D0/T_EXCITE))*
	1          (T_EXCITE/OLD_T(J))**1.5D0
	        DELTA_T=(FX-1.0D0)*T_EXCITE/FX/(1.5D0+HDKT*EDGE(I)/T_EXCITE)
	        IF(DELTA_T .GT.  0.8D0*T_EXCITE)DELTA_T=0.8D0*T_EXCITE
	        IF(DELTA_T .LT. -0.8D0*T_EXCITE)DELTA_T=-0.8D0*T_EXCITE
	        T_EXCITE=T_EXCITE-DELTA_T
	      END DO
	      DPOP(I,J)=DLOG(T_EXCITE)
	    END DO
	  END DO
	END IF
!
! Interpolate the departure coefficients / excitation temperatures.
! As DPOP is in LOG plane, there is no need to take LOGS.
!
	DO J=1,NZ
	  DO I=1,NDOLD
	    TB(I)=DPOP(J,I)
	  END DO
	  CALL LINPOP(NEW_X(NX_ST),TA(NX_ST),NX,OLD_X,TB,NDOLD)
	  DO I=NX_ST,NX_END
	    DHEN(J,I)=TA(I)
	  END DO
	END DO
!
! Convert, if necessary, back from T_EXCITE to log(b).
!
	IF(INTERP_OPTION .EQ. 'RTX')THEN
	  DO I=1,NDOLD
	    TA(I)=LOG(OLD_T(I))
	  END DO
	  CALL LINPOP(NEW_X(NX_ST),NEW_T(NX_ST),NX,OLD_X,TA,NDOLD)
	  DO I=NX_ST,NX_END
	    NEW_T(I)=EXP(NEW_T(I))
	  END DO
          DO I=NX_ST,NX_END
            T1=T(I)+(EXP(DHEN(J,I))-NEW_T(I))
            DHEN(J,I)=HDKT*EDGE(J)*(1.0D0/T1-1.0D0/T(I)) + 1.5D0*LOG(T(I)/T1)
          END DO
	END IF
!
! When extending to larger radii, we assume that the ionization
! state of the gas does not change. This supersedes the earlier assumption
! of a fixed departure coefficient. Excited populations are adjusted
! logarithmically. Code assumes ED is proportional to POP_ATOM, and
! hence may break down for SN models (with variable composition).
!
	DO I=1,NX_ST-1
	  DO J=1,NZ
	    T1=POPATOM(NX_ST)/POPATOM(I)
	    T2=MIN(1.0D0, ABS(DHEN(J,NX_ST))/ABS(DHEN(1,NX_ST)) )
	    DHEN(J,I)=DHEN(J,NX_ST)+T2*LOG(T1)
	  END DO
	END DO
	DO I=NX_END+1,ND
	  DO J=1,NZ
	    DHEN(J,I)=DHEN(J,NX_END)
	  END DO
	END DO
!
! Compute departure coefficients for N>NZ. These levels are set to have these
! same excitation temperature as the highest level.
!
	NO_TX_CONV=.FALSE.
	BAD_TX_CONV=.FALSE.
	IF(N .GT. NZ)THEN
	  T_EXCITE=T(ND)
	  DO I=ND,1,-1
!
! We first compute the excitation temperature on level NZ.
!
	    DELTA_T=100.0D0
	    COUNT=0
	    COUNT=1000
	    DO WHILE(ABS(DELTA_T) .GT. 1.0D-08 .AND. COUNT .LT. 100)
	      FX=EXP(DHEN(NZ,I)+HDKT*EDGE(NZ)*(1.0D0/T(I)-1.0D0/T_EXCITE))*(T_EXCITE/T(I))**1.5D0
	      DELTA_T=(FX-1.0D0)*T_EXCITE/FX/(1.5D0+HDKT*EDGE(NZ)/T_EXCITE)
	      IF(DELTA_T .GT.  0.8D0*T_EXCITE)DELTA_T=0.8D0*T_EXCITE
	      IF(DELTA_T .LT. -0.8D0*T_EXCITE)DELTA_T=-0.8D0*T_EXCITE
	      T_EXCITE=T_EXCITE-DELTA_T
	      COUNT=COUNT+1
	    END DO
!
! We can now compute the Departure coefficients.
!
	    IF(COUNT .GE. 100)THEN
	      NO_TX_CONV(I)=.TRUE.
	      BAD_TX_CONV=.TRUE.
	      DO J=NZ+1,N
	        DHEN(J,I)=DHEN(NZ,I)
	      END DO
	    ELSE
	      DO J=NZ+1,N
	        DHEN(J,I)=HDKT*EDGE(J)*(1.0D0/T_EXCITE-1.0D0/T(I))+1.5D0*LOG(T(I)/T_EXCITE)
	      END DO
	    END IF
	  END DO
	  IF(BAD_TX_CONV)THEN
	    WRITE(LUER,*)'Error in REGRID_LOG_DC_V1 - TX did not converge in 100 iterations'
	    WRITE(LUER,*)'File is ',TRIM(FILE_NAME)
	    WRITE(LUER,*)'It occurred at the following depths:'
	    DO I=1,ND
	      IF(NO_TX_CONV(I))WRITE(LUER,'(I5)',ADVANCE='NO')I
	      IF(MOD(I,10) .EQ. 0)WRITE(LUER,'(A)')' '
	    END DO
          END IF
	END IF
!
	CLOSE(UNIT=8)
!
! Ensure that all levels belonging to the same level have the
! same DC. We take the DC of the lowest state.
!
	IF(SUM(INT_SEQ) .EQ. 0)THEN
	  DO I=1,ND
	    TA(:)=0.0D0
	    DO J=1,N
	      IF(TA(F_TO_S(J)) .EQ. 0.0D0)TA(F_TO_S(J))=DHEN(J,I)
	    END DO
	    DO J=1,N
	      DHEN(J,I)=TA(F_TO_S(J))
	    END DO
	  END DO
	END IF
!
	WRITE(170,'(A)')FILE_NAME
	WRITE(170,'(5ES16.6)')(DHEN(1:N,1))
!
! Free up memory.
!
	DEALLOCATE (DPOP)
	DEALLOCATE (OLD_CLUMP_FAC)
	DEALLOCATE (OLD_ED)
	DEALLOCATE (OLD_DI)
	DEALLOCATE (OLD_R)
	DEALLOCATE (OLD_T)
	DEALLOCATE (OLD_TAU)
	DEALLOCATE (OLD_X)
	DEALLOCATE (TA)
	DEALLOCATE (TB)
!
	RETURN
	END
