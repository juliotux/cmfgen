!
! Program reads in population estimates from a file and uses
! linear interpolation in the log plane to lay these estimates
! on the "new" radius grid. 
!
	SUBROUTINE REGRID_TX_R(DHEN,DI,EDGE,F_TO_S,INT_SEQ,
	1                POPATOM,T_MOD,R,N,ND,FILNAME)
	IMPLICIT NONE
!
	INTEGER N,ND
	REAL*8 DHEN(N,ND)
	REAL*8 DI(ND)
!
	REAL*8 R(ND)
	REAL*8 T_MOD(ND)
	REAL*8 EDGE(N)
	REAL*8 POPATOM(ND)
	INTEGER F_TO_S(N)
	INTEGER INT_SEQ(N)
!
	REAL*8 T(ND)
	REAL*8 ED(ND)
!
	REAL*8, ALLOCATABLE :: DPOP(:,:)
	REAL*8, ALLOCATABLE :: OLDED(:)
	REAL*8, ALLOCATABLE :: OLDDI(:)
	REAL*8, ALLOCATABLE :: OLDT(:)
	REAL*8, ALLOCATABLE :: RLOG(:)
	REAL*8, ALLOCATABLE :: TA(:)
	REAL*8, ALLOCATABLE :: TB(:)
	REAL*8, ALLOCATABLE :: OLDR(:)
!
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
! Local Variables.
!
	INTEGER I,J,NOLD,NDOLD,COUNT,NX,NXST,NZ,IOS
	REAL*8 RPOLD,TX,T1,T2
	REAL*8 T_EXCITE
	REAL*8 FX
	REAL*8 DELTA_T
	LOGICAL CHECK_DC
	CHARACTER*80 STRING
	CHARACTER*(*) FILNAME
!
	LUER=ERROR_LU()
!
! Read in values from previous model.
!
	OPEN(UNIT=8,STATUS='OLD',FILE=FILNAME,IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening ',TRIM(FILNAME),' in REGRIDWSC_V3'
	  WRITE(LUER,*)'IOS=',IOS
	  STOP
	END IF
!
! Check whether the file has a record containing 'Format date'. Its presence
! effects the way we read the file. CLUMP_FAC is presently not uses in
! this routine, as we regrid in R.
!
	I=0
	STRING=' '
	DO WHILE(INDEX(STRING,'!Format date') .EQ. 0 .AND. I .LE. 10)
	  I=I+1
	  READ(8,'(A)')STRING
	END DO
	CHECK_DC=.FALSE.
	IF( INDEX(STRING,'!Format date') .EQ. 0)THEN
	  REWIND(8)
	  CHECK_DC=.TRUE.
	END IF
!
	READ(8,*)RPOLD,T1,NOLD,NDOLD
!
! Allocate required memory.
!
	I=MAX(NDOLD,NOLD,N,ND)
	ALLOCATE (DPOP(MAX(NOLD,N),NDOLD))
	ALLOCATE (OLDED(I))
	ALLOCATE (OLDDI(I))
	ALLOCATE (OLDT(I))
	ALLOCATE (RLOG(I))
	ALLOCATE (TA(I))
	ALLOCATE (TB(I))
	ALLOCATE (OLDR(I))
!
! NZ defines the number of atomic levels which can be found by
! direct interpolation
!
	NZ=N
	IF(N .GT. NOLD)NZ=NOLD
	DO I=1,NDOLD
	  READ(8,*)OLDR(I),OLDDI(I),OLDED(I),OLDT(I)
	  READ(8,*)(TA(J),J=1,NOLD)
	  DO J=1,NZ
	    DPOP(J,I)=TA(J)
	  END DO
	END DO
!
! Using the highest level read in, decide whether b, or b-1 have been read in.
!
	IF( DABS( TA(NOLD) ) .LT. 0.2D0 .AND. CHECK_DC)DPOP=DPOP+1.0D0
!
	IF(DABS(OLDR(NDOLD)/R(ND)-1.0D0) .GT. 0.0001D0)THEN
	  WRITE(LUER,*)'Warning - core radius not identical in REGRIDWSC_V3'
	  WRITE(LUER,*)'Rescaling to make Rcore identical'
	  DO I=1,NDOLD
	    OLDR(I)=R(ND)*( OLDR(I)/OLDR(NDOLD) )
	  END DO
	  OLDR(NDOLD)=R(ND)
	ELSE
	  OLDR(NDOLD)=R(ND)
	END IF
!
	NXST=1
	DO WHILE(R(NXST) .GT. OLDR(1))
	  NXST=NXST+1
	END DO
	NX=ND-NXST+1
!
! Compute excitation temperatures.
!
	DO I=1,NZ
	  T_EXCITE=OLDT(NDOLD)
	  DO J=NDOLD,1,-1
	    DELTA_T=100
	    DO WHILE(ABS(DELTA_T/T_EXCITE) .GT. 1.0E-08)
	      FX=DPOP(I,J)*EXP(HDKT*EDGE(I)*(1.0D0/OLDT(J)-1.0D0/T_EXCITE))*
	1         (T_EXCITE/OLDT(J))**1.5D0
	      DELTA_T=(FX-1.0D0)*T_EXCITE/FX/(1.5D0+HDKT*EDGE(I)/T_EXCITE)
	      IF(DELTA_T .GT.  0.8D0*T_EXCITE)DELTA_T=0.8D0*T_EXCITE
	      IF(DELTA_T .LT. -0.8D0*T_EXCITE)DELTA_T=-0.8D0*T_EXCITE
	      T_EXCITE=T_EXCITE-DELTA_T
	    END DO
	    DPOP(I,J)=T_EXCITE
	  END DO
	END DO
!
! If levels aren't present in input file, set excitation level to that
! of the last known level.
!
	IF(N .GT. NZ)THEN
	   DO I=NZ+1,N
	     DPOP(I,1:ND)=DPOP(NZ,1:ND)
	   END DO
	END IF
!
!Interpolations are performed in the log plane.
!
	DO I=1,ND
	  RLOG(I)=LOG(R(I))
	END DO
	DO I=1,NDOLD
	  OLDR(I)=LOG(OLDR(I))
	  OLDT(I)=LOG(OLDT(I))
	  OLDED(I)=LOG(OLDED(I))
	  OLDDI(I)=LOG(OLDDI(I))
	END DO
!
	CALL LINPOP(RLOG(NXST),T(NXST),NX,OLDR,OLDT,NDOLD)
	CALL LINPOP(RLOG(NXST),ED(NXST),NX,OLDR,OLDED,NDOLD)
	CALL LINPOP(RLOG(NXST),DI(NXST),NX,OLDR,OLDDI,NDOLD)

	DO I=NXST,ND
	  T(I)=EXP(T(I))
	  ED(I)=EXP(ED(I))
	  DI(I)=EXP(DI(I))
	END DO
!
	DO I=1,NXST-1
	  T(I)=T(NXST)
	  ED(I)=POPATOM(I)*ED(NXST)/POPATOM(NXST)
	  DI(I)=POPATOM(I)*DI(NXST)/POPATOM(NXST)
	END DO
!
! When extending to larger radii, we assume that the excitation
! temperature is constant.
!
	DO J=1,NZ
	  DO I=1,NDOLD
	    TA(I)=LOG(DPOP(J,I))
	  END DO
	  CALL LINPOP(RLOG(NXST),TB(NXST),NX,OLDR,TA,NDOLD)
	  DO I=NXST,ND
	    T1=T_MOD(I)+(EXP(TB(I))-T(I))
	    DHEN(J,I)=EXP(HDKT*EDGE(J)*(1.0D0/T1-1.0D0/T_MOD(I)))*(T_MOD(I)/T1)**1.5D0
	  END DO
	  DO I=1,NXST-1
	    T1=T_MOD(I)+(EXP(TB(NXST))-T(I))
	    DHEN(J,I)=EXP(HDKT*EDGE(J)*(1.0D0/T1-1.0D0/T_MOD(I)))*(T_MOD(I)/T1)**1.5D0
	  END DO
	END DO
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
! Free up memory.
!
	DEALLOCATE (DPOP)
	DEALLOCATE (OLDED)
	DEALLOCATE (OLDDI)
	DEALLOCATE (OLDT)
	DEALLOCATE (RLOG)
	DEALLOCATE (TA,TB)
	DEALLOCATE (OLDR)
!
	RETURN
	END
