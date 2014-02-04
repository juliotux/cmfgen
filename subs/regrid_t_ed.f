!
! Program reads in population estimates from a normal DC file,
! and uses linear interpolation in the log plane to deduce T and ED.
! If R(model) > R(model), R is scaled.
!
	SUBROUTINE REGRID_T_ED(R,ED,T,POPATOM,ND,FILNAME)
	IMPLICIT NONE
!
! Created 20-Sep-2006 : Based on REGRIDWSC.
!
	INTEGER ND
	REAL*8 R(ND)
	REAL*8 T(ND)
	REAL*8 ED(ND)
	REAL*8 POPATOM(ND)
!
	REAL*8, ALLOCATABLE :: TA(:)
	REAL*8, ALLOCATABLE :: OLDED(:)
	REAL*8, ALLOCATABLE :: OLDT(:)
	REAL*8, ALLOCATABLE :: RLOG(:)
	REAL*8, ALLOCATABLE :: OLDR(:)
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
! Local Variables.
!
	INTEGER I,J,IOS
	INTEGER NOLD,NDOLD
	INTEGER NX,NXST
	REAL*8 RPOLD,T1
	CHARACTER*80 STRING
	CHARACTER*(*) FILNAME
!
	LUER=ERROR_LU()
!
! Read in values from previous model.
!
	OPEN(UNIT=8,STATUS='OLD',FILE=FILNAME,IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening ',TRIM(FILNAME),' in REGRIDWSC'
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
	  IF( INDEX(STRING,'!Format date') .EQ. 0)REWIND(8)
!
	  READ(8,*)RPOLD,T1,NOLD,NDOLD
!
! Allocate required memory.
!
	  I=MAX(NDOLD,NOLD,ND)
	  ALLOCATE (OLDED(I))
	  ALLOCATE (OLDT(I))
	  ALLOCATE (RLOG(I))
	  ALLOCATE (TA(I))
	  ALLOCATE (OLDR(I))
!
	  DO I=1,NDOLD
	    READ(8,*)OLDR(I),TA(I),OLDED(I),OLDT(I)
	    READ(8,*)(TA(J),J=1,NOLD)
	  END DO
	CLOSE(UNIT=8)
!
!
	IF(DABS(OLDR(NDOLD)/R(ND)-1.0D0) .GT. 0.0001D0)THEN
	  WRITE(LUER,*)'Warning - core radius not identical in REGRIDWSC'
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
! Interpolations are performed in the log plane.
!
	DO I=1,ND
	  RLOG(I)=LOG(R(I))
	END DO
	DO I=1,NDOLD
	  OLDR(I)=LOG(OLDR(I))
	  OLDT(I)=LOG(OLDT(I))
	  OLDED(I)=LOG(OLDED(I))
	END DO
!
	CALL LINPOP(RLOG(NXST),T(NXST),NX,OLDR,OLDT,NDOLD)
	CALL LINPOP(RLOG(NXST),ED(NXST),NX,OLDR,OLDED,NDOLD)

	DO I=NXST,ND
	  T(I)=EXP(T(I))
	  ED(I)=EXP(ED(I))
	END DO
	DO I=1,NXST-1
	  T(I)=T(NXST)
	  ED(I)=ED(NXST)*POPATOM(I)/POPATOM(NXST)
	END DO
!
! Free up memory.
!
	DEALLOCATE (TA)
	DEALLOCATE (OLDED)
	DEALLOCATE (OLDT)
	DEALLOCATE (RLOG)
	DEALLOCATE (OLDR)
!
	RETURN
	END
