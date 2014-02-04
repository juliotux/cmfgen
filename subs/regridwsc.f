C
C Program reads in population estimates from a file and uses
C linear interpolation in the log plane to lay these estimates
C on the "new" radius grid. 
C
	SUBROUTINE REGRIDWSC(DHEN,R,ED,T,DI,EDGE,N,ND,FILNAME)
	IMPLICIT NONE
C
C Altered 15-Feb-2002 - Altered so that we can input files where the 
C                          outer radius ig reater than the outer radius of
C                          the new model.
C Altered 14-Jun-2000 - Error check when file opened.
C Altered 07-Jul-1997 - We check the file for a '!Format date' as this
C                          effect our free-format read.
C Altered 13-Nov-1996 - Bug fix. Typo in setting maximum allocation for
C                         vectors --- error occrrs when NOLD is the
C                         largest value.
C
C Altered 26-May-1996 - Access to scratch block removed. Dynamic memory
C                         allocation used.
C                       ERROR_LU installed.
C                       Generic call for OOG and EXP.
C
C Altered 02-Aug-1994 - PROGRDESC change from INTEGER to REAL*8. Corrects
C                       alignment problem. PROGDESC is now check for
C                       corruption, even if too many iterations occur.
C Altered 3-May-1989 - NV and check installed.
C
C Altered 29-Nov-1989 - EDGE installed to enable departure coefficients for
C                        unknow levels to be computed assuming the same
C                        excitation temperature as the last known level.
C                        NB - Call was changed. Implicit none installed.
C Altered 08-Oct-1987 - Handles departure coefficients, OR b-1.
C Altered 11-Mar-1987 - Reading sequence split into two calls so that V and ION
C                         may appear on the same line as R.
C Created 28-SEP-1985 - (Based on REGRID - scratch block added so the 
C                         NDOLD can be grater than ND)
C Altered 30-APR-1985 - (Bug fixed in error messages - import when ND .ne. NDOLD)
C Altered 28-FEB-1984
C
	INTEGER N,ND
	REAL*8 DHEN(N,ND),R(ND),T(ND),ED(ND),DI(ND),EDGE(N)
C
	REAL*8, ALLOCATABLE :: DPOP(:,:)
	REAL*8, ALLOCATABLE :: OLDED(:)
	REAL*8, ALLOCATABLE :: OLDDI(:)
	REAL*8, ALLOCATABLE :: OLDT(:)
	REAL*8, ALLOCATABLE :: RLOG(:)
	REAL*8, ALLOCATABLE :: TA(:)
	REAL*8, ALLOCATABLE :: OLDR(:)
C
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
C Local Variables.
C
	INTEGER I,J,NOLD,NDOLD,COUNT,NX,NXST,NZ,IOS
	REAL*8 RPOLD,TX,DELTA_T,T1,ADD1
	CHARACTER*80 STRING
	CHARACTER*(*) FILNAME
C
	LUER=ERROR_LU()
C
C Read in values from previous model.
C
	OPEN(UNIT=8,STATUS='OLD',FILE=FILNAME,IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening ',TRIM(FILNAME),' in REGRIDWSC'
	  WRITE(LUER,*)'IOS=',IOS
	  STOP
	END IF
C
C Check whether the file has a record containing 'Format date'. Its presence
C effects the way we read the file. CLUMP_FAC is presently not uses in
C this routine, as we regrid in R.
C
	I=0
	STRING=' '
	DO WHILE(INDEX(STRING,'!Format date') .EQ. 0 .AND. I .LE. 10)
	  I=I+1
	  READ(8,'(A)')STRING
	END DO
	IF( INDEX(STRING,'!Format date') .EQ. 0)REWIND(8)
C
	READ(8,*)RPOLD,T1,NOLD,NDOLD
C
C Allocate required memory.
C
	I=MAX(NDOLD,NOLD,N,ND)
	ALLOCATE (DPOP(NOLD,NDOLD))
	ALLOCATE (OLDED(I))
	ALLOCATE (OLDDI(I))
	ALLOCATE (OLDT(I))
	ALLOCATE (RLOG(I))
	ALLOCATE (TA(I))
	ALLOCATE (OLDR(I))
C
C NZ defines the number of atomic levels which can be found by
C direct interpolation
C
	NZ=N
	IF(N .GT. NOLD)NZ=NOLD
	DO I=1,NDOLD
	  READ(8,*)OLDR(I),OLDDI(I),OLDED(I),OLDT(I)
	  READ(8,*)(TA(J),J=1,NOLD)
	  DO J=1,NZ
	    DPOP(J,I)=TA(J)
	  END DO
	END DO
C
C Using the highest level read in, decide whether b, or b-1 have been read in.
C
	ADD1=0.0D0
	IF( DABS( TA(NOLD) ) .LT. 0.2 )ADD1=1.0D0
C
	IF(DABS(OLDR(NDOLD)/R(ND)-1.0D0) .GT. 0.0001)THEN
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
!	IF(R(1)/OLDR(1)-1.0D0 .GT. 0.0001)THEN
!	  WRITE(LUER,*)'Error in REGRIDWSC'
!	  WRITE(LUER,*)'Model boundary radius greater than that of input model'
!	  STOP
!	ELSE IF(R(2) .GT. OLDR(1))THEN
!	  WRITE(LUER,*)'Error in REGRIDWSC'
!	  WRITE(LUER,*)'Boundary radius R(2) too large'
!	  STOP
	IF(R(1) .GT. OLDR(1))THEN
	  OLDR(1)=R(1)
	END IF
	NXST=1
	NX=ND
C
C Interpolations are performed in the log plane.
C
	DO I=1,ND
	  RLOG(I)=LOG(R(I))
	END DO
	DO I=1,NDOLD
	  OLDR(I)=LOG(OLDR(I))
	  OLDT(I)=LOG(OLDT(I))
	  OLDED(I)=LOG(OLDED(I))
	  OLDDI(I)=LOG(OLDDI(I))
	END DO
C
	CALL LINPOP(RLOG(NXST),T(NXST),NX,OLDR,OLDT,NDOLD)
	CALL LINPOP(RLOG(NXST),ED(NXST),NX,OLDR,OLDED,NDOLD)
	CALL LINPOP(RLOG(NXST),DI(NXST),NX,OLDR,OLDDI,NDOLD)

	DO I=NXST,ND
	  T(I)=EXP(T(I))
	  ED(I)=EXP(ED(I))
	  DI(I)=EXP(DI(I))
	END DO
C
C Interpolate the populations assuming that the departure
C coefficients remain constant.
C
	DO J=1,NZ
	  DO I=1,NDOLD
	    OLDT(I)=LOG(DPOP(J,I)+ADD1)	!insures that never neg.
	  END DO
	  CALL LINPOP(RLOG(NXST),TA(NXST),NX,OLDR,OLDT,NDOLD)
	  DO I=NXST,ND
	    DHEN(J,I)=EXP(TA(I))
	  END DO
	END DO
C
C Compute departure coefficints for N>NZ . New levels are set to have thes
C same excitation temperature as the highest level/
C
	IF(N .GT. NZ)THEN
	  DO I=1,ND
C
C We first compute the excitation temperature on level NZ.
C
	    TX=T(I)
	    DELTA_T=10.0D0
	    COUNT=0
	    DO WHILE( ABS(DELTA_T) .GT. 1.0E-06 .AND. COUNT .LT. 100 )
	      T1=( (T(I)/TX)**(1.5) )*EXP(HDKT*EDGE(NZ)*(T(I)-TX)/T(I)/TX)
	      COUNT=COUNT+1
	      DELTA_T= (DHEN(NZ,I)-T1)*TX/T1/(1.0D0+HDKT*EDGE(NZ)/TX)
	      IF(ABS(DELTA_T) .GT. 0.8*TX)DELTA_T=0.5D0*DELTA_T
	      TX=TX-DELTA_T
	    END DO
	    IF(COUNT .EQ. 100)THEN
	      WRITE(LUER,*)'Error in REGRIDWSC - TX didnt converge ',
	1               'in 100 iterations'
	      WRITE(LUER,*)'I,J=',I,J
	      GO TO 1000
	    END IF
	    DO J=NZ+1,N
	      DHEN(J,I)=DHEN(NZ,I)*EXP( HDKT*(EDGE(NZ)-EDGE(J))*
	1                   (TX-T(I))/TX/T(I) )
	    END DO
	  END DO
	END IF
C
1000	CONTINUE
	CLOSE(UNIT=8)
C
C Free up memory.
C
	DEALLOCATE (DPOP)
	DEALLOCATE (OLDED)
	DEALLOCATE (OLDDI)
	DEALLOCATE (OLDT)
	DEALLOCATE (RLOG)
	DEALLOCATE (TA)
	DEALLOCATE (OLDR)
C
	RETURN
	END
