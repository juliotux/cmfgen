!
! Program reads in population estimates from a file and uses
! linear interpolation in the log plane to lay these estimates
! on the "new" radius grid. 
!
	SUBROUTINE REGRIDWSC_V2(DHEN,R,ED,T,DI,EDGE,F_TO_S,INT_SEQ,
	1                             N,ND,FILNAME)
	IMPLICIT NONE
!
! Altered 23-Jun-2005 - Fixed computation of the excitation temperature.
! Altered 20-Feb-2005 - CHECK_DC installed.
! Altered 04-Jul-2004 - Altered so that levels in the same super level
!                          all have the same departure coeficient.
! Altered 15-Feb-2002 - Altered so that we can input files where the 
!                          outer radius ig reater than the outer radius of
!                          the new model.
! Altered 14-Jun-2000 - Error check when file opened.
! Altered 07-Jul-1997 - We check the file for a '!Format date' as this
!                          effect our free-format read.
! Altered 13-Nov-1996 - Bug fix. Typo in setting maximum allocation for
!                         vectors --- error occrrs when NOLD is the
!                         largest value.
!
! Altered 26-May-1996 - Access to scratch block removed. Dynamic memory
!                         allocation used.
!                       ERROR_LU installed.
!                       Generic call for OOG and EXP.
!
! Altered 02-Aug-1994 - PROGRDESC change from INTEGER to REAL*8. Corrects
!                       alignment problem. PROGDESC is now check for
!                       corruption, even if too many iterations occur.
! Altered 3-May-1989 - NV and check installed.
!
! Altered 29-Nov-1989 - EDGE installed to enable departure coefficients for
!                        unknow levels to be computed assuming the same
!                        excitation temperature as the last known level.
!                        NB - Call was changed. Implicit none installed.
! Altered 08-Oct-1987 - Handles departure coefficients, OR b-1.
! Altered 11-Mar-1987 - Reading sequence split into two calls so that V and ION
!                         may appear on the same line as R.
! Created 28-SEP-1985 - (Based on REGRID - scratch block added so the 
!                         NDOLD can be grater than ND)
! Altered 30-APR-1985 - (Bug fixed in error messages - import when ND .ne. NDOLD)
! Altered 28-FEB-1984
!
	INTEGER N,ND
	REAL*8 DHEN(N,ND),R(ND),T(ND),ED(ND),DI(ND),EDGE(N)
	INTEGER F_TO_S(N)
	INTEGER INT_SEQ(N)
!
	REAL*8, ALLOCATABLE :: DPOP(:,:)
	REAL*8, ALLOCATABLE :: OLDED(:)
	REAL*8, ALLOCATABLE :: OLDDI(:)
	REAL*8, ALLOCATABLE :: OLDT(:)
	REAL*8, ALLOCATABLE :: RLOG(:)
	REAL*8, ALLOCATABLE :: TA(:)
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
	REAL*8 RPOLD,TX,DELTA_T,T1,ADD1,DC_MOD
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
	ALLOCATE (DPOP(NOLD,NDOLD))
	ALLOCATE (OLDED(I))
	ALLOCATE (OLDDI(I))
	ALLOCATE (OLDT(I))
	ALLOCATE (RLOG(I))
	ALLOCATE (TA(I))
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
	ADD1=0.0D0
	IF( DABS( TA(NOLD) ) .LT. 0.2 .AND. CHECK_DC)ADD1=1.0D0
!
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
! Interpolate the populations assuming that the departure
! coefficients remain constant.
!
	DO J=1,NZ
	  DO I=1,NDOLD
	    OLDT(I)=LOG(DPOP(J,I)+ADD1)	!insures that never neg.
	  END DO
	  CALL LINPOP(RLOG(NXST),TA(NXST),NX,OLDR,OLDT,NDOLD)
	  DO I=NXST,ND
	    DHEN(J,I)=EXP(TA(I))
	  END DO
	END DO
!
! Compute departure coefficints for N>NZ . New levels are set to have thes
! same excitation temperature as the highest level/
!
	IF(N .GT. NZ)THEN
	  DO I=1,ND
!
! We first compute the excitation temperature on level NZ.
!
	    TX=T(I)
	    DELTA_T=10.0D0
	    COUNT=0
	    DC_MOD=LOG(DHEN(NZ,I)) + HDKT*EDGE(NZ)/T(I) - 1.5D0*LOG(T(I))
	    DO WHILE( ABS(DELTA_T) .GT. 1.0E-06 .AND. COUNT .LT. 100 )
	     T1=TX
             IF(HDKT*EDGE(NZ)/TX .GT. 1.5D0)THEN
               TX=HDKT*EDGE(NZ)/(DC_MOD+1.5D0*LOG(TX))
             ELSE
               TX= ( EXP(HDKT*EDGE(NZ)/TX-DC_MOD) )**(2.0D0/3.0D0)
             END IF
             COUNT=COUNT+1
             DELTA_T=ABS(TX-T1)
	    END DO
!
!	    DO WHILE( ABS(DELTA_T) .GT. 1.0E-06 .AND. COUNT .LT. 100 )
!	      T1=( (T(I)/TX)**(1.5) )*EXP(HDKT*EDGE(NZ)*(T(I)-TX)/T(I)/TX)
!	      COUNT=COUNT+1
!	      DELTA_T= (DHEN(NZ,I)-T1)*TX/T1/(1.0D0+HDKT*EDGE(NZ)/TX)
!	      IF(ABS(DELTA_T) .GT. 0.8*TX)DELTA_T=0.5D0*DELTA_T
!	      TX=TX-DELTA_T
!	    END DO
!
	    IF(COUNT .EQ. 100)THEN
	      WRITE(LUER,*)'Error in REGRIDWSC - TX didnt converge ',
	1               'in 100 iterations'
	      WRITE(LUER,*)'Depth=',I
	      DO J=NZ+1,N
	        DHEN(J,I)=DHEN(NZ,I)
	      END DO
	    ELSE
	      DO J=NZ+1,N
	         DHEN(J,I)=DHEN(NZ,I)*EXP( HDKT*(EDGE(NZ)-EDGE(J))*
	1                   (TX-T(I))/TX/T(I) )
	      END DO
	    END IF
	  END DO
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
! Free up memory.
!
	DEALLOCATE (DPOP)
	DEALLOCATE (OLDED)
	DEALLOCATE (OLDDI)
	DEALLOCATE (OLDT)
	DEALLOCATE (RLOG)
	DEALLOCATE (TA)
	DEALLOCATE (OLDR)
!
	RETURN
	END
