!
! Program reads in population estimates from a file and uses
! linear interpolation in the log plane to lay these estimates
! on the "new" radius grid. 
!
	SUBROUTINE REGRIDWSC_V3(DHEN,R,ED,T,DI,EDGE,F_TO_S,INT_SEQ,
	1                             POPATOM,N,ND,FILNAME)
	IMPLICIT NONE
!
! Altered 16-Oct-2012 - Added BA_TX_CONV and NO_TX_CONV
! Altered 26-Jan-2009 - Section involving TX fixeda
!                         ADD1 variable removed (correction done on spot).
! Altered 23-Nov-2007 - Procedure for extending populations to larger radii changes.
!                         We now fiddle departure coefficient to try to ensure
!                         constant ionization (actual alteration date 14-Nov-2007). 
! Altered 20-Sep-2007 - Bug fix. T not being correctly sent when outer boundary 
!                         extends beyond grid boundary.
! Altered 20-Sep-2006 - Changed so that Ne and DI are scaled by the atom
!                         density when extrapolating at the outer boundary.
!                         T is held fixed.
! Altered 23-Jun-2005 - Fixed computation of the excitation temperature.
! Altered 20-Feb-2005 - CHECK_DC installed.
! Altered 04-Jul-2004 - Altered so that levels in the same super level
!                          all have the same departure coefficient.
! Altered 15-Feb-2002 - Altered so that we can input files where the 
!                          outer radius is greater than the outer radius of
!                          the new model.
! Altered 14-Jun-2000 - Error check when file opened.
! Altered 07-Jul-1997 - We check the file for a '!Format date' as this
!                          effect our free-format read.
! Altered 13-Nov-1996 - Bug fix. Typo in setting maximum allocation for
!                         vectors --- error occurs when NOLD is the
!                         largest value.
!
! Altered 26-May-1996 - Access to scratch block removed. Dynamic memory
!                         allocation used.
!                       ERROR_LU installed.
!                       Generic call for LOG and EXP.
!
! Altered 02-Aug-1994 - PROGRDESC change from INTEGER to REAL*8. Corrects
!                       alignment problem. PROGDESC is now check for
!                       corruption, even if too many iterations occur.
! Altered 3-May-1989 - NV and check installed.
!
! Altered 29-Nov-1989 - EDGE installed to enable departure coefficients for
!                        unknown levels to be computed assuming the same
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
	REAL*8 POPATOM(ND)
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
	REAL*8 RPOLD,TX,FX,DELTA_T,T1,T2
	LOGICAL BAD_TX_CONV
	LOGICAL NO_TX_CONV(ND)
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
!	IF(R(1)/OLDR(1)-1.0D0 .GT. 0.0001)THEN
!	  WRITE(LUER,*)'Error in REGRIDWSC_V3'
!	  WRITE(LUER,*)'Model boundary radius greater than that of input model'
!	  STOP
!	ELSE IF(R(2) .GT. OLDR(1))THEN
!	  WRITE(LUER,*)'Error in REGRIDWSC_V3'
!	  WRITE(LUER,*)'Boundary radius R(2) too large'
!	  STOP
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
! Interpolate the populations assuming that the departure
! coefficients remain constant.
!
! When extending to larger radii, we assume that the ionization
! state of the gas does not change. This supersedes the earlier assumption
! of a fixed departure coefficient. Excited populations are adjusted
! logarithmically. Code assumes ED is proportional to POP_ATOM, and
! hence may break down for SN models (with variable composition).
!
	DO J=1,NZ
	  DO I=1,NDOLD
	    OLDT(I)=LOG(DPOP(J,I))
	  END DO
	  CALL LINPOP(RLOG(NXST),TA(NXST),NX,OLDR,OLDT,NDOLD)
	  DO I=NXST,ND
	    DHEN(J,I)=EXP(TA(I))
	  END DO
	  DO I=1,NXST-1
	    T1=POPATOM(NXST)/POPATOM(I)
	    T2=MIN(1.0D0, ABS(LOG(DPOP(J,1)))/ABS(DLOG(DPOP(1,1))) )
	    DHEN(J,I)=DPOP(J,1)*(T1**T2)
	  END DO
	END DO
!
! Compute departure coefficients for N>NZ. New levels are set to have these
! same excitation temperature as the highest level.
!
	NO_TX_CONV=.FALSE.
	BAD_TX_CONV=.FALSE.
	IF(N .GT. NZ)THEN
	  TX=T(ND)
	  DO I=ND,1,-1
!
! We first compute the excitation temperature on level NZ.
!
	    DELTA_T=100
	    COUNT=0
	    DO WHILE(ABS(DELTA_T) .GT. 1.0D-08 .AND. COUNT .LT. 100)
	      FX=DHEN(NZ,I)*EXP(HDKT*EDGE(NZ)*(1.0D0/T(I)-1.0D0/TX))*
	1           (TX/T(I))**1.5D0
	      DELTA_T=(FX-1.0D0)*TX/FX/(1.5D0+HDKT*EDGE(NZ)/TX)
	      IF(DELTA_T .GT.  0.8D0*TX)DELTA_T=0.8D0*TX
	      IF(DELTA_T .LT. -0.8D0*TX)DELTA_T=-0.8D0*TX
	      TX=TX-DELTA_T
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
	        DHEN(J,I)=EXP(HDKT*EDGE(J)*(1.0D0/TX-1.0D0/T(I)))*(T(I)/TX)**1.5D0
	      END DO
	    END IF
	  END DO
!
	  IF(BAD_TX_CONV)THEN
	    WRITE(LUER,*)'Error in REGRIDWSC_V3 - TX did not converge in 100 iterations'
	    WRITE(LUER,*)'File is ',TRIM(FILNAME)
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
	IF(SUM(INT_SEQ) .EQ. 0.0D0)THEN
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
