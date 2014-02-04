!
! Subroutine to rewrite the free-format SCRTEMP file. The user
! input how many of the final iterations are to be output.
! NEW_SCRTEMP, NEW_POINT1, and NEW_POINT1 are created. These need
! to be renamed. This allos space to be saved, and at the same time 
! allows a new model to be run with exactly the same input populations.
!
	PROGRAM REWRITE_SCR
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Revised: 29-May-2009 : Bug fixed -- didn't output correct record when IREC NE NITSF.
! Revised: 20-Jan-2006 : Updated so that R, V, & SIGMA writing option can be altered.
! Revised: 06-May-2004 : Updated reading of SCRTEMP to SCR_READ_V2
! Revised: 13-Feb-2002
!
	INTEGER ND,NT,NIT
!
	REAL*8, ALLOCATABLE :: POPS(:,:)		!NT,ND
	REAL*8, ALLOCATABLE :: R(:)			!ND
	REAL*8, ALLOCATABLE :: V(:)			!ND
	REAL*8, ALLOCATABLE :: SIGMA(:)			!ND
!
	INTEGER, PARAMETER :: T_OUT=6
!
	INTEGER NPLTS
	INTEGER IREC
	INTEGER IREC_RD
	INTEGER IREC_WR
	INTEGER NIT_WR
	INTEGER NITSF
	INTEGER LST_NG
	INTEGER WRITEN_N_TIMES
	INTEGER RITE_N_TIMES
	INTEGER LUSCR
	INTEGER IOS
	INTEGER I,K,J,L
!
	LOGICAL NEWMOD
	LOGICAL READ_WRITE_RVSIG	!Indicates if R, V & SIGMA were output for each iteration.
	LOGICAL WRITE_RVSIG             !As above, but for new file.
!
	CHARACTER*132 STRING
!
	LUSCR=26
	RITE_N_TIMES=1
!
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'This routine should be run from the data directory'
	WRITE(T_OUT,*)'It expects to find the following files:'
	WRITE(T_OUT,*)'                                       POINT1.DAT'
	WRITE(T_OUT,*)'                                       SCRTEMP.DAT'
	WRITE(T_OUT,*)'                                       MODEL.DAT'
	WRITE(T_OUT,*)' '
!
!
! Get NT (# of levels) and ND (number of depth points) from MODEL_SPEC.
!
	OPEN(UNIT=12,FILE='MODEL',STATUS='OLD',ACTION='READ')
	  STRING=' '
	  DO WHILE(INDEX(STRING,'!Number of depth') .EQ. 0)
	    READ(12,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)GOTO 100
	  END DO
	  READ(STRING,*)ND
	  DO WHILE(INDEX(STRING,'!Total number of variables') .EQ. 0)
	    READ(12,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)GOTO 100
	  END DO
	  READ(STRING,*)NT
!
100	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Unable to read MODEL file'
	    CALL GEN_IN(NT,'Total number of levels')
	    CALL GEN_IN(ND,'Number of depth points')
	  END IF 
	CLOSE(UNIT=12)
!
! Read record information from SCRTEMP information file POINT1.
!
	OPEN(UNIT=12,FILE='POINT1',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .EQ. 0)READ(12,'(A)',IOSTAT=IOS)STRING
	  IF(IOS .EQ. 0)THEN
	    IF(INDEX(STRING,'!Format date') .EQ. 0)THEN
	      READ(STRING,*,IOSTAT=IOS)IREC_RD,NIT,WRITEN_N_TIMES
	    ELSE
	      READ(12,*,IOSTAT=IOS)IREC_RD,NIT
	    END IF
	  END IF
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Possible error reading POINT1'
	    CALL GEN_IN(NIT,'Number of iterations')
	  END IF
	CLOSE(UNIT=12)
!
	ALLOCATE (POPS(NT,ND))
	ALLOCATE (R(ND))
	ALLOCATE (V(ND))
	ALLOCATE (SIGMA(ND))
!
!
! NIT_WR=2 means that we write out the LAST 2 iterations only.
!
	NIT_WR=2
	CALL GEN_IN(NIT_WR,'Number of (FINAL) iterations to write')
	IF(NIT_WR .GT. NIT)THEN
	   WRITE(6,*)'Error -- asking to write more iterations than available'
	   STOP
	END IF
!
	IREC_WR=0
	DO I=1,NIT_WR
	   IREC=IREC_RD-(NIT_WR-I)
	   CALL SCR_READ_V2(R,V,SIGMA,POPS,IREC,NITSF,
	1              WRITEN_N_TIMES,LST_NG,READ_WRITE_RVSIG,
	1              NT,ND,LUSCR,NEWMOD)
	   IF(NEWMOD)THEN
	      WRITE(6,'(A,I4,A)')' Unable to read iteration # ',IREC,' from SCRTEMP'
	      STOP
	   END IF
	   IF(I .EQ. 1)THEN
	     WRITE_RVSIG=READ_WRITE_RVSIG
	     CALL GEN_IN(WRITE_RVSIG,'Write R, V & SIGMA for each iteration')
	   END IF
!
	   LST_NG=-1000
	   CALL SCR_RITE_NAM_V2(R,V,SIGMA,POPS,IREC_WR,I,
	1              RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              NT,ND,LUSCR,NEWMOD)
!
	END DO
!
	STOP
	END
