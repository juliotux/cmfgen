	PROGRAM PLT_NRR
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
	INTEGER, PARAMETER :: ND_MAX=100
!
	REAL*8 R(ND_MAX)
	REAL*8 T(ND_MAX)
	REAL*8 ED(ND_MAX)
	REAL*8 LUM(ND_MAX)
!
	REAL*8 NET(ND_MAX)
	REAL*8 NET_PER(ND_MAX)
!
	CHARACTER*80 FILENAME
	CHARACTER*80 STRING
!
	INTEGER I,IBEG
	INTEGER ND
	INTEGER LUM_STAR 
!
	CALL GEN_IN(ND,'Number of depth points in model')
!
	DO WHILE(1 .EQ. 1)
 	  FILENAME=' '
	  CALL GEN_IN(FILENAME,'File with data to be plotted')
	  IF(FILENAME .EQ. ' ')GOTO 1000
!
	  OPEN(UNIT=11,FILE=FILENAME,STATUS='OLD',ACTION='READ')
!
	    DO IBEG=1,ND,10
	      READ(11,'(A)')STRING
	      READ(11,'(A)')STRING
	      READ(11,*)(R(I),I=IBEG,MIN(IBEG+9,ND))
!
	      READ(11,'(A)')STRING
	      READ(11,'(A)')STRING
	      READ(11,*)(T(I),I=IBEG,MIN(IBEG+9,ND))
!
	      READ(11,'(A)')STRING
	      READ(11,'(A)')STRING
	      READ(11,*)(ED(I),I=IBEG,MIN(IBEG+9,ND))
!
	      DO WHILE (1 .EQ. 1)
	        READ(11,'(A)')STRING
	        IF(INDEX(STRING,'Net Cooling') .NE. 0)EXIT
	      END DO
	     
!
	      READ(11,*)(NET(I),I=IBEG,MIN(IBEG+9,ND))
!
	      READ(11,'(A)')STRING
	      READ(11,'(A)')STRING
	      READ(11,*)(NET_PER(I),I=IBEG,MIN(IBEG+9,ND))
	      IF(IBEG+9 .LT. ND)READ(11,'(A)')STRING
	    END DO
	  CLOSE(UNIT=11)
!
	  ED(1:ND)=DLOG10(ED(1:ND))
	  CALL DP_CURVE(ND,ED,NET_PER)
	END DO
!
1000	CONTINUE
!
	CALL GEN_IN(LUM_STAR,'Stellar luminosity')

	  OPEN(UNIT=11,FILE='OBSFLUX',STATUS='OLD',ACTION='READ')
!
	      DO WHILE (1 .EQ. 1)
	        READ(11,'(A)')STRING
	        IF(INDEX(STRING,'Total (Rad. + Mech.)') .NE. 0)EXIT
	      END DO
	      READ(11,*)(LUM(I),I=1,ND)
	  CLOSE(UNIT=11)
!
	  LUM(1:ND)=100.0D0*(LUM(1:ND)-LUM_STAR)/LUM_STAR
	  CALL DP_CURVE(ND,ED,LUM)
!
	CALL GEN_IN(LUM_STAR,'Stellar luminosity')

	CALL GRAMON_PGPLOT('Log(Ne)','Net % cooling rate',' ',' ')
!
	STOP
	END