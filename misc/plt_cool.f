	PROGRAM PLT_NRR
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered: 01-Mar-2016 : Use MODEL etc to get ND [24-Feb-2016].
!                          Changed reading -- string is 'Luminosity Check' rather than
!                               'Total (Rad. + Mech.'. May need updating.
!
	INTEGER, PARAMETER :: ND_MAX=200
	INTEGER, PARAMETER :: LU_RD=7
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
	INTEGER IOS
	INTEGER I,IBEG
	INTEGER ND
	INTEGER LUM_STAR 
!
	WRITE(6,'(A)')' '
        ND=0
	OPEN(UNIT=LU_RD,FILE='MODEL',STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    DO WHILE(1 .EQ. 1)
	      READ(LU_RD,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .NE. 0)EXIT
	      IF(INDEX(STRING,'!Number of depth points') .NE. 0)THEN
	        READ(STRING,*)ND
	        WRITE(6,'(A,I4)')' Number of depth points in the model is:',ND
	      END IF
	    END DO
	    CLOSE(LU_RD)
	  END IF
	IF(ND .EQ. 0)CALL GEN_IN(ND,'Number of depth points in model')
!
	DO WHILE(1 .EQ. 1)
 	  FILENAME='GENCOOL'
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
	  CALL DP_CURVE(ND,R,NET_PER)
	END DO
!
1000	CONTINUE
!
	CALL GEN_IN(LUM_STAR,'Stellar luminosity')

	  OPEN(UNIT=11,FILE='OBSFLUX',STATUS='OLD',ACTION='READ')
!
	      DO WHILE (1 .EQ. 1)
	        READ(11,'(A)')STRING
	        IF(INDEX(STRING,'Luminosity Check') .NE. 0)EXIT
!	        IF(INDEX(STRING,'Total (Rad. + Mech.)') .NE. 0)EXIT
	      END DO
	      READ(11,*)(LUM(I),I=1,ND)
	  CLOSE(UNIT=11)
!
	  LUM(1:ND)=100.0D0*(LUM(1:ND)-LUM_STAR)/LUM_STAR
	  CALL DP_CURVE(ND,R,LUM)
!
	CALL GEN_IN(LUM_STAR,'Stellar luminosity')

!	CALL GRAMON_PGPLOT('Log(Ne)','Net % cooling rate',' ',' ')
	CALL GRAMON_PGPLOT('R','Net % cooling rate',' ',' ')
!
	STOP
	END
