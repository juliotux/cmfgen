!
! Auxilary program designed to modify the HYDRO file output from CMFGEN.
! Progam modifies the adopted stellar mass. Not that the percentage error
! os now defined so that it has a rang of pm 200%.
!
	PROGRAM PLT_GF
	USE GEN_IN_INTERFACE
!
! Cleaned: 07-Nov-200
!
	IMPLICIT NONE
!
	INTEGER I,J,IOS
	INTEGER NSTR
	INTEGER ND
	INTEGER, PARAMETER :: LU_OUT=11
	INTEGER, PARAMETER :: T_OUT=6
!
	INTEGER, PARAMETER :: ND_MAX=200
	REAL*8 V(ND_MAX)
	REAL*8 INDX(ND_MAX)
	REAL*8 FRAC(ND_MAX,12)
	REAL*8 LAM(12)
!
	CHARACTER*132 STRING
	CHARACTER*132 FMT
	CHARACTER*132 FILENAME
!
	FILENAME='HYDRO'
	CALL GEN_IN(FILENAME,'Input hydro file')
	OPEN(UNIT=10,FILE='HYDRO',ACTION='READ',STATUS='OLD')
	  NSTR=0        
	  DO WHILE(1 .EQ. 1)
	    READ(10,'(A)',END=1000)STRING
	    IF(STRING .EQ. ' ' .AND. ND .EQ. 0)THEN
	      ND=NSTR-1
	      EXIT
	    END IF
	    NSTR=NSTR+1
	  END DO         
1000	CONTINUE
!
	WRITE(T_OUT,*)'Number of depth points is',ND
        I=160
!
	DO WHILE(1 .EQ. 1)
	   READ(10,'(A)',END=2000)STRING
	   IF( INDEX(STRING,'Stars mass') .NE. 0)EXIT
	END DO         
	DO WHILE(1 .EQ. 1)
	   READ(10,'(A)',END=2000)STRING
	   IF( INDEX(STRING,' V ') .NE. 0)EXIT
	END DO
2000	CONTINUE
!         
	DO I=1,ND
	  READ(10,*)INDX(I),V(I),(FRAC(I,J),J=1,12)
	   READ(10,'(A)',END=3000)STRING
	   READ(10,'(A)',END=3000)STRING
	END DO
3000	CONTINUE
	FRAC(1:ND,1:12)=FRAC(1:ND,1:12)*100.0D0
!
	DO J=1,12
	  CALL DP_CURVE(ND,V,FRAC(1,J))
	END DO
	CALL GRAMON_PGPLOT('V(km/s)','%Force',' ',' ')
!
	STOP
	END
