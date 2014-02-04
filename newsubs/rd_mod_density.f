!
! Routine to read in the density and clumping factor for a variable mass-loss
! rate model.
!
	SUBROUTINE RD_MOD_DENSITY(DENSITY,CLUMP_FAC,R,ND,DENSITY_FILE)
!
	IMPLICIT NONE
!
	INTEGER ND
	REAL*8 R(ND)
	REAL*8 DENSITY(ND)
	REAL*8 CLUMP_FAC(ND)
	CHARACTER*(*) DENSITY_FILE
!
! Local variables.
!
	REAL*8 T1,T2,T3
	INTEGER LUIN
	INTEGER LUER
	INTEGER I
	INTEGER IOS
	INTEGER ND_LOC
	INTEGER ERROR_LU
	INTEGER, PARAMETER :: IZERO=0
	CHARACTER*132 STRING
!
	LUER=ERROR_LU()
!
! R, V, and SIGMA in column format, with simple header.
! As output by NEWRG in DISPGEN
!
	LUIN=7
	CALL GEN_ASCI_OPEN(LUIN,DENSITY_FILE,'OLD',' ','READ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error in RD_MOD_CLUMP: IOSTAT=',IOS
	    WRITE(LUER,*)'Unable to open ',TRIM(DENSITY_FILE)
	    STOP
	  END IF
	  STRING=' '
	  DO WHILE( INDEX(STRING,'!Number of depth points').EQ. 0)
	    READ(LUIN,'(A)')STRING
	  END DO
	  READ(STRING,*)ND_LOC
	  IF(ND_LOC .NE. ND)THEN
	    WRITE(LUER,*)'Error in ',TRIM(DENSITY_FILE)
	    WRITE(LUER,*)'Routine can''t yet handle a differnet number of depth points'
	    STOP
	  END IF
!
! Skip any further blank strings or comments.
!
	  STRING=' '
	  DO WHILE(STRING .EQ. ' ' .OR. STRING(1:1) .EQ. '!')
	    READ(LUIN,'(A)')STRING
	  END DO
	  BACKSPACE(LUIN)
!
	  DO I=1,ND
	    READ(LUIN,*)T1,T2,T3,DENSITY(I),CLUMP_FAC(I)
	    IF( ABS(T1-R(I))/R(I) .GT. 1.0D-06)THEN
	      WRITE(LUER,*)'Error in RD_MOD_DENSITY'
	      WRITE(LUER,*)'R scales don''t agree'
	      WRITE(LUER,*)I,R(I),T1
	      STOP
	    END IF
	  END DO
	CLOSE(LUIN)
!
	RETURN
	END
