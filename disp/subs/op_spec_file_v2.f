C
C Open population file. Get model parameters. Check that model and main times
C are compatible. Return an error message depending on the error.
C
	SUBROUTINE OP_SPEC_FILE_V2(FILNAME,LU,ABUND_SPEC,POP_VEC,ND,
	1            FORMAT_DATE,IOS,TIME,DESC)
	IMPLICIT NONE
C
	INTEGER LU,IOS,ND
	REAL*8 ABUND_SPEC
	REAL*8 POP_VEC(ND)
	CHARACTER*(*) DESC,FILNAME,TIME,FORMAT_DATE
C
C Local variables
C
	CHARACTER*132 STRING
	CHARACTER*20 TIMECHK
	INTEGER I
C
C New asci formatted files.
C
	OPEN(UNIT=LU,FILE=FILNAME,STATUS='OLD',IOSTAT=IOS,ACTION='READ')
	  READ(LU,'(A)')STRING
	  IF( INDEX(STRING,'Output format date:') .NE. 0)THEN
	    READ(STRING,'(T30,A11)')FORMAT_DATE
	  ELSE
	    WRITE(6,*)'Error getting formate date for '//DESC
	    IOS=2
	    RETURN
	  END IF
	  READ(LU,'(A)')STRING
	  IF( INDEX(STRING,'Completion of Model:') .NE. 0)THEN
	    READ(STRING,'(T30,A20)')TIMECHK
	  ELSE
	    WRITE(6,*)'Error getting TIMECHK for '//DESC
	    IOS=3
	    RETURN
	  END IF
	  READ(LU,'(A)')STRING
	  IF( INDEX(STRING,'ND:') .EQ. 0)THEN
	    WRITE(6,*)'Error getting # of depth points in '//DESC
	    IOS=4
	    RETURN
	  END IF
	  READ(LU,'(A)')STRING
	  IF( INDEX(STRING,'abundance:') .NE. 0)THEN
	    READ(STRING,'(T30,BN,F16.0)')ABUND_SPEC
C
C Read in population 
C
	    READ(LU,*)(POP_VEC(I),I=1,ND)
	  ELSE
	    WRITE(6,*)'Error getting formate date for '//DESC
	    IOS=5
	    RETURN
	  END IF
C
	IF(TIME .NE. TIMECHK)THEN
	  WRITE(6,860)DESC,TIME,DESC,TIMECHK
860	  FORMAT(' WARNING ',A,' and main model times unequal',
	1   /,1X,'MAIN TIME= ',A20,/X,A20,'Time=',A20)
	END IF
	WRITE(6,'(A,A,T31,A,1P,E14.6)')
	1            ' Fractional abundance of ',TRIM(DESC),' is',ABUND_SPEC
C
	RETURN
	END
