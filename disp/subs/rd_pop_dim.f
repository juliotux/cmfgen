C
C Routine to read in populations for a particular ionization stage (eg C2) of
C a given species (eg C). The ion population (eg DC2). oscilator strengths,
C photoionization cross sections, and the LTE populations are also returned.
C
C The species in the file must be ordered from lowest to highest ioization.
C
C NB. CIII is used as the dummy species.
C
	SUBROUTINE RD_POP_DIM(NCIII,CIII_PRES,DESC,FORMAT_DATE,LUIN)
	IMPLICIT NONE
C
C Altered 23-Jan-2013 : Changed location of "Number of levels" warning to overcome
C                         an issue with GFORTRAN which ignored the "END=" statement.
C Altered 04-Jun-1998 : Routine checks for occurrence of DESC and upper case
C                         version of DESC. To handle change from HE2 to He2
C                         and HEI to HeI.
C Altered 27-Jan-1996 : Number of lvels descriptor no longer has to be present
C                         for a species not included in model.
C
	INTEGER NCIII
	INTEGER LUIN
	LOGICAL CIII_PRES
	CHARACTER*(*) DESC,FORMAT_DATE
	CHARACTER*30 UC
	EXTERNAL UC
C
C Local variables.
C
	CHARACTER*132 STRING
	INTEGER IOS,LUER,ERROR_LU,INDx
	EXTERNAL ERROR_LU
C
	LUER=ERROR_LU()
C
	IF(FORMAT_DATE .NE. '27-JAN-1992' .AND.
	1         FORMAT_DATE .NE. '30-NOV-1991')THEN
	  WRITE(LUER,*)'Error : invalid FORMAT DATE for '//DESC
	  WRITE(LUER,*)'Error occured in RD_POP_DIM'
	  STOP
	END IF
C
C Get number of CIII levels. If zero, we exit, otherwise we input the
C population levels.
C
	INDX=0
	NCIII=0
	DO WHILE(INDX .EQ. 0)
	  READ(LUIN,'(A)',IOSTAT=IOS,END=1000)STRING
C
C The test against 5001 is to overcome an issue with GFROTRAN which does not detect the
C end of the file. This should cause no issues with othr compilers.
C
	  IF(IS_IOSTAT_END(IOS) .OR. IOS .EQ. 5001)THEN
	    GOTO 1000
	  ELSE IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'***************************************************'
	    WRITE(LUER,*)'***************************************************'
	    WRITE(LUER,*)'Warning : Unable to get Number of levels record'
	    WRITE(LUER,*)'Species is: ', DESC
	    WRITE(LUER,*)'Warning ocurred in RD_POP_DIM'
	    WRITE(LUER,*)'IOS=',IOS
	    WRITE(LUER,*)'***************************************************'
	    WRITE(LUER,*)'***************************************************'
	    CIII_PRES=.FALSE.
	    NCIII=1		!So NCIII can be used as dimension limit.
	    RETURN
	  END IF
!
	  INDX=INDEX(STRING,'Number of ')
	END DO
C
	IF( INDEX(STRING,'Number of '//DESC//' levels:') .NE. 0 .OR.
	1       INDEX(STRING,'Number of '//TRIM(UC(DESC))//' levels:')
	1                                     .NE. 0 )THEN
	  READ(STRING,'(T30,BN,I8)')NCIII
	  CIII_PRES=.TRUE.
	ELSE
C
C We backspace so that the record with the numer of levels is available for
C the next call to RD_POP_DIM.
C
	  BACKSPACE(LUIN)
!	  WRITE(LUER,*)'  ',DESC//' data unavaialable'
	  CIII_PRES=.FALSE.
	  NCIII=1		!So NCIII can be used as dimension limit.
	END IF
C
1000	CONTINUE
	IF(NCIII .EQ. 0)THEN
!	  WRITE(LUER,*)'  ',DESC//' data unavaialable'
	  CIII_PRES=.FALSE.
	  NCIII=1		!So NCIII can be used as dimension limit.
	END IF
C
	RETURN
	END
 
