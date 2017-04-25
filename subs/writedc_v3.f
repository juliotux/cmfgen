	SUBROUTINE WRITEDC_V3(HYD,LOG_HYDLTE,NHYD,DHYD,NION,R,T,ED,V,CLUMP_FAC,
	1                         LUM,ND,FILENAME,OPTION,FORM)
	IMPLICIT NONE
!
! Altered 06-Sep-2016 : Now output 10 digits for R.
! Altered 21-Sep-2011 : Only outputs LOG(DC) when minimum dc < 10^{-290}
! Altered 05-Apr-2011 : Based on WRITEDC_V3 (10-Dec-2010)
!                         LOG_HYDLTE instead of HYDLTE passed in call.
!                         Now ouput Log(DCs) instead of DCs.
! Altered 04-Oct-2004 : Output space after depth index so easier to search for.
! Altered 24-Feb-2004 : Now depth index to first line of output. Should not effect
!                         any input files.
! Altered 07-Jul-1997 : CLUMP_FAC inserted in call (now _V2), and now output
!                         as last argument. R now written out with a precision
!                         of 7 decimal digits.
! Altered 26-Jun-1996 : CALL GEN_ASCI_OPEN installed.
! Altered 28-May-1996 : Removed for [jdh.disp]SETVEC routine
!                       DOUBLE PRECISION declaration removed.
! Altered  4-Aug-1988 : Write Departure coefficients out - not b-1.
!
	INTEGER NHYD,NION,ND,FORM,I,J,IOS
	INTEGER LU
	INTEGER, PARAMETER :: IZERO=0
	REAL*8 HYD(NHYD,ND),LOG_HYDLTE(NHYD,ND),DHYD(NION,ND)
	REAL*8 R(ND),T(ND),ED(ND),V(ND),CLUMP_FAC(ND)
	REAL*8 LUM,T1,T2
	REAL*8 LOG_TEN
	LOGICAL LOG_OUTPUT
	CHARACTER*(*)FILENAME,OPTION
	CHARACTER*30 NEWNAME
	CHARACTER*90 FMT
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
	NEWNAME=FILENAME
	LOG_TEN=LOG(10.0D0)
	LU=9
!
! Most of these formats are no longer in use.
!
! 1 = H, HeII
! 2 = HeI Sing
! 3 = HeI triplets, CIV,NV
! 4 = CIII,NIV
! 5 = HeI (Singlets and Triplets)
!
	FMT='(1X,5ES16.7)'		!Default format
	IF(FORM .EQ. 2)FMT='(1X,1P,1E15.5,:/,1X,2E15.5,:/,1X,3E15.5,:/,(1X,4E15.5))'
	IF(FORM .EQ. 3)FMT='(1X,1P,2E15.5,:/,1X,3E15.5,:/,(1X,4E15.5))'
	IF(FORM .EQ. 4)FMT='(1X,1P,1E15.5,:/,1X,2E15.5,:/,1X,3E15.5,:/,(1X,6E15.5))'
	IF(FORM .EQ. 5)FMT='(1X,1P,1E15.5,:/,1X,2E15.5,:/,1X,2E15.5,:/,1X,
	1                           3E15.5,:/,1X,3E15.5,/:,(1X,6E15.5))'
!
!
	IF(DHYD(1,ND) .NE. 0)THEN
!
	  CALL GEN_ASCI_OPEN(LU,NEWNAME,'REPLACE',' ',' ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error opening D.C, file',NEWNAME
	    WRITE(LUER,*)'IOSTAT=',IOS
	    RETURN
	  END IF
!
	  T1=1.0D0
	  DO I=1,ND
	    DO J=1,NHYD
	      T2=LOG(HYD(J,I))-LOG_HYDLTE(J,I)
	      T1=MIN(T1,T2)
	    END DO
	  END DO
	  T1=T1/LOG_TEN
	  IF(T1 .LT. -290.0D0)THEN
	    LOG_OUTPUT=.TRUE.
	    WRITE(LU,'(/,1X,A,T40,A)')'10-Dec-2010','!Format date'
	  ELSE
	    LOG_OUTPUT=.FALSE.
	    WRITE(LU,'(/,1X,A,T40,A)')'24-FEB-2004','!Format date'
	  END IF
!
	  WRITE(LU,2120)R(ND),LUM,NHYD,ND
	  IF(OPTION(1:2) .EQ. 'DC' .AND. LOG_OUTPUT)THEN
	    DO I=1,ND
	      T1=0.0D0
	      T2=0.0D0
	      DO J=1,NHYD
	        T1=T1+HYD(J,I)
	      END DO
	      DO J=1,NION
	        T2=T2+DHYD(J,I)
	      END DO
	      T1=T1/T2
	      WRITE(LU,2122)R(I),DHYD(1,I),ED(I),T(I),T1,V(I),CLUMP_FAC(I),I,' '
	      WRITE(LU,FMT)((LOG(HYD(J,I))-LOG_HYDLTE(J,I))/LOG_TEN,J=1,NHYD)
	    END DO
	  ELSE IF(OPTION(1:2) .EQ. 'DC')THEN
	    DO I=1,ND
	      T1=0.0D0
	      T2=0.0D0
	      DO J=1,NHYD
	        T1=T1+HYD(J,I)
	      END DO
	      DO J=1,NION
	        T2=T2+DHYD(J,I)
	      END DO
	      T1=T1/T2
	      WRITE(LU,2122)R(I),DHYD(1,I),ED(I),T(I),T1,V(I),CLUMP_FAC(I),I,' '
	      WRITE(LU,FMT)(EXP(LOG(HYD(J,I))-LOG_HYDLTE(J,I)),J=1,NHYD)
	    END DO
	  ELSE
	    DO I=1,ND
	      T1=0.0D0
	      T2=0.0D0
	      DO J=1,NHYD
	        T1=T1+HYD(J,I)
	      END DO
	      DO J=1,NION
	        T2=T2+DHYD(J,I)
	      END DO
	      T1=T1/T2
	      WRITE(LU,2122)R(I),DHYD(1,I),ED(I),T(I),T1,V(I),CLUMP_FAC(I),I,' '
	      WRITE(LU,FMT)(HYD(J,I),J=1,NHYD)
	    END DO
	  END IF
	  CLOSE(UNIT=LU)
	 END IF
!
2120	 FORMAT(/,ES17.10,4X,1PE11.4,5X,0P,I4,5X,I4)
2122	 FORMAT(/,ES17.10,6ES16.7,2X,I4,A1)
!
	 RETURN
	 END
