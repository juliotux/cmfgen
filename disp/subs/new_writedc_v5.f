!
	SUBROUTINE NEW_WRITEDC_V5(HYD,LOG_HYDLTE,WHYD,
	1               EDGEHYD,GHYD,NHYD,
	1               DHYD,GION,NION,R,T,ED,V,CLUMP_FAC,
	1               DO_DPTH,LUM,ND,FILENAME,OPTION,FORM)
	IMPLICIT NONE
!
! 25-Sep-2011 - Changed to output Log10(DC's) when DCs are very small.
!                 Based on earlier changes which always output log(DCs).
! 11-Jan-2006 - Improved output format for R (useful for SN).
! 20-jan-2004 - Changed so that inner boundary radius can also be omitted.
! 28-May-2003 - Minor bug fix: Now correctly omits depths for TX and POP options.
! 20-Jul-2002 - '' quotes installed to avoid breaking strings across lines (FORM).
! 19-Jun-2000 - DO_DPTH inserted (changed to V4)
!               TEXCITE calculation improved.
! 05-Mar-1999 - TEXCITE now implicitly dimensioned by NHYD, not 200
!               WHYD (Level dissolution coeffiecients) inclued in call for
!                 forrect evaluation of TX.
! 07-Jul-1997 - CLUMP_FAC installed in call, and called _V2.
!                 R now written with 7 digits of precision.
! 17-Nov-1989 -Based on WRITEDC : TX option installed.
!
	INTEGER NHYD,NION,ND,FORM
	REAL*8 HYD(NHYD,ND),LOG_HYDLTE(NHYD,ND),WHYD(NHYD,ND)
	REAL*8 EDGEHYD(NHYD),GHYD(NHYD)
	REAL*8 DHYD(NION,ND),R(ND),T(ND),ED(ND),V(ND),CLUMP_FAC(ND)
	REAL*8 LUM
	REAL*8 GION
!
! DPTH_IND indicates which depths are to be output. Useful for
! debugging purposes when running a new model.
!
	LOGICAL DO_DPTH(ND)
	CHARACTER*(*)FILENAME,OPTION
	CHARACTER*30 NEWNAME
	CHARACTER*90 FMT
!
	REAL*8 LOG_TEN
	REAL*8 RCORE
	REAL*8 T1,T2,T3,CONST,DELTA_T
	REAL*8 TEXCITE(NHYD)
	INTEGER I,J,COUNT,ND_CNT
	LOGICAL LOG_OUTPUT
	INTEGER, PARAMETER :: LU_OUT=9
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
!
	NEWNAME=TRIM(FILENAME)
!
	LOG_TEN=LOG(10.0D0)
!
	  IF(DHYD(1,ND) .NE. 0)THEN
!
! 1 = H, HeII
! 2 = HeI Sing
! 3 = HeI triplets, CIV,NV
! 4 = CIII,NIV
! 5 = HeI (Singlets and Triplets)
!
	    IF(FORM .EQ. 1)FMT='(1X,1P,5E17.7)'
	    IF(FORM .EQ. 2)FMT='(1X,1P,1E15.5,:/,1X,2E15.5,:,/,1X,3E15.5,'//
	1                        ':,/,(1X,4E15.5))'
	    IF(FORM .EQ. 3)FMT='(1X,1P,2E15.5,:,/,1X,3E15.5,:,/,(1X,4E15.5))'
	    IF(FORM .EQ. 4)FMT='(1X,1P,1E15.5,:,/,1X,2E15.5,:,/,1X,3E15.5,'//
	1                       ':,/,(1X,6E15.5))'
	    IF(FORM .EQ. 5)FMT='(1X,1P,1E15.5,:,/,1X,2E15.5,:,/,1X,2E15.5,:/,1X,'//
	1            '3E15.5,:,/,1X,3E15.5,/:,(1X,6E15.5))'
!
	    OPEN(UNIT=LU_OUT,STATUS='NEW',FILE=NEWNAME)
!
! Determine number points that wll be output.
!
	    ND_CNT=0
	    DO I=1,ND
	      IF(DO_DPTH(I))THEN
	        ND_CNT=ND_CNT+1
	        RCORE=R(I)
	      END IF
	    END DO
!
	    IF(OPTION .EQ. 'DC')THEN
!
! Determine if we will output DC's in LOG format.
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
	        WRITE(LU_OUT,'(/,1X,A,T40,A)')'10-Dec-2010','!Format date'
	      ELSE
	        LOG_OUTPUT=.FALSE.
	        WRITE(LU_OUT,'(/,1X,A,T40,A)')'24-FEB-2004','!Format date'
	      END IF
	      WRITE(LU_OUT,2120)RCORE,LUM,NHYD,ND_CNT
!
	      DO I=1,ND
	        IF(DO_DPTH(I))THEN
	          T1=0.0D0
	          T2=0.0D0
	          DO J=1,NHYD
	            T1=T1+HYD(J,I)
	          END DO
	          DO J=1,NION
	            T2=T2+DHYD(J,I)
	          END DO
	          T1=T1/T2
	          WRITE(LU_OUT,2122)R(I),DHYD(1,I),ED(I),T(I),T1,V(I),CLUMP_FAC(I),I
	          IF(LOG_OUTPUT)THEN
	             WRITE(LU_OUT,FMT)( (LOG(HYD(J,I))-LOG_HYDLTE(J,I))/LOG_TEN,J=1,NHYD )
	          ELSE
	             WRITE(LU_OUT,FMT)( EXP(LOG(HYD(J,I))-LOG_HYDLTE(J,I)),J=1,NHYD )
	          END IF
	        END IF
	      END DO
	    ELSE IF(OPTION .EQ. 'TX')THEN
	      WRITE(LU_OUT,'(/,1X,A,T40,A)')'07-Jul_1997','!Format date'
	      WRITE(LU_OUT,2120)RCORE,LUM,NHYD,ND_CNT
	      DO I=1,ND
	        IF(DO_DPTH(I))THEN
	          T1=0.0D0
	          T2=0.0D0
	          DO J=1,NHYD
	            T1=T1+HYD(J,I)
	          END DO
	          DO J=1,NION
	            T2=T2+DHYD(J,I)
	          END DO
	          T1=T1/T2
	          WRITE(LU_OUT,2122)R(I),DHYD(1,I),ED(I),T(I),T1,V(I),CLUMP_FAC(I),I
	          DO J=1,NHYD
	            T2=HDKT*EDGEHYD(J)
	            CONST=HYD(J,I)*(T2**1.5)*GION/GHYD(J)/WHYD(J,I)/2.07078D-22/ED(I)/DHYD(1,I)
	            IF(CONST .LT. 2.8)THEN
	               TEXCITE(J)=CONST**(0.67)
	            ELSE
	               TEXCITE(J)=LOG(CONST)
	            END IF
	            COUNT=0
	            DELTA_T=1.0D+10
	            DO WHILE( ABS(DELTA_T/TEXCITE(J)) .GT. 1.0E-08 .AND. COUNT .LT. 100 )
	              COUNT=COUNT+1
	              T1=SQRT(TEXCITE(J))
	              T3=EXP(TEXCITE(J))
	              DELTA_T=(T1*TEXCITE(J)*T3 - CONST)/T3/T1/(1.5D0+TEXCITE(J))
	              IF(DELTA_T .GT. 0.8*TEXCITE(J))DELTA_T=0.8*TEXCITE(J)
	              IF(DELTA_T .LT. -0.8*TEXCITE(J))DELTA_T=-0.8*TEXCITE(J)
	              TEXCITE(J)=TEXCITE(J)-DELTA_T
	              IF(I .EQ. ND .AND. J .EQ. 1)THEN
	                 WRITE(101,'(6ES14.4)')T2,CONST,TEXCITE(J),T1,T2,DELTA_T
	              END IF
	            END DO
	            TEXCITE(J)=T2/TEXCITE(J)
	            IF(COUNT .EQ. 100)THEN
	              WRITE(6,*)'Error - TEXC didnt converge in 100 iterations'
	              WRITE(6,*)'I,J=',I,J
	              RETURN
	            END IF
	          END DO
	          WRITE(LU_OUT,FMT)(TEXCITE(J),J=1,NHYD)
	        END IF
	      END DO
	    ELSE
	      WRITE(LU_OUT,'(/,1X,A,T40,A)')'07-Jul-1997','!Format date'
	      WRITE(LU_OUT,2120)RCORE,LUM,NHYD,ND_CNT
	      DO I=1,ND
	        IF(DO_DPTH(I))THEN
	          T1=0.0D0
	          T2=0.0D0
	          DO J=1,NHYD
	            T1=T1+HYD(J,I)
	          END DO
	          DO J=1,NION
	            T2=T2+DHYD(J,I)
	          END DO
	          T1=T1/T2
	          WRITE(LU_OUT,2122)R(I),DHYD(1,I),ED(I),T(I),T1,V(I),CLUMP_FAC(I),I
	          WRITE(LU_OUT,FMT)(HYD(J,I),J=1,NHYD)
	        END IF
	      END DO
	    END IF
	    CLOSE(UNIT=9)
	  END IF
!
2120	  FORMAT(/,1X,ES16.8,5X,ES12.6,5X,I4,5X,I4)
2122	  FORMAT(/,1X,ES16.8,6ES17.8,3X,I4)
!
	  RETURN
	  END
