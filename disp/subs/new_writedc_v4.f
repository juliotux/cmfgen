C
	SUBROUTINE NEW_WRITEDC_V4(HYD,HYDLTE,WHYD,
	1               EDGEHYD,GHYD,NHYD,
	1               DHYD,GION,NION,R,T,ED,V,CLUMP_FAC,
	1               DO_DPTH,LUM,ND,FILENAME,OPTION,FORM)
	IMPLICIT NONE
C
C 11-Jan-2006 - Improved output format for R (useful for SN).
C 20-jan-2004 - Changed so that inner boundary radius can also be omitted.
C 28-May-2003 - Minor bug fix: Now correctly omits depths for TX and POP options.
C 20-Jul-2002 - ' qutes installed to avoid breaking strings across lines (FORM).
C 19-Jun-2000 - DO_DPTH inserted (changed to V4)
C               TEXCITE calculation improved.
C 05-Mar-1999 - TEXCITE now implicitly dimensioned by NHYD, not 200
C               WHYD (Level dissolution coeffiecients) inclued in call for
C                 forrect evaluation of TX.
C 07-Jul-1997 - CLUMP_FAC installed in call, and called _V2.
C                 R now written with 7 digits of precision.
C 17-Nov-1989 -Based on WRITEDC : TX option installed.
C
	INTEGER NHYD,NION,ND,FORM
	REAL*8 HYD(NHYD,ND),HYDLTE(NHYD,ND),WHYD(NHYD,ND)
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
	REAL*8 RCORE
	REAL*8 T1,T2,T3,CONST,DELTA_T
	REAL*8 TEXCITE(NHYD)
	INTEGER I,J,COUNT,ND_CNT
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
C
	NEWNAME=TRIM(FILENAME)
C
	  IF(DHYD(1,ND) .NE. 0)THEN
C
C 1 = H, HeII
C 2 = HeI Sing
C 3 = HeI triplets, CIV,NV
C 4 = CIII,NIV
C 5 = HeI (Singlets and Triplets)
C
	    IF(FORM .EQ. 1)FMT='(1X,1P,5E17.7)'
	    IF(FORM .EQ. 2)FMT='(1X,1P,1E15.5,:/,1X,2E15.5,:,/,1X,3E15.5,'//
	1                        ':,/,(1X,4E15.5))'
	    IF(FORM .EQ. 3)FMT='(1X,1P,2E15.5,:,/,1X,3E15.5,:,/,(1X,4E15.5))'
	    IF(FORM .EQ. 4)FMT='(1X,1P,1E15.5,:,/,1X,2E15.5,:,/,1X,3E15.5,'//
	1                       ':,/,(1X,6E15.5))'
	    IF(FORM .EQ. 5)FMT='(1X,1P,1E15.5,:,/,1X,2E15.5,:,/,1X,2E15.5,:/,1X,'//
	1            '3E15.5,:,/,1X,3E15.5,/:,(1X,6E15.5))'
C
	    OPEN(UNIT=9,STATUS='NEW',FILE=NEWNAME)
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
	    WRITE(9,'(/,1X,A,T40,A)')'07-Jul-1997','!Format date'
	    WRITE(9,2120)RCORE,LUM,NHYD,ND_CNT
	    IF(OPTION .EQ. 'DC')THEN
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
	          WRITE(9,2122)R(I),DHYD(1,I),ED(I),T(I),T1,V(I),CLUMP_FAC(I),I
	          WRITE(9,FMT)((HYD(J,I)/HYDLTE(J,I)),J=1,NHYD)
	        END IF
	      END DO
	    ELSE IF(OPTION .EQ. 'TX')THEN
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
	          WRITE(9,2122)R(I),DHYD(1,I),ED(I),T(I),T1,V(I),CLUMP_FAC(I),I
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
	          WRITE(9,FMT)(TEXCITE(J),J=1,NHYD)
	        END IF
	      END DO
	    ELSE
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
	          WRITE(9,2122)R(I),DHYD(1,I),ED(I),T(I),T1,V(I),CLUMP_FAC(I),I
	          WRITE(9,FMT)(HYD(J,I),J=1,NHYD)
	        END IF
	      END DO
	    END IF
	    CLOSE(UNIT=9)
	  END IF
C
2120	  FORMAT(/,1X,ES16.8,5X,ES12.6,5X,I4,5X,I4)
2122	  FORMAT(/,1X,ES16.8,6ES17.8,3X,I4)
C
	  RETURN
	  END
