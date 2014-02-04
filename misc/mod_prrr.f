!
! Program to put the GENCOOL file into a more user friendly format.
! Two files are created:
!
!     GENCOOL_SUM:   Same as GENCOOL but we have summed up over all bound-free rates.
!     GENSCOOL_SORT: Only the top rates are printed: Sorted using depths 1, 11, 21 etc. 
!
	PROGRAM MOD_PRRR
	USE MOD_COLOR_PEN_DEF
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered: 17-Nov-2009: Now read in charge exchange cooling.
!                         Slight format change.
! Altered: 29-Jan-2009: ND is now read in from MODEL (if it exists).
! Altered: 08-Feb-2008: Extra terms (such as V term) sheck and output.
!
!
	INTEGER, PARAMETER :: ND_MAX=200
	REAL*8 PHOT(ND_MAX,2000)
	REAL*8 RECOM(ND_MAX)
	REAL*8 RECOM_SUM(ND_MAX)
	REAL*8 PHOT_SUM(ND_MAX)
	REAL*8 COL_IR(ND_MAX)
	REAL*8 CHG_IR(ND_MAX)
	REAL*8 NT_IR(ND_MAX)
	REAL*8 COL_RR(ND_MAX)
	REAL*8 CHG_RR(ND_MAX)
	REAL*8 ADVEC_RR(ND_MAX)
!
	REAL*8 R(ND_MAX)
	REAL*8 V(ND_MAX)
	REAL*8 TEMP(ND_MAX)
	REAL*8 ED(ND_MAX)
	REAL*8 DI(ND_MAX)
!
	REAL*8 XVEC(ND_MAX)
	REAL*8 YSUM(ND_MAX)
!
	INTEGER ND
	INTEGER NV
	INTEGER I,J,K,L
	INTEGER N_INIT_RECS
	INTEGER NRECS
	INTEGER IOS
	INTEGER IST,IEND,NLEV
	REAL*8 T1
	REAL*8 MIN_VAL
	LOGICAL FILE_OPEN
	LOGICAL NET_RECOM_PER_LEVEL
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
        CHARACTER(LEN=2), PARAMETER :: FORMFEED=' '//CHAR(12)
!
	CHARACTER(LEN=132) TMP_STR
	CHARACTER(LEN=132) STRING
	CHARACTER(LEN=132) FILE_NAME
	CHARACTER(LEN=10) SPECIES
	CHARACTER(LEN=2) XAX_OPT
	CHARACTER(LEN=30) XLABEL
	CHARACTER(LEN=30) YLABEL
        CHARACTER(LEN=30) UC
        EXTERNAL UC
!
	OPEN(UNIT=20,FILE='MODEL',STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    DO WHILE(1 .EQ. 1)
	      READ(20,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .NE. 0)EXIT
	      IF(INDEX(STRING,'!Number of depth points') .NE. 0)THEN
	        READ(STRING,*)ND
	        WRITE(6,'(A,I4)')' Number of depth points in the model is:',ND
	        EXIT
	      END IF
	    END DO
	  END IF
	  INQUIRE(UNIT=20,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=20)
!
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)' Unable to open MODEL file to get # of depth points'
	  CALL GEN_IN(ND,'Number of depth points')
	END IF
!
	IOS=0
	OPEN(UNIT=20,FILE='RVTJ',STATUS='OLD',IOSTAT=IOS)
	  DO WHILE(INDEX(STRING,'Velocity (km/s)') .EQ. 0)
            READ(20,'(A)',IOSTAT=IOS)STRING
	  END DO
	  READ(20,*,IOSTAT=IOS)(V(I),I=1,ND)
	CLOSE(UNIT=20)
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)'Unable to get V(km/s) from RVTJ'
	  V(1:ND)=1.0
	END IF
!
	NET_RECOM_PER_LEVEL=.FALSE.
	CALL GEN_IN(NET_RECOM_PER_LEVEL,'Ouput net recombination rate to each level?')
	SPECIES='FeI'
!
1000	CONTINUE
!
	CHG_IR(1:ND)=0.0D0
	COL_IR(1:ND)=0.0D0
	NT_IR(1:ND)=0.0D0
	CHG_RR(1:ND)=0.0D0
	COL_RR(1:ND)=0.0D0
	ADVEC_RR(1:ND)=0.0D0
	RECOM_SUM(1:ND)=0.0D0
	PHOT_SUM(1:ND)=0.0D0
!
	IOS=1
	DO WHILE(IOS .NE. 0)
	  CALL GEN_IN(SPECIES,'File is assumed to be SPECIES//PRRR- EX to exit')
	  IF(UC(SPECIES) .EQ. 'EX')STOP
	  FILE_NAME=TRIM(SPECIES)//'PRRR'
	  OPEN(UNIT=20,FILE=FILE_NAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,'(2A,I5,A)')' Error -- unable to open file: '//TRIM(FILE_NAME),' (IOS=',IOS,')'
	    WRITE(6,'(A)')' Check capatilzation of ION name'
	  END IF
	END DO
	FILE_NAME=TRIM(FILE_NAME)//'_SUM'
	OPEN(UNIT=21,FILE=FILE_NAME,STATUS='UNKNOWN',ACTION='WRITE')
!
	STRING=' '
	DO L=1,ND,10
	   IST=L; IEND=MIN(ND,IST+9)
	   WRITE(6,*)IST,IEND
!
	   DO WHILE(INDEX(STRING,'Photoionization Rate') .EQ. 0)
	     READ(20,'(A)')STRING
	     IF(INDEX(STRING,'Photoionization Rate') .EQ. 0)THEN
	       WRITE(21,'(A)')TRIM(STRING)
	     END IF
	     IF(INDEX(STRING,'Radius') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,'(A)')TRIM(STRING)
	       READ(STRING,*)(R(I),I=IST,IEND)
	     ELSE IF(INDEX(STRING,'Temperature') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,'(A)')TRIM(STRING)
	       READ(STRING,*)(TEMP(I),I=IST,IEND)
	     ELSE IF(INDEX(STRING,'Electron Density') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,'(A)')TRIM(STRING)
	       READ(STRING,*)(ED(I),I=IST,IEND)
	     ELSE IF(INDEX(STRING,'Ion Density') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,'(A)')TRIM(STRING)
	       READ(STRING,*)(DI(I),I=IST,IEND)
	     END IF
	   END DO
!
	   NLEV=0
	   DO WHILE(1 .EQ. 1)
	     READ(20,'(A)')STRING
	     IF(STRING .NE. ' ')THEN
	       NLEV=NLEV+1
	       READ(STRING,*)(PHOT(I,NLEV),I=IST,IEND)
	       PHOT_SUM(IST:IEND)=PHOT_SUM(IST:IEND)+PHOT(IST:IEND,NLEV)
	     ELSE
	       EXIT
	     END IF
	   END DO
!
	   WRITE(21,'(A,I4,A,I4,A)')'   Total photoionization rate (d=',IST,' to',IEND,'):'
	   WRITE(21,'(X,10ES12.4)')(PHOT_SUM(I),I=IST,IEND)
	   WRITE(21,'(A)')
!
	   DO WHILE(INDEX(STRING,'Recombination Rates') .EQ. 0)
	     READ(20,'(A)')STRING
	     IF(INDEX(STRING,'Recombination Rates') .EQ. 0)THEN
	       WRITE(21,'(A)')TRIM(STRING)
	     END IF
	     IF(INDEX(STRING,'Colisional Ionization Rate') .NE. 0)THEN
	        READ(20,'(A)')STRING
	        WRITE(21,'(A)')TRIM(STRING)
	        READ(STRING,*)(COL_IR(I),I=IST,IEND)
	     ELSE IF(INDEX(STRING,'Charge Transfer Ionization Rate') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,'(A)')TRIM(STRING)
	       READ(STRING,*)(CHG_IR(I),I=IST,IEND)
	     END IF
	   END DO
!
	   IF(NET_RECOM_PER_LEVEL)THEN
	     WRITE(21,'(A,A)')'   Net ',TRIM(STRING(4:))
	     DO J=1,NLEV
	       READ(20,*)(RECOM(I),I=IST,IEND)
	       RECOM_SUM(IST:IEND)=RECOM_SUM(IST:IEND)+RECOM(IST:IEND)
	       WRITE(21,'(X,10ES12.4)')(RECOM(I)-PHOT(I,J),I=IST,IEND)
	     END DO
	   ELSE
	     DO J=1,NLEV
	       READ(20,*)(RECOM(I),I=IST,IEND)
	       RECOM_SUM(IST:IEND)=RECOM_SUM(IST:IEND)+RECOM(IST:IEND)
	     END DO
	     WRITE(21,'(A,I4,A,I4,A)')'   Total recombination rate (d=',IST,' to',IEND,'):'
	     WRITE(21,'(X,10ES12.4)')(RECOM_SUM(I),I=IST,IEND)
!	     WRITE(21,'(A)')
	   END IF
!
	   DO WHILE(INDEX(STRING,'Net Recombination Rate') .EQ. 0)
	     READ(20,'(A)')STRING
	     IF(INDEX(STRING,'Net Recombination Rate') .EQ. 0)THEN
	       WRITE(21,'(A)')TRIM(STRING)
	     END IF
	     IF(INDEX(STRING,'Colisional Recombination Rate') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,'(A)')TRIM(STRING)
	       READ(STRING,*)(COL_RR(I),I=IST,IEND)
	     ELSE IF(INDEX(STRING,'Charge Transfer Recombination Rate') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,'(A)')TRIM(STRING)
	       READ(STRING,*)(CHG_RR(I),I=IST,IEND)
	     ELSE IF(INDEX(STRING,'Effective Advection Recombination Rate') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,'(A)')TRIM(STRING)
	       READ(STRING,*)(ADVEC_RR(I),I=IST,IEND)
	     ELSE IF(INDEX(STRING,'Non-Thermal Ionization') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,'(A)')TRIM(STRING)
	       READ(STRING,*)(NT_IR(I),I=IST,IEND)
	     END IF
	   END DO
!
	   WRITE(21,'(A)')TRIM(STRING)
	   FLUSH(21)
	   DO WHILE(1 .EQ. 1)
	     READ(20,'(A)',END=200)STRING
	     IF(STRING(1:1) .EQ. '1' .OR. STRING(1:1) .EQ. FORMFEED)THEN
	       EXIT
	       WRITE(21,'(A)')FORMFEED
	     END IF
	     WRITE(21,'(A)')TRIM(STRING)
	   END DO
	   FLUSH(21)
	END DO
200	CONTINUE
!
	WRITE(6,'(8ES14.4)')PHOT_SUM(1),RECOM_SUM(1),
	1             COL_IR(1),CHG_IR(1),COL_RR(1),CHG_RR(1),
	1             ADVEC_RR(1),NT_IR(1)
	WRITE(6,'(8ES14.4)')PHOT_SUM(ND),RECOM_SUM(ND),
	1             COL_IR(ND),CHG_IR(ND),COL_RR(ND),CHG_RR(ND),
	1             ADVEC_RR(ND),NT_IR(ND)
	WRITE(6,'(8ES14.4)')R(1),V(1),TEMP(1),R(ND),V(ND),TEMP(ND)
!
	DO I=1,ND
	   YSUM(I)=( PHOT_SUM(I)+RECOM_SUM(I) +
	1             COL_IR(I)+CHG_IR(I) +
	1             COL_RR(I)+CHG_RR(I) +
	1             ABS(NT_IR(I)) +
	1             ABS(ADVEC_RR(I)) )/2.0D0
	   PHOT_SUM(I)=PHOT_SUM(I)/YSUM(I)
	   COL_IR(I)=COL_IR(I)/YSUM(I)
	   CHG_IR(I)=CHG_IR(I)/YSUM(I)
	   NT_IR(I)=NT_IR(I)/YSUM(I)
!
	   RECOM_SUM(I)=-RECOM_SUM(I)/YSUM(I)
	   COL_RR(I)=-COL_RR(I)/YSUM(I)
	   CHG_RR(I)=-CHG_RR(I)/YSUM(I)
	   ADVEC_RR(I)=-ADVEC_RR(I)/YSUM(I)
	END DO
	YLABEL='Normalized rate'
!
2000	CONTINUE
	XAX_OPT='XN'
	CALL GEN_IN(XAX_OPT,'X axis option: XN, R, V, T, ED, EX(stop), NS (new species)')
	XAX_OPT=UC(XAX_OPT)
	IF(XAX_OPT(1:2) .EQ. 'XN')THEN
	  DO I=1,ND
	    XVEC(I)=I
	  END DO
	  XLABEL='Depth index'
	ELSE IF(XAX_OPT(1:2) .EQ. 'NS')THEN
	  GOTO 1000
	ELSE IF(XAX_OPT(1:1) .EQ. 'R')THEN
	  XVEC(1:ND)=R(1:ND)
	  XLABEL='Radius(10\u10\d cm)'
	ELSE IF(XAX_OPT(1:1) .EQ. 'V')THEN
	  XVEC(1:ND)=V(1:ND)
	  XLABEL='V(km/s)'
	ELSE IF(XAX_OPT(1:1) .EQ. 'T')THEN
	  XVEC(1:ND)=TEMP(1:ND)
	  XLABEL='T(10\u4\dK)'
	ELSE IF(XAX_OPT(1:2) .EQ. 'ED')THEN
	  XVEC(1:ND)=ED(1:ND)
	  XLABEL='Ne(cm\u-3\d)'
	ELSE IF(XAX_OPT(1:2) .EQ. 'EX')THEN
	  STOP
	ELSE
	  WRITE(6,*)'Unrecognized X-axis option'
	  GOTO 2000
	END IF
!
	MIN_VAL=0.01D0
	WRITE(6,*)' '
	WRITE(6,*)'Calling DP_CURVE'
	WRITE(6,*)'                : +ve means ionizing'
	WRITE(6,*)'                : -ve means recombination'
	WRITE(6,*)' '
!
	WRITE(6,'(A,I2,2A)')PG_PEN(2)//'   Curve ',1,': Photoionization rate',DEF_PEN
	CALL DP_CURVE(ND,XVEC,PHOT_SUM)
	I=2
	IF(MINVAL(RECOM_SUM(1:ND)) .LE. -MIN_VAL)THEN
	  CALL DP_CURVE(ND,XVEC,RECOM_SUM)
	  I=I+1
	  WRITE(6,'(A,I2,2A)')PG_PEN(I)//'   Curve ',I-1,': Radiative recombination rate',DEF_PEN
	END IF
!
	IF(MAXVAL(NT_IR(1:ND)) .GE. MIN_VAL)THEN
	  CALL DP_CURVE(ND,XVEC,NT_IR)
	  I=I+1
	  WRITE(6,'(A,I2,2A)')PG_PEN(I)//'   Curve ',I-1,': Non-thermal ionization rate',DEF_PEN
	END IF
	IF(MAXVAL(ABS(ADVEC_RR(1:ND))) .GE. MIN_VAL)THEN
	  CALL DP_CURVE(ND,XVEC,ADVEC_RR)
	  I=I+1
	  WRITE(6,'(A,I2,2A)')PG_PEN(I)//'   Curve ',I-1,': Advection rate',DEF_PEN
	END IF
!
	IF(MAXVAL(CHG_IR(1:ND)) .GE. MIN_VAL)THEN
	  CALL DP_CURVE(ND,XVEC,CHG_IR)
	  I=I+1
	  WRITE(6,'(A,I2,2A)')PG_PEN(I)//'   Curve ',I-1,': Charge exch. ionization rate',DEF_PEN
	END IF
	IF(MINVAL(CHG_RR(1:ND)) .LE. -MIN_VAL)THEN
	  CALL DP_CURVE(ND,XVEC,CHG_RR)
	  I=I+1
	  WRITE(6,'(A,I2,2A)')PG_PEN(I)//'   Curve ',I-1,': Charge exch. recombination rate',DEF_PEN
	END IF
!
	IF(MAXVAL(COL_IR(1:ND)) .GE. MIN_VAL)THEN
	  CALL DP_CURVE(ND,XVEC,COL_IR)
	  I=I+1
	  WRITE(6,'(A,I2,2A)')PG_PEN(I)//'   Curve ',I-1,': Collsional ionization rate',DEF_PEN
	END IF
	IF(MINVAL(COL_RR(1:ND)) .LE. -MIN_VAL)THEN
	  CALL DP_CURVE(ND,XVEC,COL_RR)
	  I=I+1
	  WRITE(6,'(A,I2,2A)')PG_PEN(I)//'   Curve ',I-1,': Collisional recombination rate',DEF_PEN
	END IF
!
	WRITE(6,*)' '
	CALL GRAMON_PGPLOT(XLABEL,YLABEL,' ',' ')
	GOTO 2000
!
	END
