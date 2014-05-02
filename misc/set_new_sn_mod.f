        PROGRAM NEW_SN_MOD
        USE GEN_IN_INTERFACE
        IMPLICIT NONE
! 
! Altered: 29-Jan-2014 : Extensive improvements 
!
! Simple program to modify the VADAT for a new MODEL at an
! advanced time step.
!
! Two cases:
!           (1) Old model has TS_NO > 1
!               Routine changes AGE, TS_NO, LIN_INT, and FIX_T only.
!
!           (2) Old models has TS_NO=1
!                 Routine changes all switches necesary to go from static
!                 to a fully time-dependent model.
!
! Some othe common parameters which often vary for a promted for, or, in
! the case of parameters associated with NG aceleration, set to default values.
!
	REAL*8 FRAC_TIME_STEP
	REAL*8 AGE
	REAL*8 T1
!
        INTEGER, PARAMETER :: NMAX=2000
        INTEGER I,J
        INTEGER L
	INTEGER KB	
	INTEGER KT	
        INTEGER NREC
        INTEGER IOS
        INTEGER TS_NO
!
	LOGICAL ANSWER
!
	CHARACTER(LEN=1) TAB
        CHARACTER(LEN=132) LINE(NMAX)
        CHARACTER(LEN=132) STRING
        CHARACTER(LEN=20) TMP_STR
!
	TAB=CHAR(9)
        OPEN(UNIT=10,FILE='VADAT',STATUS='OLD')
        I=0
        DO WHILE(I .LT. NMAX-1)
          READ(10,'(A)',IOSTAT=IOS)LINE(I+1)
          IF(IOS .NE. 0)EXIT
          I=I+1
        END DO
        NREC=I
        REWIND(UNIT=10)
        WRITE(6,'(A,I4,A)')'Read',NREC,'records from VADAT'
        WRITE(6,*)'Successfully read VADAT file'
!
        DO L=1,NREC
!
! We first remove TABS, and line up the [.
!
	  KT=1
	  DO WHILE(KT .NE. 0)
	    KT=INDEX(LINE(L),TAB)
	    IF(KT .NE. 0)LINE(L)(KT:KT)=' '
	  END DO
!
	  KB=INDEX(LINE(L),'[')
	  IF(KB .LT. 14 .AND. KB .NE. 0)THEN
	    STRING=LINE(L)(KB:)
	    LINE(L)(KB-1:)=' '
	    LINE(L)(14:)=STRING
	    KB=14
	  END IF
!
	  J=INDEX(LINE(L),'!')
	  IF(J .NE. 0)THEN
	    STRING=LINE(L)(J:)
	    LINE(L)(J:)=' '
	    LINE(L)(35:)=STRING
	  END IF
!
          IF(INDEX(LINE(L),'[SN_AGE]') .NE. 0)THEN
            FRAC_TIME_STEP=1.1D0
            CALL GEN_IN(FRAC_TIME_STEP,'Fractional time increment -- e.g. 1.10 for a 10% increase')
            READ(LINE(L),*)AGE
            AGE=AGE*FRAC_TIME_STEP
            WRITE(TMP_STR,'(F8.4,A)')AGE,'D0'
            TMP_STR=ADJUSTL(TMP_STR)
	    LINE(L)(1:KB-1)=' '
	    J=LEN_TRIM(TMP_STR)
            LINE(L)(1:J)=TRIM(TMP_STR)
!
          ELSE IF(INDEX(LINE(L),'[TS_NO]') .NE. 0)THEN
            READ(LINE(L),*)TS_NO
            TS_NO=TS_NO+1
            WRITE(TMP_STR,'(I3)')TS_NO
            TMP_STR=ADJUSTL(TMP_STR)
	    LINE(L)(1:KB-1)=' '
	    LINE(L)(1:3)=TMP_STR(1:3)
!
          ELSE IF(INDEX(LINE(L),'[FIX_T]') .NE. 0) THEN
	    LINE(L)(1:KB-1)=' '
            LINE(L)(1:1) = 'T'
!
          ELSE IF(INDEX(LINE(L),'[LIN_INT]') .NE. 0)THEN
	    LINE(L)(1:KB-1)=' '
            LINE(L)(1:1)='F'
	  END IF
	END DO
!
! Changes made by Kevin (DEC. 17, 2013)
! The following changes only need to be made when go from TIME STEP 1 to TIME STEP 2.
!
	IF(TS_NO .EQ. 2)THEN
	  DO L=1,NREC
	    LINE(L)=ADJUSTL(LINE(L))
	    KB=INDEX(LINE(L),'[')
            IF(INDEX(LINE(L),'[DO_DDT]') .NE. 0) THEN
	      LINE(L)(1:KB-1)=' '
              LINE(L)(1:1) = 'T'
            ELSE IF(INDEX(LINE(L),'[USE_J_REL]') .NE. 0) THEN
	      LINE(L)(1:KB-1)=' '
              LINE(L)(1:1) = 'F'
            ELSE IF(INDEX(LINE(L),'[INCL_REL]') .NE. 0) THEN
	      LINE(L)(1:KB-1)=' '
              LINE(L)(1:1) = 'F'
            ELSE IF(INDEX(LINE(L),'[INCL_ADV_TRANS]') .NE. 0) THEN
	      LINE(L)(1:KB-1)=' '
              LINE(L)(1:1) = 'F'
            ELSE IF(INDEX(LINE(L),'[INCL_DJDT]') .NE. 0) THEN
	      LINE(L)(1:KB-1)=' '
              LINE(L)(1:1) = 'T'
            ELSE IF(INDEX(LINE(L),'[USE_DJDT_RTE]') .NE. 0) THEN
	      LINE(L)(1:KB-1)=' '
              LINE(L)(1:1) = 'T'
            ELSE IF(INDEX(LINE(L),'[DC_METH]') .NE. 0) THEN
	      LINE(L)(1:KB-1)=' '
              LINE(L)(1:1) = 'R'
            ELSE IF(INDEX(LINE(L),'[INC_AD]') .NE. 0) THEN
	      LINE(L)(1:KB-1)=' '
              LINE(L)(1:1) = 'T'
            ELSE IF(INDEX(LINE(L),'[LTE_EST]') .NE. 0) THEN
	      LINE(L)(1:KB-1)=' '
              LINE(L)(1:1) = 'F'
            END IF
	  END DO
	END IF
!
! The following options generally do need not change. Current values are
! used as defaults.
!
	ANSWER=.FALSE.
	CALL GEN_IN(ANSWER,'Change other parameters:')
	IF(ANSWER)THEN
	  DO L=1,NREC
	    KB=INDEX(LINE(L),'[')
            IF(INDEX(LINE(L),'[GF_CUT]') .NE. 0) THEN
              READ(LINE(L),*)T1
	      CALL GEN_IN(T1,'GF_CUT -- 1.0D-03')
	      LINE(L)(1:KB-1)=' '
              WRITE(TMP_STR,'(ES8.2)')T1; TMP_STR=ADJUSTL(TMP_STR)
              LINE(L)(1:8)=TMP_STR(1:8)
            ELSE IF(INDEX(LINE(L),'[T_MIN]') .NE. 0) THEN
              READ(LINE(L),*)T1
	      CALL GEN_IN(T1,'T_MIN -- 0.30')
	      LINE(L)(1:KB-1)=' '
              WRITE(TMP_STR,'(F5.3,A2)')T1,'D0'; TMP_STR=ADJUSTL(TMP_STR)
              LINE(L)(1:7)=TMP_STR(1:7)
            ELSE IF(INDEX(LINE(L),'[T_INIT_TAU]') .NE. 0) THEN
              READ(LINE(L),*)T1
	      CALL GEN_IN(T1,'T_INIT_TAU -- 10.0')
	      LINE(L)(1:KB-1)=' '
              WRITE(TMP_STR,'(F5.2,A2)')T1,'D0'; TMP_STR=ADJUSTL(TMP_STR)
              LINE(L)(1:7)=TMP_STR(1:7)
            ELSE IF(INDEX(LINE(L),'[GREY_TAU]') .NE. 0) THEN
              READ(LINE(L),*)T1
	      LINE(L)(1:KB-1)=' '
	      CALL GEN_IN(T1,'GREY_TAU -- 5.0')
              WRITE(TMP_STR,'(F5.2,A2)')T1,'D0'; TMP_STR=ADJUSTL(TMP_STR)
              LINE(L)(1:7)=TMP_STR(1:7)
	    ELSE IF(INDEX(LINE(L),'[AMASS_DOP]') .NE. 0) THEN
              WRITE(6,'(A)')LINE(L)
	      READ(LINE(L),*)T1
	      CALL GEN_IN(T1,'AMASS_DOP -- 1.0D+10')
	      LINE(L)(1:KB-1)=' '
              WRITE(TMP_STR,'(ES7.1)')T1; TMP_STR=ADJUSTL(TMP_STR)
              LINE(L)(1:7)=TMP_STR(1:7)
            ELSE IF(INDEX(LINE(L),'[N_PAR]') .NE. 0) THEN
              READ(LINE(L),*)I
	      CALL GEN_IN(I,'N_PAR: 5000')
	      LINE(L)(1:KB-1)=' '
              WRITE(TMP_STR,'(I5)')I; TMP_STR=ADJUSTL(TMP_STR)
              LINE(L)(1:5)=TMP_STR(1:5)
            ELSE IF(INDEX(LINE(L),'[FIX_T_AUTO]') .NE. 0) THEN
	      LINE(L)(1:KB-1)=' '
              LINE(L)(1:1) = 'F'
            END IF
	  END DO
	END IF
!
	ANSWER=.FALSE.
	CALL GEN_IN(ANSWER,'Set NG parameter to default values?')
	IF(ANSWER)THEN
	  DO L=1,NREC
	    KB=INDEX(LINE(L),'[')
	    IF(INDEX(LINE(L),'[DO_NG]') .NE. 0) THEN
	      LINE(L)(1:KB-1)=' '
              LINE(L)(1:1) = 'T'
            ELSE IF(INDEX(LINE(L),'[BEG_NG]') .NE. 0) THEN
	      LINE(L)(1:KB-1)=' '
              LINE(L)(1:5) = '5.0D0'
            ELSE IF(INDEX(LINE(L),'[IBEG_NG]') .NE. 0) THEN
	      LINE(L)(1:KB-1)=' '
              LINE(L)(1:2) = '20'
            ELSE IF(INDEX(LINE(L),'[BW_NG]') .NE. 0) THEN
	      LINE(L)(1:KB-1)=' '
              LINE(L)(1:3) = '200'
            ELSE IF(INDEX(LINE(L),'[ITS/NG]') .NE. 0) THEN
	      LINE(L)(1:KB-1)=' '
              LINE(L)(1:2) = '20'
            ELSE IF(INDEX(LINE(L),'[DO_UNDO]') .NE. 0) THEN
	      LINE(L)(1:KB-1)=' '
              LINE(L)(1:1) = 'F'
	    END IF
	  END DO
	END IF
!
	DO L=1,NREC
          LINE(L)=ADJUSTL(LINE(L))
          WRITE(10,'(A)')TRIM(LINE(L))
        END DO
	CLOSE(UNIT=10)
        WRITE(6,'(/,A,I4,A,/)')' Written ',NREC,' records  to VADAT'
!
	ANSWER=.TRUE.
	CALL GEN_IN(ANSWER,'Generate default IN_ITS file')
	IF(ANSWER)THEN
          OPEN(UNIT=10,FILE='IN_ITS',STATUS='UNKNOWN')
	    WRITE(10,'(A)')'200            [NUM_ITS]               !Numer of iterations?'
	    WRITE(10,'(A)')'T              [DO_LAM_IT]             !Do lambda iterations?'
	    WRITE(10,'(A)')'T              [DO_LAM_AUTO]           !Switch automatically from lambda its.'
	    WRITE(10,'(A)')'T              [DO_T_AUTO]             !Automatically switch to variable T'
	    WRITE(10,'(A)')'F              [DO_GT_AUTO]            !'
	    WRITE(10,'(A)')'F              [D2_EQ_D1]              !Only for testing'
	  CLOSE(UNIT=10)
	END IF
        STOP
        END
