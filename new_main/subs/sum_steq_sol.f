!
! This routine ouputs information regarding the levels which
! show the largest corrections. Ouput is to CORRECTION_LINK.
! Changes at 5 (up to 10 if max and min depths differ) are output.
! Information output includes % change, species, and level.
!
	SUBROUTINE SUM_STEQ_SOL(SOL,NT,ND,LUOUT)
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Altered:   05-Aug-2011 : Also output level ID in STEQ/SOL array.
! Altered:   05-Apr-2011 : cleaning
! Finalized: 01-Feb-2011
!
	INTEGER LUOUT
	INTEGER NT
	INTEGER ND
	REAL*8 SOL(NT,ND)
!
	INTEGER INDX(NT)
	INTEGER VEC_SL(NT)
	CHARACTER(LEN=12) VEC_DESC(NT)
!
	REAL*8 MIN_CHANGE,MAX_CHANGE
	INTEGER ID,I,J,K,L
	INTEGER LMAX,LMIN
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
!
	VEC_SL=0
	DO ID=1,NUM_IONS-1
	  IF(ATM(ID)%XzV_PRES)THEN
	    DO J=ATM(ID)%EQXzV,ATM(ID)%EQXzV+ATM(ID)%NXzV-1
	      VEC_DESC(J)=ION_ID(ID)
	      VEC_SL(J)=J+1-ATM(ID)%EQXzV
	    END DO
	    J=ATM(ID)%EQXzV+ATM(ID)%NXzV
	    VEC_DESC(J)=ION_ID(ID+1)
	    VEC_SL(J)=1
	  END IF
	END DO
	VEC_DESC(NT-1)='Ne'
	VEC_DESC(NT)='T'
!
	OPEN(FILE='CORRECTION_LINK',STATUS='UNKNOWN',ACTION='WRITE',UNIT=LUOUT)
!
	DO K=1,5
	  MAX_CHANGE=0.0D0
	  MIN_CHANGE=0.0D0
	  LMAX=K
	  LMIN=5+K
!
	  DO L=1,ND
	    DO I=1,NT
	      IF(SOL(I,L) .LT. MAX_CHANGE)THEN
	        MAX_CHANGE=SOL(I,L)
	        LMAX=L
	      ELSE IF(SOL(I,L) .GT. MIN_CHANGE)THEN
	        MIN_CHANGE=SOL(I,L)
	        LMIN=L
	      END IF
	    END DO
	  END DO
!
! We now print out the 5 largest corrections at each of these depths.
!
	  L=LMAX
	  WRITE(LUOUT,'(A)')' '
	  DO WHILE(1 .EQ. 1)
	    WRITE(LUOUT,'(A)')' '
	    WRITE(LUOUT,'(A,I4)')' 5 largest reductions at depth:',L
	    WRITE(LUOUT,'(A)')' '
	    WRITE(LUOUT,'(8X,A,3X,A,6X,A,5X,A)')'SOL(J,L)','Species','SL','I(STEQ)'
	    CALL INDEXX(NT,SOL(1,L),INDX,L_TRUE)
	    DO I=NT,NT-4,-1
	      J=INDX(I)
	      WRITE(LUOUT,'(2X,ES14.4,A10,I8,6X,I5)')SOL(J,L),TRIM(VEC_DESC(J)),VEC_SL(J),J
	    END DO
!
	    WRITE(LUOUT,'(A)')' '
	    WRITE(LUOUT,'(A,I4)')' 5 largest increases at depth:',L
	    WRITE(LUOUT,'(A)')' '
	    DO I=1,5
	      J=INDX(I)
	      WRITE(LUOUT,'(2X,ES14.4,A10,I8,6X,I5)')SOL(J,L),TRIM(VEC_DESC(J)),VEC_SL(J),J
	    END DO
	    IF(L .EQ. LMIN)EXIT
	    L=LMIN 
	  END DO
	  SOL(:,LMIN)=0.0D0
	  SOL(:,LMAX)=0.0D0
	  WRITE(LUOUT,'(A)')' '
	  WRITE(LUOUT,'(A)')' '
!
	END DO
	CLOSE(UNIT=LUOUT)
!
	RETURN
	END
