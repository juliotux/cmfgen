!
! Small subroutine to read in BA_ASCI_N_D? file. These files are output
! by CMFGEN, and allow the matrix soultions and corrections to be checked.
!
! This routine is called by CHK_BA_EST.
!
	SUBROUTINE RD_BA_MAT(CMAT,STEQ,POPS,N,FILENAME,LU)
!
! Altered: 6-Dec-2009
!
	INTEGER N
	INTEGER LU
	REAL*8 POPS(N)
	REAL*8 STEQ(N)
	REAL*8 CMAT(N,N)
	CHARACTER(LEN=*)FILENAME
!
	INTEGER I
	CHARACTER(LEN=132)STRING
!
	OPEN(UNIT=LU, FILE=FILENAME,STATUS='OLD',ACTION='READ')
	  STRING=' '
	  DO WHILE(INDEX(STRING,'POP') .EQ. 0)		
	    READ(LU,'(A)')STRING
	  END DO
	  READ(LU,'(A)')STRING
	  DO WHILE(STRING .EQ. ' ')
	    READ(LU,'(A)')STRING
	  END DO
	  DO I=1,N
	    READ(STRING(10:),*)POPS(I)
	    READ(LU,'(A)')STRING
	  END DO
	  WRITE(6,*)'Successfully read POPS'
!
	  STRING=' '
	  DO WHILE(INDEX(STRING,'STEQ') .EQ. 0)		
	    READ(LU,'(A)')STRING
	  END DO
	  READ(LU,'(A)')STRING
	  DO WHILE(STRING .EQ. ' ')
	    READ(LU,'(A)')STRING
	  END DO
	  DO I=1,N
	    READ(STRING(10:),*)STEQ(I)
	    READ(LU,'(A)')STRING
	  END DO
	  WRITE(6,*)'Successfully read STEQ'
!
	  STRING=' '
	  DO WHILE(INDEX(STRING,'C_MAT') .EQ. 0)		
	    READ(LU,'(A)')STRING
	  END DO
!
	  DO K=1,N,5
	    STRING=' '
	    DO WHILE(STRING .EQ. ' ')
	      READ(LU,'(A)')STRING
	    END DO
	    DO I=1,N
	      IF(I .NE. 1)READ(LU,'(A)')STRING
	      READ(STRING(10:),*)(CMAT(I,J),J=K,MIN(K+4,N))
	    END DO
	  END DO
	  WRITE(6,*)'Successfully read CMAT'
!
	CLOSE(UNIT=LU)
!
	RETURN
	END
