!
! Subroutine to determine those frequencies at which the continuum opacity 
! will be evaluated. At other frequenecies it can be held fixed.
!
	SUBROUTINE DET_MAIN_CONT_FREQ_V2(NU,NCF,NU_CONT,NCF_CONT,
	1                NU_EVAL,MID_POINT_EVAL,COMPUTE_ALL_CROSS)
	IMPLICIT NONE
!
! NCF represents the toal number of frequencies at which we solve the radiative
! transfer equation. The frequencies are stored in the vector NU which
! has been chosen to sample all lines and continuum edges adequately.
!
	INTEGER*4 NCF
	REAL*8 NU(NCF)
!
! For each frequncy NU(I), the continuum cross-section will be evaluated at
! freqency NU_EVAL(I).
!
	REAL*8 NU_EVAL(NCF)
	REAL*8 NU_EVAL_SAV(NCF)
!
! NCF_CONT represents the number of continuum frequencies before line 
! insertion. NU_CONT contains these frequencies, and has been chosen to 
! sample continuum cross-sections (with allowance for level dissolution and 
! important bound-free edges) adequately. 
!
	INTEGER*4 NCF_CONT
	REAL*8 NU_CONT(NCF_CONT)
!
	LOGICAL COMPUTE_ALL_CROSS
	LOGICAL MID_POINT_EVAL
!
	INTEGER*4 ERROR_LU
	EXTERNAL ERROR_LU
!
! Local variables:
!
	INTEGER*4 LU_OUT,L,K,ML,ML_ST,ML_END,LST_COMP
	REAL*8 T1,T2,DOP_RAT,VRAT,C_KMS
!
	LU_OUT=ERROR_LU()
	IF(COMPUTE_ALL_CROSS)THEN
	  NU_EVAL(:)=NU(:)
	  WRITE(LU_OUT,'(A)')'The continuum will be evaluated at all',
	1                        ' frequencies.'
	  RETURN
	END IF
!
	WRITE(6,*)'Entering Det'
!
! We set NU_EVAL to zero. Subsequently a NU_EVAL value of zero for any index
! indicates that NU_EVAL still needs to be set.
!
	NU_EVAL(:)=0.0D0
!
! We first ensure that the continuum opacity is evaluated at all important
! bound_free edges etc. These are located in NU_CONT. We use DOPV to do this
! since all continuum edges should be in NU unless they were in DOPV km/s
! of an inserted line frequency.
!
	K=1
	DO ML=2,NCF_CONT-1
	  DO WHILE( NU(K)-NU_CONT(ML) .GT. 0)
	    K=K+1
	  END DO
	  NU_EVAL(K-1)=NU(K-1)
	END DO
	NU_EVAL(1)=NU(1)
	NU_EVAL(NCF)=NU(NCF)
	NU_EVAL_SAV(1:NCF)=NU_EVAL(1:NCF)
!
! Continuum will be evaluated at the mid point of each frequency band.
! The highest frequecny in the band is not adjusted. 
!
	IF(MID_POINT_EVAL)THEN
          ML_ST=1
	  DO WHILE(ML_ST .LT. NCF)
	    IF(NU_EVAL(ML_ST+1) .EQ. 0)THEN
	      ML=ML_ST+1
	      DO WHILE (NU_EVAL(ML+1) .EQ. 0)
	        ML=ML+1
	      END DO
              T1=0.5D0*(NU(ML_ST)+NU(ML+1))
	      DO K=ML_ST+1,ML+1
	        NU_EVAL(K)=T1
	      END DO
	      ML_ST=ML+1
	    ELSE
	      ML_ST=ML_ST+1
	    END IF
	  END DO 
	ELSE
!
! Continuum is evaluated at highest frequency in band. This is used
! at all internal band points, until new continuum edge is reached.
!
	  DO ML=1,NCF
	    IF(NU_EVAL(ML) .EQ. 0)NU_EVAL(ML)=NU_EVAL(ML-1)
	  END DO 
	END IF
!
	K=1
	DO ML=2,NCF
	  IF(NU_EVAL(ML) .NE. NU_EVAL(ML-1))K=K+1
	END DO
!
	WRITE(LU_OUT,'(A,I8,A)')' The continuum will be evaluated at ',
	1                        K,' frequencies'
	WRITE(LU_OUT,'(A)')' '
!
	C_KMS=2.998D+05
	OPEN(UNIT=27,FILE='CONT_FREQ',STATUS='UNKNOWN')
	   WRITE(27,*)NCF,'        !NCF'
	   T1=NU(1)
	   DO ML=1,NCF
	     WRITE(27,'(X,I6,3X,3ES14.5)')ML,NU(ML),NU_EVAL_SAV(ML),NU_EVAL(ML)
!	     IF(NU_EVAL(ML) .EQ. NU(ML))THEN
!	       WRITE(27,'(X,I6,3X,2ES14.5,F10.2)')ML,NU(ML),0.01D0*C_KMS/NU(ML),C_KMS*(NU(ML)-T1)/NU(ML)
!	       T1=NU(ML)
!	     END IF
	   END DO
	CLOSE(UNIT=27)
!
	RETURN
	END
