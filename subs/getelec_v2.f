!
	SUBROUTINE GETELEC_V2(SPEC_DENS,AT_NO,N_SPEC,ED,ND,LU,FILNAME)
	IMPLICIT NONE
!
! Altered 24-Apr-2004 : Check monoticty of ED. Bug fixed on 24th.
! Altered 04-Mar-2003 : Use grid index, rather than average gamma, to provide
!                         first estimated for the electron density. This is
!                         taken as the final estimate if the electron density
!                         is not monotonic.
! Altered 11-Jun-1996 : Populations of all specied no passed through matrix
!                         SPEC_DENS. AT_NO must also be passed.
!                        
! Altered 24-May-1996 : Dynamic memmoray allocation now used.
!                       ERROR_LU and LUER installed.
! Altered 18-Apr-1990 - Change correction had bug. Comparing ED(log)
!                       with non log ED.
!
	INTEGER N_SPEC,ND,LU
	REAL*8 ED(ND)
	REAL*8 SPEC_DENS(ND,N_SPEC)
	REAL*8 AT_NO(N_SPEC)
	CHARACTER*(*) FILNAME
C
	REAL*8, ALLOCATABLE :: OLD_ED(:)		!NGAM
	REAL*8, ALLOCATABLE :: OLD_GAM(:,:)		!NGAM,N_SPEC
	REAL*8, ALLOCATABLE :: NEW_GAM(:)		!ND
	REAL*8, ALLOCATABLE :: NEW_ED(:)		!ND
C
	CHARACTER*132 STRING
	INTEGER I,J,K,NGAM
	INTEGER ID,XST,XEND,NX
	REAL*8 CHANGE,T1,T2
	REAL*8 LOC_AT_NO
C
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
C
	LUER=ERROR_LU()
C
	OPEN(UNIT=LU,FILE=FILNAME,STATUS='OLD')
	  STRING=' '
	  DO WHILE( INDEX(STRING,'!Number of depth points') .EQ. 0)
	    READ(LU,'(A)')STRING
	  END DO
	  READ(STRING,*)NGAM
C
C Now allocate required stroage.
C
	  ALLOCATE (OLD_ED(NGAM))
	  ALLOCATE (OLD_GAM(NGAM,N_SPEC))
	  ALLOCATE (NEW_GAM(ND))
	  ALLOCATE (NEW_ED(ND))
C
	  OLD_GAM(:,:)=0.0D0
C
	  STRING=' '
	  DO WHILE( INDEX(STRING,'!Electron density') .EQ. 0)
	    READ(LU,'(A)')STRING
	  END DO
	  READ(LU,*)(OLD_ED(J),J=1,NGAM)
C
C Read in the gammas for species present in the input file. We link
C them with the populations stored in SPEC_DENS via their atomic number.
C The red continues until an end of file condition occurs.
C
	  DO WHILE(1 .EQ. 1)
	    STRING=' '
	    DO WHILE( INDEX(STRING,'!Atomic N') .EQ. 0)
	      READ(LU,'(A)',END=100)STRING
	    END DO
	    READ(STRING,*)LOC_AT_NO
	    ID=0
	    DO J=1,N_SPEC
	      IF(LOC_AT_NO .EQ. AT_NO(J))ID=J
	    END DO
	    IF(ID .NE. 0)THEN
	      READ(LU,*)(OLD_GAM(J,ID), J=1,NGAM)
	    END IF
	  END DO
100	  CONTINUE		!End of file
	CLOSE(UNIT=LU)
!
! Estimate the GAMMAS by simply using the grid index.
!
	ED(:)=0.0D0	  
	DO ID=1,N_SPEC
	  DO I=1,ND
	    J=I+(NGAM-ND)
	    IF(J .LT. 1)J=1
	    IF(J .GT. NGAM)J=NGAM
	    ED(I)=ED(I)+OLD_GAM(J,ID)*SPEC_DENS(I,ID)
	  END DO
	END DO
!
! Ie ED is non-monotnoic, we have no method implemented to improve the
! ED estimate.
!
	DO I=1,NGAM-2
	  IF( (OLD_ED(I)-OLD_ED(I+1))*(OLD_ED(I+1)-OLD_ED(I+2)) .LE. 0 )THEN
	    WRITE(LUER,*)'Warning ED in GETELEC_V2 is non-monotonic'
	    WRITE(LUER,*)'We estimated the gammas using the grid index.'
	    RETURN
	  END IF
	END DO
!
	OLD_ED(:)=LOG(OLD_ED(:))
	ED(:)=LOG(ED(:))
C
C Improve estimate of the electron density. We do this loop a maximum of
C 5 times. Convergence should be fast for a H or He domiated atmosphere.
C
	DO K=1,5
C
	  XST=1
	  XEND=ND
	  DO I=1,ND
	    J=ND-I+1
	    IF(ED(I) .LT. OLD_ED(1))XST=I+1
	    IF(ED(J) .GT. OLD_ED(NGAM))XEND=J-1
	  END  DO
	  NX=XEND-XST+1
C
	  IF(K .EQ. 1)THEN
	    WRITE(LUER,*)'Estimating electron density as a function of depth.'
	    WRITE(LUER,*)'Illustration of ED convergence:'
	  END IF
	  WRITE(LUER,*)EXP(ED(XST)),EXP(ED(XEND))
C
	  NEW_ED(1:ND)=0.0D0
	  DO J=1,N_SPEC
	    CALL LIN_INTERP(ED(XST),NEW_GAM(XST),NX,
	1             OLD_ED,OLD_GAM(1,J),NGAM)
	    DO I=XST,XEND
	      NEW_ED(I)=NEW_ED(I)+NEW_GAM(I)*SPEC_DENS(I,J)
	    END DO
	  END DO
C
C For the extrapolation, we assume constant ionization.
C
	  DO I=1,XST-1
	    T1=0.0D0
	    T2=0.0D0
	    DO J=1,N_SPEC
	      T1=T1+SPEC_DENS(I,J)
	      T2=T2+SPEC_DENS(XST,J)
	    END DO
	    NEW_ED(I)=NEW_ED(XST)*T1/T2
	  END DO
C
	  DO I=XEND+1,ND
	    T1=0.0D0
	    T2=0.0D0
	    DO J=1,N_SPEC
	      T1=T1+SPEC_DENS(I,J)
	      T2=T2+SPEC_DENS(XEND,J)
	    END DO
	    NEW_ED(I)=NEW_ED(XEND)*T1/T2
	  END DO
C
	  CHANGE=0.0D0
	  DO I=1,ND
	    T1=LOG(NEW_ED(I))
	    CHANGE=MAX( CHANGE,ABS(ED(I)-T1) )
	    ED(I)=T1
	  END DO
	  IF(CHANGE .LT. 0.01)GOTO 500
	END DO
	WRITE(LUER,*)'Warning ED calculation has not converged',
	1               ' % change is ',CHANGE*LOG(10.0)
C
500	CONTINUE
C
        DO I=1,ND
	  ED(I)=NEW_ED(I)
	END DO
C
C Deallocate requested storage.
C
	DEALLOCATE (OLD_ED)
	DEALLOCATE (OLD_GAM)
	DEALLOCATE (NEW_GAM)
	DEALLOCATE (NEW_ED)
C
	RETURN
	END
