C
C Routine to solve for the perturbations to the population
C parameters using gaussian elimination. Due to numerical
C instabilities it is necessary to rewrite the equations such
C that the 'depth parameters' ocurr first in the system of
C equations.
C
	SUBROUTINE NAGMINSOL(BA,STEQ,ABT,FQ,NV,ND)
	IMPLICIT NONE
C
C Altered 24-May-1996 - IMPLICIT NONE installed.
C                       IONE passed in calls
C                       ERROR_LU installed.
C Altered 27-Nov-1987 - Altered so that matrix can be restructured in situ.
C                        Dimension of ABT was changed to NV*NV*ND (from
C                        (NV*NV*ND*ND), and ABT was replaced BA in solution
C                        calls.
C Altered 16-APRIL-85 - Nag rouitines called to solve the similtaneous
C                        equations. This routine assumes a paged enviroment.
C Changed 22-JUL-82 -total number of variables passed
C
	INTEGER ND,NV
	REAL*8 BA(NV,NV,ND,ND),STEQ(NV,ND)
	REAL*8 ABT(NV*NV*ND),FQ(NV*ND)
C
C Local variables.
C
	REAL*8 DP
C
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
C
	INTEGER I,J,NU,IFAIL
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
	NU=ND*NV
	LUER=ERROR_LU()
C
C Restructure the order of the eqations.
C
	CALL TUNE(IONE,'RESTR')
	CALL TRANSPOSEBA(BA,BA,ABT,ABT,NV,ND)
C	DO K=1,ND
C	  DO J=1,NV
C	    DO L=1,ND
C	      DO I=1,NV
C	        ABT(NU+1-(I+(L-1)*NV),NU+1-(J+(K-1)*NV))=BA(I,J,K,L)
C	      END DO
C	    END DO
C	  END DO
C	END DO
	CALL TUNE(ITWO,'RESTR')
C
	DO J=1,ND
	  DO I=1,NV
	    FQ(NU+1-(J-1)*NV-I)=STEQ(I,J)
	  END DO
	END DO
C
C
C Solve for the perturbations. STEQ is used as the pivot array. FQ contains
C the single RHS on entry to F04AYF and on exit contains the solution vector.
C
C	CALL GAUSEL(BA,FQ,STEQ,NU,IONE,KS)	!Work area BA changed
	IFAIL=0
	CALL F01BTF(NU,BA,NU,STEQ,DP,IFAIL)
	IF(IFAIL .NE. 0)THEN
	  WRITE(LUER,*)' ERROR in F01BTF (Nagminsol)'
	  WRITE(LUER,*)' IFAIL=',IFAIL
	  STOP
	END IF
	CALL F04AYF(NU,IONE,BA,NU,STEQ,FQ,NU,IFAIL)
	IF(IFAIL .NE. 0)THEN
	  WRITE(LUER,*)' ERROR in F04AYF (Nagminsol)'
	  WRITE(LUER,*)' IFAIL=',IFAIL
	  STOP
	END IF
C
C Undo restructuring.
C
	DO I=1,ND
	  DO J=1,NV
	    STEQ(J,I)=FQ(NV*(ND-I+1)-J+1)
	  END DO
	END DO
C
	RETURN
	END
