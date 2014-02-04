C
C This routine is used to compute the perturbations to
C Jv (1 ... ND) as a function of the population levels and T .
C Uses Schuster or diffusion approximation for lower boundary
C condition. Subroutine may be used with or without a variable
C temperature. The FEAUTRIER technique is used.
C
	SUBROUTINE PERTJFEAU_IBC(F2DA,FC,FA,DTAU,CHI,R,ZETA,THETA,
	1         RJ,Q,F,dCHIdR,TA,TB,TC,HBC_J,HBC_S,
	1         INBC,DBB,DIFF,THK,ND,METHOD)
	IMPLICIT NONE
C
C Altered 25-May-1996 - Call to DP_ZERO removed.
C                        IONE insertedin call to SIMPTH
C
C Altered 12-JUN-1991 - HBC replaced by HBC_J and HBC_S. Name changed from
C                     PERTJFEAUNEW.
C Altered 26-FEB-1986 - Bug fixed)
C Created 18-FEB-1986
C
	INTEGER ND
	LOGICAL DIFF,THK
	CHARACTER*6 METHOD
	REAL*8 F2DA(ND,ND),FC(ND,ND),FA(ND),DTAU(ND),CHI(ND),R(ND)
	REAL*8 ZETA(ND),THETA(ND),RJ(ND),Q(ND),F(ND)
	REAL*8 dCHIdR(ND),TA(ND),TB(ND),TC(ND)
	REAL*8 DBB,HBC_J,HBC_S,INBC
C
C Functions called.
C
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
C
C Local variables.
C
	INTEGER, PARAMETER :: IONE=1
	INTEGER I,LUER
	REAL*8 T1
C
C Form "SPHERICAL" optical depth scale.
C Compute dCHI/dR for use in computation of optical depth scale.
C Compute d[dCHI/dR]/dCHI for use in linearization.
C
	DO I=1,ND
	  TA(I)=Q(I)*CHI(I)
	END DO
	CALL DERIVCHI(dCHIdR,TA,R,ND,METHOD)
	CALL d_DERIVCHI_dCHI(dCHIdR,TA,R,ND,METHOD)
	CALL NORDTAU(DTAU,TA,R,R,dCHIdR,ND)
C
C Compute T ( a tridiagonal matrix) and store it as three vectors
C TA,TB and TC .
C
	T1=HBC_J-HBC_S*THETA(1)
	CALL TFEAU(TA,TB,TC,R,Q,F,THETA,DTAU,T1,INBC,DIFF,ND)
C
C Compute WM matrix (variation of eta).
C
	FC(:,:)=0.0D0
	FC(1,1)=-HBC_S*R(1)*R(1)/CHI(1)
	DO I=2,ND-1
	  FC(I,I)=R(I)*R(I)/Q(I)/CHI(I)
	END DO
C
C Compute VK matrix (multiply's variation of CHI:-see notes)
C
	CALL VKIFEAU_IBC(F2DA,DTAU,CHI,RJ,TA,TB,TC,R,Q,F,
	1                ZETA,THETA,HBC_S,DIFF,DBB,ND)
C
C Solve the tridiagonal system of equations for the matrix illustrating
C the variation of J with ETA.
C
	CALL THOMAS(TA,TB,TC,FC,ND,ND)
	CALL SIMPTH(TA,TB,TC,F2DA,ND,ND)
C
C Compute &W vector if diffusion approximation. Note that we could
C use the FC matrix since it has two columns which are zero's
C but for simplicity we have retained a formulation consistent with
C that of PERTJD.
C
	IF(DIFF)THEN
	  CALL DP_ZERO(FA,ND)
	  FA(ND)=R(ND)*R(ND)/CHI(ND)/3.0D0
	  CALL SIMPTH(TA,TB,TC,FA,ND,IONE)
	ELSE
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error - only diffusion approximation installed'
	  STOP
	END IF
C
	RETURN
	END
