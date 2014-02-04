!
! Routine to compute the opacity & emissivity variation matrices for
! the case with lines. Also computes the matrix dRHS_dCHI which
! multiply's %KI in the V equation. X is the line profile. It is assumed
! that :-
!				VK( , ,1)=dCHI
!				VK( , ,2)=dETA
!  				dRHS_dCHI( , ,)=dCHI
!
! Routines is designed for a plane-parellel atmosphere with a
! velocity field.
!
	SUBROUTINE PP_EDD_VAR_CMF_V1(VK,RHS_dHdCHI,dTAUdCHI,
	1                  SOURCE,CHI,ESEC,ES_COH_VEC,DTAU,R,SIGMA,
	1                  MIDF,HU,HL,HS,MID_DTAU,
	1                  W,WPREV,PSI,PSIPREV,EPS,EPS_PREV,
	1                  JNU,JNUM1,HNUM1,
	1                  DBB,DIF,ND,NM)
	IMPLICIT NONE
!
! Created 19-Mar-2006 : Based on EDD_J_VAR_V5
!
	INTEGER ND,NM
	REAL*8 VK(ND,ND,NM)
	REAL*8 RHS_dHdCHI(ND-1,ND)
	REAL*8 dTAUdCHI(ND,ND)
	REAL*8 SOURCE(ND),CHI(ND),ESEC(ND),ES_COH_VEC(ND)
	REAL*8 DTAU(ND),R(ND),SIGMA(ND)
	REAL*8 MIDF(ND),HU(ND),HL(ND),HS(ND),MID_DTAU(ND)
	REAL*8 W(ND),WPREV(ND),PSI(ND),PSIPREV(ND)
	REAL*8 EPS(ND)
	REAL*8 EPS_PREV(ND)
	REAL*8 JNU(ND),JNUM1(ND),HNUM1(ND)
	REAL*8 DBB
	LOGICAL DIF
!
! Local variables.
!
	REAL*8 dHUdCHI(ND),dHLdCHI(ND),dHSdCHI(ND)
	REAL*8 dHUdTAU(ND),dHLdTAU(ND),EPS_FAC(ND)
	REAL*8 dRHSdI(ND),dRHSdJ(ND)
	REAL*8 Q(ND)
!
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
	INTEGER I,J,K,L
	REAL*8 T1
	REAL*8 dUdCHI
	REAL*8 dTAdCHI_J,dTAdCHI_I
	REAL*8 dTCdCHI_I,dTCdCHI_K
	REAL*8 dTBdCHI_J,dTBdCHI_I,dTBdCHI_K,dTBdCHI
	REAL*8 dXM_EPS_J,dXM_EPS_I,dXM_EPS_K
!
! 
!
	IF(NM .LT. 2)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error in EDD_J_VAR_V4 - NM_KI too small'
	  WRITE(I,*)'NM_KI='
	  STOP
	END IF
	VK(:,:,:)=0.0D0
	RHS_dHdCHI(:,:)=0.0D0
!
! Compute the dTAUdCHI matrix.
!
	Q(1:ND)=1.0D0
	CALL dSPHEREdCHI(dTAUdCHI,DTAU,R,Q,ND)
!
! The following derivatives are valid for all ML.
!
	DO I=1,ND-1
	  T1=(1.0D0+W(I))*(CHI(I)+CHI(I+1))
	  dHUdCHI(I)=HU(I)*W(I)/T1
	  dHUdTAU(I)=-HU(I)/DTAU(I)
	  dHLdCHI(I)=HL(I)*W(I)/T1
	  dHLdTAU(I)=-HL(I)/DTAU(I)
	  dHSdCHI(I)=-HS(I)/T1
	  EPS_FAC(I)=-1.0D0/T1
	END DO
!
! 
!
! Firstly we compute the variation of the elements with respect to
! DTAU. We then multiply by dTAUdCHI matrix.
!
! We have 3 separate loops over I to allow vectorization.
! To improve rounding error, we note that dTBdCHI needs to be added to
! both dTBdCHI_I and dTBdCHI_J.
!
	DO I=2,ND-1
	  J=I-1
	  K=I+1
	  dTAdCHI_J=-dHLdTAU(J)
	  dTCdCHI_I=-dHUdTAU(I)
	  dTBdCHI_I=dHLdTAU(I)
	  dTBdCHI_J=dHUdTAU(J)
	  dTBdCHI=PSI(I)/(DTAU(J)+DTAU(I))+0.5D0*(1.0D0-ES_COH_VEC(I))
!
! dDELUB is use as correction because UB(I)=-TB(I)-PSI(I)-PSIPREV(I)
!
	  dUdCHI=PSIPREV(I)/(DTAU(J)+DTAU(I))
!
	  dRHSdJ(I)=  0.5D0*SOURCE(I)
	1         - (dTAdCHI_J*JNU(J)+dTBdCHI_J*JNU(I))
	1         - dTBdCHI*JNU(I)
	1         + dUdCHI*JNUM1(I)
!
	  dRHSdI(I)= 0.5D0*SOURCE(I)
	1         - (dTCdCHI_I*JNU(K)+dTBdCHI_I*JNU(I))
	1         + (dUdCHI*JNUM1(I)-dTBdCHI*JNU(I))
!
	END DO
!
! At present we sum over all depths. This could be speed up by summing
! only over the diaginal and adjacent depths.
!
	DO I=2,ND-1
	  J=I-1
	  K=I+1
	  DO L=1,ND
	    VK(I,L,1)=VK(I,L,1)+dRHSdJ(I)*dTAUdCHI(J,L)
	    VK(I,L,1)=VK(I,L,1)+dRHSdI(I)*dTAUdCHI(I,L)
	  END DO
	END DO
!
! Can now update VK for direct opacity variation.
!
	DO I=2,ND-1
	  J=I-1
	  K=I+1
!
	  dTAdCHI_J=-dHLdCHI(J)
	  dTAdCHI_I=-dHLdCHI(J)
	  dTCdCHI_I=-dHUdCHI(I)
	  dTCdCHI_K=-dHUdCHI(I)
!
	  dTBdCHI_J=dHUdCHI(J)
	  dTBdCHI_I=dHLdCHI(I)+dHUdCHI(J)
	  dTBdCHI= -PSI(I)/CHI(I)+MID_DTAU(I)*ES_COH_VEC(I)/CHI(I)
	  dTBdCHI_K=dHLdCHI(I)
!
	  dUdCHI=-PSIPREV(I)/CHI(I)
!
	  dXM_EPS_J=( EPS(J)*JNU(J)-EPS_PREV(J)*JNUM1(J) +
	1             EPS(J)*JNU(I)-EPS_PREV(J)*JNUM1(I) )*EPS_FAC(J)
	  dXM_EPS_K=( EPS_PREV(I)*JNUM1(I)-EPS(I)*JNU(I) +
	1             EPS_PREV(I)*JNUM1(I+1)-EPS(I)*JNU(I+1) )*EPS_FAC(I)
	  dXM_EPS_I=dXM_EPS_J+dXM_EPS_K
!
! NB  :  VB(I)=-HS(J) and VC(I)=HS(I)
!
	  VK(I,J,1)=VK(I,J,1)
	1             - (dTAdCHI_J*JNU(J)+dTBdCHI_J*JNU(I))
	1             - dHSDCHI(J)*HNUM1(J)
	1             + dXM_EPS_J
!
	  VK(I,K,1)=VK(I,K,1)
	1             - (dTCdCHI_K*JNU(K)+dTBdCHI_K*JNU(I))
	1             + dHSDCHI(I)*HNUM1(I)
	1             + dXM_EPS_K
!
	  T1=0.5D0*(DTAU(J)+DTAU(I))/CHI(I)
	  VK(I,I,1)=VK(I,I,1)
	1             - (dTAdCHI_I*JNU(J)+dTCdCHI_I*JNU(K)+dTBdCHI_I*JNU(I))
	1             + (dUdCHI*JNUM1(I)-dTBdCHI*JNU(I))
	1             + (dHSDCHI(I)*HNUM1(I)-dHSDCHI(J)*HNUM1(J))
	1             - T1*SOURCE(I)
	1             + dXM_EPS_I
!
	  VK(I,I,2)=T1
	END DO
!
! Now do the boundary conditions.
!
	T1=(MIDF(1)*JNU(1)-MIDF(2)*JNU(2))/DTAU(1)/DTAU(1)
	DO L=1,ND
	  VK(1,L,1)=VK(1,L,1)+T1*dTAUdCHI(1,L)
	END DO
	VK(1,1,1)=VK(1,1,1)+(PSI(1)*JNU(1)-PSIPREV(1)*JNUM1(1))/CHI(1)
!
	IF(DIF)THEN
	  T1= (MIDF(ND)*JNU(ND)-MIDF(ND-1)*JNU(ND-1))/DTAU(ND-1)/DTAU(ND-1)
	  DO L=1,ND
	    VK(ND,L,1)=VK(ND,L,1)+T1*dTAUdCHI(ND-1,L)
	  END DO
	  VK(ND,ND,1)=VK(ND,ND,1)-DBB/3.0D0/CHI(ND)/CHI(ND)
	ELSE
	  T1= (MIDF(ND)*JNU(ND)-MIDF(ND-1)*JNU(ND-1))/DTAU(ND-1)/DTAU(ND-1)
	  DO L=1,ND
	    VK(ND,L,1)=VK(ND,L,1)+T1*dTAUdCHI(ND-1,L)
	  END DO
	END IF
!
! 
!
! Now we compute the varaition of the equation which updates the
! flux variation. We don't correct for the line profile, since
! we would the require two matrices. Note that HU(I), HL(I) and
! HS(I) depend directly on CHI(I) and CHI(I+1).
!
	DO I=1,ND-1
	  T1=dHUdTAU(I)*JNU(I+1)-dHLdTAU(I)*JNU(I)
	  DO L=1,ND
	    RHS_dHdCHI(I,L)=RHS_dHdCHI(I,L)+T1*dTAUdCHI(I,L)
	  END DO
	  T1=dHUdCHI(I)*JNU(I+1) - dHLdCHI(I)*JNU(I)
	1                       + dHSdCHI(I)*HNUM1(I) +
	1     EPS_FAC(I)*( EPS_PREV(I)*JNUM1(I)-EPS(I)*JNU(I) +
	1                  EPS_PREV(I)*JNUM1(I+1)-EPS(I)*JNU(I+1) )
	  RHS_dHdCHI(I,I)=RHS_dHdCHI(I,I) + T1
	  RHS_dHdCHI(I,I+1)=RHS_dHdCHI(I,I+1) + T1
	END DO
!
	RETURN
	END
