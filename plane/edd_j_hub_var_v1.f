!
! Routine to compute the opacity & emissivity variation matrices for
! the case with lines. Also computes the matrix dRHS_dCHI which
! multiply's %KI in the V equation. X is the line profile. It is assumed
! that :-
!				VK( , ,1)=dCHI
!				VK( , ,2)=dETA
!  				dRHS_dHdCHI( , ,)  !RHS of H equation.
!
	SUBROUTINE EDD_J_HUB_VAR_V1(VK,RHS_dHdCHI,dTAUdCHI,
	1                  SOURCE,CHI,ESEC,ES_COH_VEC,DTAU,R,
	1                  EDDF,Q,HU,HL,HS,HT,RSQ_DTAUONQ,
	1                  W,WPREV,PSI,PSIPREV,DJDt,DJDt_OLDT,
	1                  JNU,JNUM1,JNU_OLDT,
	1                  RSQ_HNUM1,RSQ_HNU_OLDT,
	1                  dRHSdCHI_DIF_BC,DIF,ND,NM)
	IMPLICIT NONE
!
! Altered 25-Apr-2009 : Bug fig: dHS replaced by dHT in one expression.
! Altered 21-Feb-2007 : DBB replaced by dRHSdCHI_DIF_BC. Fixes small bug in EDD_J_HUB_VAR_V1.
! Created 19-July-2006
!
	INTEGER ND,NM
	REAL*8 VK(ND,ND,NM)
	REAL*8 RHS_dHdCHI(ND-1,ND)
	REAL*8 dTAUdCHI(ND,ND)
	REAL*8 SOURCE(ND)
	REAL*8 CHI(ND)
	REAL*8 ESEC(ND)
	REAL*8 ES_COH_VEC(ND)
	REAL*8 DTAU(ND)
	REAL*8 R(ND)
	REAL*8 EDDF(ND)
	REAL*8 Q(ND)
	REAL*8 HU(ND),HL(ND),HS(ND),HT(ND)
	REAL*8 RSQ_DTAUONQ(ND)
	REAL*8 W(ND),WPREV(ND)
	REAL*8 DJDt(ND),DJDt_OLDT(ND)
	REAL*8 PSI(ND),PSIPREV(ND)
	REAL*8 JNU(ND),JNUM1(ND),JNU_OLDT(ND)
	REAL*8 RSQ_HNUM1(ND),RSQ_HNU_OLDT(ND)
	REAL*8 dRHSdCHI_DIF_BC
	LOGICAL DIF
!
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
!
! Local vectors.
!
	REAL*8 dHUdCHI(ND),dHLdCHI(ND)
	REAL*8 dHSdCHI(ND),dHTdCHI(ND)
	REAL*8 dHUdTAU(ND),dHLdTAU(ND)
	REAL*8 dRHSdI(ND),dRHSdJ(ND)
!
! Local varoables.
!
	INTEGER I,J,K,L
	REAL*8 T1
	REAL*8 dUdCHI
	REAL*8 dTAdCHI_J,dTAdCHI_I
	REAL*8 dTCdCHI_I,dTCdCHI_K
	REAL*8 dTBdCHI_J,dTBdCHI_I,dTBdCHI_K,dTBdCHI
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
	  dHTdCHI(I)=-HT(I)/T1
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
	  T1=0.5D0*R(I)*R(I)/Q(I)
	  dTBdCHI=(PSI(I)+DJDT(I))/(DTAU(J)+DTAU(I))+T1*(1.0D0-ES_COH_VEC(I))
!
	  dUdCHI=(PSIPREV(I)*JNUM1(I)+DJDT_OLDT(I)*JNU_OLDT(I))/(DTAU(J)+DTAU(I))
	  dRHSdJ(I)=  T1*SOURCE(I)
	1         - (dTAdCHI_J*JNU(J)+dTBdCHI_J*JNU(I))
	1         - dTBdCHI*JNU(I)
	1         + dUdCHI
!
	  dRHSdI(I)= T1*SOURCE(I)
	1         - (dTCdCHI_I*JNU(K)+dTBdCHI_I*JNU(I))
	1         + dUdCHI-dTBdCHI*JNU(I)
!
	END DO
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
	  T1=0.5D0*R(I)*R(I)/Q(I)
!
	  dTAdCHI_J=-dHLdCHI(J)
	  dTAdCHI_I=-dHLdCHI(J)
	  dTCdCHI_I=-dHUdCHI(I)
	  dTCdCHI_K=-dHUdCHI(I)
!
	  dTBdCHI_J=dHUdCHI(J)
	  dTBdCHI_I=dHLdCHI(I)+dHUdCHI(J)
	  dTBdCHI= RSQ_DTAUONQ(I)*ES_COH_VEC(I)/CHI(I)-(PSI(I)+DJDT(I))/CHI(I)
	  dTBdCHI_K=dHLdCHI(I)
!
! Note: We have changed sign of dUdCHI reative to an earlier version.
!
	  dUdCHI=(PSIPREV(I)*JNUM1(I)+DJDT_OLDT(I)*JNU_OLDT(I))/CHI(I)
!
! NB  :  VB(I)=-HS(J) and VC(I)=HS(I)
!
	  VK(I,J,1)=VK(I,J,1)
	1             - (dTAdCHI_J*JNU(J)+dTBdCHI_J*JNU(I))
	1             - dHSdCHI(J)*RSQ_HNUM1(J)
	1             - dHTdCHI(J)*RSQ_HNU_OLDT(J)
!
	  VK(I,K,1)=VK(I,K,1)
	1             - (dTCdCHI_K*JNU(K)+dTBdCHI_K*JNU(I))
	1             + dHSDCHI(I)*RSQ_HNUM1(I)
	1             + dHTdCHI(I)*RSQ_HNU_OLDT(I)
!
	  T1=T1*(DTAU(J)+DTAU(I))/CHI(I)
	  VK(I,I,1)=VK(I,I,1)
	1             - (dTAdCHI_I*JNU(J)+dTCdCHI_I*JNU(K)+dTBdCHI_I*JNU(I))
	1             - (dUdCHI+dTBdCHI*JNU(I))
	1             + (dHSDCHI(I)*RSQ_HNUM1(I)-dHSDCHI(J)*RSQ_HNUM1(J))
	1             + (dHTDCHI(I)*RSQ_HNU_OLDT(I)-dHTDCHI(J)*RSQ_HNU_OLDT(J))
	1             - T1*SOURCE(I)
!
	  VK(I,I,2)=T1
	END DO
!
! Now do the boundary conditions.
!
	T1=  ( EDDF(1)*Q(1)*JNU(1)*R(1)*R(1)
	1      - EDDF(2)*Q(2)*JNU(2)*R(2)*R(2) )/DTAU(1)/DTAU(1)
	DO L=1,ND
	  VK(1,L,1)=VK(1,L,1)+T1*dTAUdCHI(1,L)
	END DO
	VK(1,1,1)=VK(1,1,1) +
	1          ( PSI(1)*JNU(1)- PSIPREV(1)*JNUM1(1) )/CHI(1) +
	1          ( DJDt(1)*JNU(1)- DJDt_OLDT(1)*JNU_OLDT(1) )/CHI(1)
!
! Inner boundary --- diffusion approximation.
!
	IF(DIF)THEN
	  T1= ( R(ND)*R(ND)*EDDF(ND)*JNU(ND) -
	1          R(ND-1)*R(ND-1)*EDDF(ND-1)*Q(ND-1)*JNU(ND-1) )
	1           / DTAU(ND-1)/DTAU(ND-1)
	  DO L=1,ND
	    VK(ND,L,1)=VK(ND,L,1)+T1*dTAUdCHI(ND-1,L)
	  END DO
	  VK(ND,ND,1)=VK(ND,ND,1)+dRHSdCHI_DIF_BC
	ELSE
	  I=ERROR_LU()
	  WRITE(I,*)'Non-diffusion boundary condition not yet implemented'
	  WRITE(I,*)'Error occured in EDD_J_HIB_VAR_V1'
	  STOP
	END IF
!
! 
!
! We compute the variation of the equation which updates the
! flux variation. We don't correct for the line profile, since
! we would the require two matrices. Note that HU(I), HL(I),
! HS(I), and HT(I) depend directly on CHI(I) and CHI(I+1).
!
	DO I=1,ND-1
	  T1=dHUdTAU(I)*JNU(I+1)-dHLdTAU(I)*JNU(I)
	  DO L=1,ND
	    RHS_dHdCHI(I,L)=RHS_dHdCHI(I,L)+T1*dTAUdCHI(I,L)
	  END DO
	  T1=dHUdCHI(I)*JNU(I+1) - dHLdCHI(I)*JNU(I)
	1                       + dHSdCHI(I)*RSQ_HNUM1(I)
	1                       + dHTdCHI(I)*RSQ_HNU_OLDT(I)
	  RHS_dHdCHI(I,I)=RHS_dHdCHI(I,I) + T1
	  RHS_dHdCHI(I,I+1)=RHS_dHdCHI(I,I+1) + T1
	END DO
!
	RETURN
	END
