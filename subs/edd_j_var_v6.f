!
! Routine to compute the opacity & emissivity variation matrices for
! the case with lines. Also computes the matrix dRHS_dCHI which
! multiply's %KI in the V equation. X is the line profile. It is assumed
! that :-
!				VK( , ,1)=dCHI
!				VK( , ,2)=dETA
!
!  				dRHS_dCHI( , ,)=dCHI
!
	SUBROUTINE EDD_J_VAR_V6(VK,RHS_dHdCHI,dTAUdCHI,
	1                  SOURCE,CHI,ESEC,ES_COH_VEC,DTAU,R,SIGMA,
	1                  MIDF,Q,HU,HL,HS,RSQ_DTAUONQ,
	1                  W,WPREV,PSI,PSIPREV,
	1                  EPS_A,EPS_B,EPS_PREV_A,EPS_PREV_B,
	1                  JNU,JNUM1,RSQ_HNUM1,
	1                  DBB,DIF,HBC,OUT_BC_TYPE,ND,NM)
	IMPLICIT NONE
!
! Altered 12-Mar-2003 : Evaluation of terms rewritten for better numerical accuracy.
! Altered 05-Feb-1998 : ES_COH_VEC now passed. Logical variable coherent
!                         no longer required. Introduced to allow
!                         coherent scattering at some depths, and
!                         incoherent scattering at other depths.
!                         Equivalent to THETA in many programs.
!                         Changed to version 5.
! Altered 16-Aug-1996 : Loop splitting to improve vectorization.
! Altered 24-May-1996 : NV removed (Max valued for ND) using F90
!                       CALL to DP_ZERO replaced.
! Created 16-Sep-1994 : Based on EDDLINE_VAR
!
	INTEGER ND,NM
	REAL*8 VK(ND,ND,NM),RHS_dHdCHI(ND-1,ND)
	REAL*8 dTAUdCHI(ND,ND)
	REAL*8 SOURCE(ND),CHI(ND),ESEC(ND),ES_COH_VEC(ND)
	REAL*8 DTAU(ND),R(ND),SIGMA(ND)
	REAL*8 MIDF(ND),Q(ND),HU(ND),HL(ND),HS(ND),RSQ_DTAUONQ(ND)
	REAL*8 W(ND),WPREV(ND),PSI(ND),PSIPREV(ND)
	REAL*8 EPS_A(ND),EPS_B(ND)
	REAL*8 EPS_PREV_A(ND),EPS_PREV_B(ND)
	REAL*8 JNU(ND),JNUM1(ND),RSQ_HNUM1(ND)
	REAL*8 DBB
	REAL*8 HBC
	LOGICAL DIF
	INTEGER OUT_BC_TYPE
!
! Local variables.
!
	REAL*8 dHUdCHI(ND),dHLdCHI(ND),dHSdCHI(ND)
	REAL*8 dHUdTAU(ND),dHLdTAU(ND),EPS_FAC(ND)
	REAL*8 dRHSdI(ND),dRHSdJ(ND)
!
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
	INTEGER I,J,K,L
	REAL*8 T1,T2,T3
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
	  WRITE(I,*)'Error in EDD_J_VAR_V6 - NM_KI too small'
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
	  T1=0.5D0*R(I)*R(I)/Q(I)
	  dTBdCHI_I=dHLdTAU(I)
	  dTBdCHI_J=dHUdTAU(J)
	  dTBdCHI=PSI(I)/(DTAU(J)+DTAU(I))+0.5D0*(1.0D0-ES_COH_VEC(I))*R(I)*R(I)/Q(I)
!
! dDELUB is use as correction because UB(I)=-TB(I)-PSI(I)-PSIPREV(I)
!
	  dUdCHI=PSIPREV(I)/(DTAU(J)+DTAU(I))
!
	  dRHSdJ(I)=  T1*SOURCE(I)
	1         - (dTAdCHI_J*JNU(J)+dTBdCHI_J*JNU(I))
	1         - dTBdCHI*JNU(I)
	1         + dUdCHI*JNUM1(I)
!
	  dRHSdI(I)= T1*SOURCE(I)
	1         - (dTCdCHI_I*JNU(K)+dTBdCHI_I*JNU(I))
	1         + (dUdCHI*JNUM1(I)-dTBdCHI*JNU(I))
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
	  dTBdCHI= -PSI(I)/CHI(I)+RSQ_DTAUONQ(I)*ES_COH_VEC(I)/CHI(I)
	  dTBdCHI_K=dHLdCHI(I)
!
	  dUdCHI=-PSIPREV(I)/CHI(I)
!
	  dXM_EPS_J=( EPS_A(J)*JNU(J)-EPS_PREV_A(J)*JNUM1(J) +
	1             EPS_B(J)*JNU(I)-EPS_PREV_B(J)*JNUM1(I) )*EPS_FAC(J)
	  dXM_EPS_K=( EPS_PREV_A(I)*JNUM1(I)-EPS_A(I)*JNU(I) +
	1             EPS_PREV_B(I)*JNUM1(I+1)-EPS_B(I)*JNU(I+1) )*EPS_FAC(I)
	  dXM_EPS_I=dXM_EPS_J+dXM_EPS_K
!
! NB  :  VB(I)=-HS(J) and VC(I)=HS(I)
!
	  VK(I,J,1)=VK(I,J,1)
	1             - (dTAdCHI_J*JNU(J)+dTBdCHI_J*JNU(I))
	1             - dHSDCHI(J)*RSQ_HNUM1(J)
	1             + dXM_EPS_J
!
	  VK(I,K,1)=VK(I,K,1)
	1             - (dTCdCHI_K*JNU(K)+dTBdCHI_K*JNU(I))
	1             + dHSDCHI(I)*RSQ_HNUM1(I)
	1             + dXM_EPS_K
!
	  T1=T1*(DTAU(J)+DTAU(I))/CHI(I)
	  VK(I,I,1)=VK(I,I,1)
	1             - (dTAdCHI_I*JNU(J)+dTCdCHI_I*JNU(K)+dTBdCHI_I*JNU(I))
	1             + (dUdCHI*JNUM1(I)-dTBdCHI*JNU(I))
	1             + (dHSDCHI(I)*RSQ_HNUM1(I)-dHSDCHI(J)*RSQ_HNUM1(J))
	1             - T1*SOURCE(I)
	1             + dXM_EPS_I
!
	  VK(I,I,2)=T1
	END DO
!
! Now do the boundary conditions.
!
	IF(OUT_BC_TYPE .LE. 1)THEN
	  T1=  ( MIDF(1)*Q(1)*JNU(1)*R(1)*R(1)
	1        - MIDF(2)*Q(2)*JNU(2)*R(2)*R(2) )/DTAU(1)/DTAU(1)
	  DO L=1,ND
	    VK(1,L,1)=VK(1,L,1)+T1*dTAUdCHI(1,L)
	  END DO
	  VK(1,1,1)=VK(1,1,1) +
	1          ( PSI(1)*JNU(1)- PSIPREV(1)*JNUM1(1) )/CHI(1)
	ELSE
!
! First do the variation with respect to CHI arising from the 0.5*dR*(CHI(1)+CHI(2)) term.
!
	  T1=0.5D0*(R(2)-R(1))
	  T2=T1*(CHI(2)+CHI(1))
	  T3=( HS(1)*RSQ_HNUM1(1) + (HU(1)*JNU(2)-(HL(1)-HBC*R(1)*R(1))*JNU(1)) +
	1       EPS_A(1)*(JNUM1(1)-JNU(1)) +
	1       EPS_B(1)*(JNUM1(2)-JNU(2)) )*T1/T2/T2
	  VK(1,1,1)=VK(1,1,1) + T3
	  VK(1,2,1)=VK(1,1,1) + T3
!
! We now do the variation arising from DTAU(1) (which contains a q factor)
!
	  T3= -(dHUdTAU(1)*JNU(2)-dHLdTAU(1)*JNU(1))/T2
	  VK(1,1,1)=VK(1,1,1) + T3*dTAUdCHI(1,1)
	  VK(1,2,1)=VK(1,2,1) + T3*dTAUdCHI(1,2)
	  VK(1,3,1)=VK(1,3,1) + T3*dTAUdCHI(1,3)
!
! Now do the variation arising directly from terms containing 1/CHI(1).
!
	  T3= R(1)*R(1)*(SOURCE(1)+ES_COH_VEC(1)*JNU(1))/CHI(1) +
	1           (PSIPREV(1)*JNUM1(1)-PSI(1)*JNU(1))/CHI(1)
	  VK(1,1,1)=VK(1,1,1) + T3
!
! Now do terms containing 1/(CHI(1)+CHI(2))
!
	  T3=(dHLdCHI(1)*JNU(1)- dHUdCHI(1)*JNU(2)-dHSdCHI(1)*RSQ_HNUM1(1))/T2
	1       + ( EPS_A(1)*(JNU(1)-JNUM1(1))+EPS_B(1)*(JNU(2)-JNUM1(2)) )*EPS_FAC(1)/T2
	  VK(1,1,1)=VK(1,1,1) + T3
	  VK(1,2,1)=VK(1,2,1) + T3
!
! Now for the variation with ETA.
!
	  VK(1,1,2)=-R(1)*R(1)/CHI(1)
!
	END IF
!
	IF(DIF)THEN
	  T1= ( R(ND)*R(ND)*MIDF(ND)*JNU(ND) -
	1          R(ND-1)*R(ND-1)*MIDF(ND-1)*Q(ND-1)*JNU(ND-1) )
	1           / DTAU(ND-1)/DTAU(ND-1)
	  DO L=1,ND
	    VK(ND,L,1)=VK(ND,L,1)+T1*dTAUdCHI(ND-1,L)
	  END DO
	  VK(ND,ND,1)=VK(ND,ND,1)-
	1               DBB*R(ND)*R(ND)/3.0D0/CHI(ND)/CHI(ND)
	ELSE
	  T1= ( R(ND)*R(ND)*MIDF(ND)*JNU(ND) -
	1           R(ND-1)*R(ND-1)*MIDF(ND-1)*Q(ND-1)*JNU(ND-1) )
	1           / DTAU(ND-1)/DTAU(ND-1)
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
	1                       + dHSdCHI(I)*RSQ_HNUM1(I) +
	1     EPS_FAC(I)*( EPS_PREV_A(I)*JNUM1(I)-EPS_A(I)*JNU(I) +
	1                  EPS_PREV_B(I)*JNUM1(I+1)-EPS_B(I)*JNU(I+1) )
	  RHS_dHdCHI(I,I)=RHS_dHdCHI(I,I) + T1
	  RHS_dHdCHI(I,I+1)=RHS_dHdCHI(I,I+1) + T1
	END DO
!
	RETURN
	END
