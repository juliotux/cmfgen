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
	SUBROUTINE EDD_JREL_VAR_V3(VK,RHS_dHdCHI,
	1                  R,SIGMA,CHI,ESEC,FEDD,dIBCHI_A,dIBCHI_B,
	1                  RHS_JNU,IB_STAB_FACTOR,DBB,
	1                  INNER_BND_METH,METHOD,ND,NM)
	USE MOD_VAR_JREL_V2
	IMPLICIT NONE
!
	INTEGER ND,NM
	REAL*8 VK(ND,ND,NM)
	REAL*8 RHS_dHdCHI(ND-1,ND)
	REAL*8 ESEC(ND)
	REAL*8 CHI(ND)
	REAL*8 R(ND)
	REAL*8 SIGMA(ND)
	REAL*8 FEDD(ND)
	REAL*8 DBB
	REAL*8 RHS_JNU
	REAL*8 IB_STAB_FACTOR
	REAL*8 dIBCHI_A,dIBCHI_B
	CHARACTER(LEN=*) METHOD
	CHARACTER(LEN=*) INNER_BND_METH
!
! Local variables.
!
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
	REAL*8 dUdCHI
	REAL*8 dTAdCHI_J,dTAdCHI_I
	REAL*8 dTCdCHI_I,dTCdCHI_K
	REAL*8 dTBdCHI_J,dTBdCHI_I,dTBdCHI_K,dTBdCHI
	REAL*8 dXM_EPS_J,dXM_EPS_I,dXM_EPS_K
	REAL*8 T1,T2
	INTEGER I,J,K,L
!
! 
!
	IF(NM .LT. 2)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error in EDD_JREL_VAR_V2 - NM_KI too small'
	  WRITE(I,*)'NM_KI='
	  STOP
	END IF
	VK(:,:,:)=0.0D0
	RHS_dHdCHI(:,:)=0.0D0
!
! NB: The dTAUdCHI_J & dTAUdCHI_H matrices were computed in the calling routine..
!
! The following derivatives are valid for all ML.
!
	DO I=1,ND-1
	  T1=(P_H(I)+W(I))*(CHI_H(I)+CHI_H(I+1))
	  dHUdCHI(I)=HU(I)*W(I)/T1
	  dHUdTAU(I)=-HU(I)/DTAU_H(I)
	  dHLdCHI(I)=HL(I)*W(I)/T1
	  dHLdTAU(I)=-HL(I)/DTAU_H(I)
	  dHSdCHI(I)=-P_H(I)*HS(I)/T1
	  EPS_FAC(I)=-1.0D0/T1
	END DO
!
! 
!
! Firstly we compute the variation of the elements with respect to
! DTAU. We then multiply by dTAUdCHI matrix.
!
! DTAU_H terms
!
	DO I=2,ND-1
	  J=I-1
	  K=I+1
	  dTAdCHI_J=-dHLdTAU(J)
	  dTCdCHI_I=-dHUdTAU(I)
	  dTBdCHI_I=dHLdTAU(I)
	  dTBdCHI_J=dHUdTAU(J)
	  dRHSdJ(I)=  - dTAdCHI_J*JNU(J) - dTBdCHI_J*JNU(I)
	  dRHSdI(I)=  - dTCdCHI_I*JNU(K) - dTBdCHI_I*JNU(I)
	END DO
!
	DO I=2,ND-1
	  J=I-1
	  DO L=1,ND
	    VK(I,L,1)=VK(I,L,1)+dRHSdJ(I)*dTAUdCHI_H(J,L)
	    VK(I,L,1)=VK(I,L,1)+dRHSdI(I)*dTAUdCHI_H(I,L)
	  END DO
	END DO
!
! DTAU_J terms
!
	DO I=2,ND-1
	  J=I-1
	  T1=0.5D0*R(I)*R(I)/Q(I)
	  dTBdCHI=PSI(I)/(DTAU_J(J)+DTAU_J(I)) +
	1          T1*(P_J(I)-COH_VEC(I))*GAM_REL(I)
!
! dDELUB is use as correction because UB(I)=-TB(I)-PSI(I)-PSIPREV(I)
!
	  dUdCHI=PSIPREV(I)/(DTAU_J(J)+DTAU_J(I))
	  dRHSdJ(I)= T1*SOURCE(I) - dTBdCHI*JNU(I)  + dUdCHI*JNU_PREV(I)
	  dRHSdI(I)= T1*SOURCE(I) + dUdCHI*JNU_PREV(I) - dTBdCHI*JNU(I)
!
	END DO
!
	DO I=2,ND-1
	  J=I-1
	  DO L=1,ND
	    VK(I,L,1)=VK(I,L,1)+dRHSdJ(I)*dTAUdCHI_J(J,L)
	    VK(I,L,1)=VK(I,L,1)+dRHSdI(I)*dTAUdCHI_J(I,L)
	  END DO
	END DO
!
! Can now update VK for direct opacity variation due to CHI_H.
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
	  dTBdCHI_K=dHLdCHI(I)
!
	  dXM_EPS_J=( EPS_A(J)*JNU(J)-EPS_PREV_A(J)*JNU_PREV(J) +
	1             EPS_B(J)*JNU(I)-EPS_PREV_B(J)*JNU_PREV(I) )*EPS_FAC(J)
	  dXM_EPS_K=( EPS_PREV_A(I)*JNU_PREV(I)-EPS_A(I)*JNU(I) +
	1             EPS_PREV_B(I)*JNU_PREV(I+1)-EPS_B(I)*JNU(I+1) )*EPS_FAC(I)
	  dXM_EPS_I=dXM_EPS_J+dXM_EPS_K
!
! NB  :  VB(I)=-HS(J) and VC(I)=HS(I)
!
	  VK(I,J,1)=VK(I,J,1)
	1             - (dTAdCHI_J*JNU(J)+dTBdCHI_J*JNU(I))
	1             - dHSDCHI(J)*GAM_RSQHNU_PREV(J)
	1             + dXM_EPS_J
!
	  VK(I,K,1)=VK(I,K,1)
	1             - (dTCdCHI_K*JNU(K)+dTBdCHI_K*JNU(I))
	1             + dHSDCHI(I)*GAM_RSQHNU_PREV(I)
	1             + dXM_EPS_K
!
	  VK(I,I,1)=VK(I,I,1)
	1             - (dTAdCHI_I*JNU(J)+dTCdCHI_I*JNU(K)+dTBdCHI_I*JNU(I))
	1             + (dHSDCHI(I)*GAM_RSQHNU_PREV(I)-dHSDCHI(J)*GAM_RSQHNU_PREV(J))
	1             + dXM_EPS_I
!
	END DO
!
! Can now update VK for direct opacity variation due to CHI_J.
!
	DO I=2,ND-1
	  J=I-1
	  T1=0.5D0*R(I)*R(I)/Q(I)
!
	  dTBdCHI= -PSI(I)/CHI_J(I)+GAM_RSQ_DTAUONQ(I)*(COH_VEC(I)-VdJdR_TERM(I))/CHI_J(I)
	  dUdCHI=-PSIPREV(I)/CHI_J(I)
!
	  T1=T1*(DTAU_J(J)+DTAU_J(I))/CHI_J(I)
	  VK(I,I,1)=VK(I,I,1)
	1             + (dUdCHI*JNU_PREV(I)-dTBdCHI*JNU(I))
	1             - T1*SOURCE(I)
!
	  VK(I,I,2)=T1
	END DO
!
! Now do the boundary conditions.
!
	T1=  ( (FEDD(1)+VdHdR_TERM(1))*Q(1)*JNU(1)*GAM_RSQ(1) -
	1      (FEDD(2)+VdHdR_TERM(2))*Q(2)*JNU(2)*GAM_RSQ(2) )/DTAU_H(1)/DTAU_H(1)
	DO L=1,ND
	  VK(1,L,1)=VK(1,L,1)+T1*dTAUdCHI_H(1,L)
	END DO
	VK(1,1,1)=VK(1,1,1)+ (PSI(1)*JNU(1)- PSIPREV(1)*JNU_PREV(1))/CHI_H(1)
!
! NB: We multiply DBB by GAM_REL, as we will divid VK by GAM_REL later.
! The DBB derivative is with  respect to CHI, where as the other derivatives 
! were donie with respect to CHI_J or CHI_H.
!
	IF(INNER_BND_METH .EQ. 'JEQB')THEN
	  
	ELSE IF(INNER_BND_METH .EQ. 'DIFFUSION')THEN
	  T1= ( GAM_RSQ(ND)*(FEDD(ND)+VdHdR_TERM(ND))*JNU(ND) -
	1          GAM_RSQ(ND-1)*(FEDD(ND-1)+VdHdR_TERM(ND-1))*Q(ND-1)*JNU(ND-1) )
	1           / DTAU_H(ND-1)/DTAU_H(ND-1)
	  DO L=1,ND
	    VK(ND,L,1)=VK(ND,L,1)+T1*dTAUdCHI_H(ND-1,L)
	  END DO
	  VK(ND,ND,1)=VK(ND,ND,1)-GAM_RSQ(ND)*GAM_REL(ND)*DBB/3.0D0/CHI(ND)/CHI(ND)
!
	ELSE IF(INNER_BND_METH .EQ. 'ZERO_FLUX')THEN
	  T1 = GAM_RSQ(ND-1)*(FEDD(ND-1)+VdHdR_TERM(ND-1))*Q(ND-1)
	  T2 = GAM_RSQ(ND)*(FEDD(ND)+VdHdR_TERM(ND))
	  T2 = (T2*JNU(ND)-T1*JNU(ND-1)) + IB_STAB_FACTOR*T2*(JNU(ND)-RHS_JNU)
	  T2 = T2/ DTAU_H(ND-1)/DTAU_H(ND-1)
	  DO L=1,ND
	    VK(ND,L,1)=VK(ND,L,1)+T2*dTAUdCHI_H(ND-1,L)
	  END DO
!
	ELSE IF(INNER_BND_METH .EQ. 'HOLLOW')THEN
	    VK(ND,ND,1)=VK(ND,ND,1)+dIBCHI_A
	    VK(ND,ND-1,1)=VK(ND,ND-1,1)+dIBCHI_B
!
	ELSE
	  T1= ( GAM_RSQ(ND)*(FEDD(ND)+VdHdR_TERM(ND))*JNU(ND) -
	1           GAM_RSQ(ND-1)*(FEDD(ND-1)+VdHdR_TERM(ND-1))*Q(ND-1)*JNU(ND-1) )
	1           / DTAU_H(ND-1)/DTAU_H(ND-1)
	  DO L=1,ND
	    VK(ND,L,1)=VK(ND,L,1)+T1*dTAUdCHI_H(ND-1,L)
	  END DO
	END IF
!
! Divide VK(:,:,1) by GAM_REL since dCHI_J/dCHI=dCHI_H/dCHI=1/GAM_REL
!
	DO J=1,ND
	  DO I=1,ND
	    VK(I,J,1)=VK(I,J,1)/GAM_REL(J)
	  END DO
	END DO
!
! 
!
! Now we compute the variation of the equation which updates the
! flux variation. We don't correct for the line profile, since
! we would the require two matrices. Note that HU(I), HL(I) and
! HS(I) depend directly on CHI(I) and CHI(I+1).
!
	DO I=1,ND-1
	  T1=dHUdTAU(I)*JNU(I+1)-dHLdTAU(I)*JNU(I)
	  DO L=1,ND
	    RHS_dHdCHI(I,L)=RHS_dHdCHI(I,L)+T1*dTAUdCHI_H(I,L)
	  END DO
	  T1=dHUdCHI(I)*JNU(I+1) - dHLdCHI(I)*JNU(I)
	1                       + dHSdCHI(I)*GAM_RSQHNU_PREV(I) +
	1     EPS_FAC(I)*( EPS_PREV_A(I)*JNU_PREV(I)-EPS_A(I)*JNU(I) +
	1                  EPS_PREV_B(I)*JNU_PREV(I+1)-EPS_B(I)*JNU(I+1) )
	  RHS_dHdCHI(I,I)=RHS_dHdCHI(I,I) + T1
	  RHS_dHdCHI(I,I+1)=RHS_dHdCHI(I,I+1) + T1
	END DO
!
! Recall dCHI_H/dCHI=1/GAM_REL
!
	DO L=1,ND
	  DO I=1,ND-1
	    RHS_dHdCHI(I,L)=RHS_dHdCHI(I,L)/GAM_REL(L)
	  END DO
	END DO
! 
	RETURN
	END
