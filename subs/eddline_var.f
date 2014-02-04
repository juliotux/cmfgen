C
C Routine to compute the opacity & emissivity variation matrices for
C the case with lines. Also computes the matrix dRHS_dCHI which
C multiply's %KI in the V equation. X is the line profile. It is assumed
C that :-
C				VK( , ,1)=dCHIL
C				VK( , ,2)=dETAL
C				VK( , ,3)=dCHI (if needed)
C				VK( , ,4)=dETA (if needed)
C
C  				dRHS_dCHI( , ,)=dTCHI
C
	SUBROUTINE EDDLINE_VAR(VK,RHS_dHdCHI,
	1                  dTAUdCHI,
	1                  SOURCE,TCHI,DTAU,R,SIGMA,
	1                  MIDF,Q,HU,HL,HS,
	1                  W,WPREV,PSI,PSIPREV,
	1                  AV,AVM1,CVM1,
	1                  DEPTHPRO,DBB,DIF,
	1                  ML,NLF,ND,NM)
	IMPLICIT NONE
C
C Altered 24-May-1996 - NV removed using F90 (was max value for ND)
C Altered  2-Jun-1989 - d(HU,HL)dCHI derivative sign fixed.
C Altered 15-May-1989 - HS derivative fixed.
C Created 09-May-1989
C
	INTEGER ML,NLF,ND,NM
	REAL*8 VK(ND,ND,NM),RHS_dHdCHI(ND-1,ND)
	REAL*8 dTAUdCHI(ND,ND)
	REAL*8 SOURCE(ND),TCHI(ND),DTAU(ND),R(ND),SIGMA(ND)
	REAL*8 MIDF(ND),Q(ND),HU(ND),HL(ND),HS(ND)
	REAL*8 W(ND),WPREV(ND),PSI(ND),PSIPREV(ND)
	REAL*8 AV(ND),AVM1(ND),CVM1(ND)
	REAL*8 DEPTHPRO(ND)
	REAL*8 DBB
	LOGICAL DIF
C
C Local variables.
C
	REAL*8 dHUdCHI(ND),dHLdCHI(ND),dHSdCHI(ND)
	REAL*8 dHUdTAU(ND),dHLdTAU(ND)
C
	INTEGER I,J,K,L
	REAL*8 T1
	REAL*8 dRHSdI,dRHSdJ,dUdCHI
	REAL*8 dTAdCHI_J,dTAdCHI_I
	REAL*8 dTCdCHI_I,dTCdCHI_K
	REAL*8 dTBdCHI_J,dTBdCHI_I,dTBdCHI_K
C
C 
C
	VK(:,:,:)=0.0D0 	!Dimension ND,ND,NM
	RHS_dHdCHI(:,:)=0.0D0	!Dimensiond ND-1,ND
C
C Compute the dTAUdCHI matrix.
C
	CALL dSPHEREdCHI(dTAUdCHI,DTAU,R,Q,ND)
C
C The following derivatives are valid for all ML.
C
	DO I=1,ND-1
	  T1=(1.0D0+W(I))*(TCHI(I)+TCHI(I+1))
	  dHUdCHI(I)=HU(I)*W(I)/T1
	  dHUdTAU(I)=-HU(I)/DTAU(I)
	  dHLdCHI(I)=HL(I)*W(I)/T1
	  dHLdTAU(I)=-HL(I)/DTAU(I)
	  dHSdCHI(I)=-HS(I)/T1
	END DO
C
C 
C
C Handle the special case of ML=1. For this special case, the dCHIL
C and dETAL derivatives are zero. Also AVM1 and CVM1 are zero.
C
	IF( ML .EQ. 1 )THEN
	  IF(NM .NE. 4)RETURN
C
C Firstly we compute the variation of the elements with respect to
C DTAU. We then multiply by dTAUdCHI matrix.
C
	  DO I=2,ND-1
	    J=I-1
	    K=I+1
	    dTAdCHI_J=-dHLdTAU(J)
	    dTCdCHI_I=-dHUdTAU(I)
	    dTBdCHI_I=dHLdTAU(I)+0.5D0/Q(I)
	    dTBdCHI_J=dHUdTAU(J)+0.5D0/Q(I)
	    T1=0.5D0*R(I)*R(I)/Q(I)
C
	    dRHSdJ= T1*SOURCE(I) - dTAdCHI_J*AV(J) - dTBdCHI_J*AV(I)
	    dRHSdI= T1*SOURCE(I) - dTCdCHI_I*AV(K) - dTBdCHI_I*AV(I)
C
	    DO L=1,ND
	      VK(I,L,3)=VK(I,L,3)+dRHSdJ*dTAUdCHI(J,L)
	      VK(I,L,3)=VK(I,L,3)+dRHSdI*dTAUdCHI(I,L)
	    END DO
C
	    T1=T1*(DTAU(J)+DTAU(I))/TCHI(I)
	    VK(I,I,3)=VK(I,I,3) - T1*SOURCE(I)
	    VK(I,I,4)=T1
C
	  END DO
C
C Now do the boundary conditions.
C
	  T1= ( MIDF(1)*Q(1)*AV(1) -
	1         MIDF(2)*Q(2)*AV(2) )/DTAU(1)/DTAU(1)
	  DO L=1,ND
	    VK(1,L,3)=VK(1,L,3)+T1*dTAUdCHI(1,L)
	  END DO
C
	  IF(DIF)THEN
	    T1=  ( MIDF(ND)*AV(ND) - MIDF(ND-1)*Q(ND-1)*AV(ND-1) )
	1            /DTAU(ND-1)/DTAU(ND-1)
	    DO L=1,ND
	      VK(ND,L,3)=VK(ND,L,3)+T1*dTAUdCHI(ND-1,L)
	    END DO
	    VK(ND,ND,3)=VK(ND,ND,3)-
	1                 DBB*R(ND)*R(ND)/3.0D0/TCHI(ND)/TCHI(ND)
	  ELSE
	    T1=  ( MIDF(ND)*AV(ND) - MIDF(ND-1)*Q(ND-1)*AV(ND-1) )
	1            /DTAU(ND-1)/DTAU(ND-1)
	    DO L=1,ND
	      VK(ND,L,3)=VK(ND,L,3)+T1*dTAUdCHI(ND-1,L)
	    END DO
	  END IF
C
C Now we compute the varaition of the equation which updates the
C flux variation. We don't correct for the line profile, since
C we would the require two matrices.
C
	  DO I=1,ND-1
	    T1=dHUdTAU(I)*AV(I+1)-dHLdTAU(I)*AV(I)
	    DO L=1,ND
	      RHS_dHdCHI(I,L)=RHS_dHdCHI(I,L)+T1*dTAUdCHI(I,L)
	    END DO
	  END DO
C
	  RETURN
	END IF
C
C 
C
C Firstly we compute the variation of the elements with respect to
C DTAU. We then multiply by dTAUdCHI matrix.
C
	DO I=2,ND-1
	 J=I-1
	 K=I+1
	 dTAdCHI_J=-dHLdTAU(J)
	 dTCdCHI_I=-dHUdTAU(I)
	 dTBdCHI_I=dHLdTAU(I)+PSI(I)/(DTAU(J)+DTAU(I))+0.5D0/Q(I)
	 dTBdCHI_J=dHUdTAU(J)+PSI(I)/(DTAU(J)+DTAU(I))+0.5D0/Q(I)
	 T1=0.5D0*R(I)*R(I)/Q(I)
C
C dDELUB is use as correction because UB(I)=-TB(I)-PSI(I)-PSIPREV(I)
C
	 dUdCHI=PSIPREV(I)/(DTAU(J)+DTAU(I))
C
	 dRHSdJ=  T1*SOURCE(I)
	1         - dTAdCHI_J*AV(J)-dTBdCHI_J*AV(I)
	1         + dUdCHI*AVM1(I)
C
	 dRHSdI= T1*SOURCE(I)
	1         - dTCdCHI_I*AV(K)-dTBdCHI_I*AV(I)
	1         + dUdCHI*AVM1(I)
C
	DO L=1,ND
	  VK(I,L,1)=VK(I,L,1)+dRHSdJ*dTAUdCHI(J,L)
	  VK(I,L,1)=VK(I,L,1)+dRHSdI*dTAUdCHI(I,L)
	END DO
C
C Can now update VK for direct opacity variation.
C
	 dTAdCHI_J=-dHLdCHI(J)
	 dTAdCHI_I=-dHLdCHI(J)
	 dTCdCHI_I=-dHUdCHI(I)
	 dTCdCHI_K=-dHUdCHI(I)
C
	 dTBdCHI_J=dHUdCHI(J)
	 dTBdCHI_I=dHLdCHI(I)+dHUdCHI(J)-PSI(I)/TCHI(I)
	 dTBdCHI_K=dHLdCHI(I)
C
	 dUdCHI=-PSIPREV(I)/TCHI(I)
C
C NB  :  VB(I)=-HS(J) and VC(I)=HS(I)
C
	 VK(I,J,1)=VK(I,J,1)
	1             - dTAdCHI_J*AV(J)
	1             - dTBdCHI_J*AV(I)
	1             - dHSDCHI(J)*CVM1(J)
C
	 VK(I,K,1)=VK(I,K,1)
	1             - dTCdCHI_K*AV(K)
	1             - dTBdCHI_K*AV(I)
	1             + dHSDCHI(I)*CVM1(I)
C
	 T1=T1*(DTAU(J)+DTAU(I))/TCHI(I)
	 VK(I,I,1)=VK(I,I,1)
	1             - dTAdCHI_I*AV(J)
	1             - dTCdCHI_I*AV(K)
	1             - dTBdCHI_I*AV(I)
	1             + dUdCHI*AVM1(I)
	1             - dHSDCHI(J)*CVM1(J)
	1             + dHSDCHI(I)*CVM1(I)
	1             - T1*SOURCE(I)
C
	 VK(I,I,2)=T1
C
	END DO
C
C Now do the boundary conditions.
C
	T1=  ( MIDF(1)*Q(1)*AV(1)
	1      - MIDF(2)*Q(2)*AV(2) )/DTAU(1)/DTAU(1)
	DO L=1,ND
	  VK(1,L,1)=VK(1,L,1)+T1*dTAUdCHI(1,L)
	END DO
	VK(1,1,1)=VK(1,1,1) +
	1          ( PSI(1)*AV(1)- PSIPREV(1)*AVM1(1) )/TCHI(1)
C
	IF(DIF)THEN
	  T1= ( MIDF(ND)*AV(ND) - MIDF(ND-1)*Q(ND-1)*AV(ND-1) )
	1           / DTAU(ND-1)/DTAU(ND-1)
	  DO L=1,ND
	    VK(ND,L,1)=VK(ND,L,1)+T1*dTAUdCHI(ND-1,L)
	  END DO
	  VK(ND,ND,1)=VK(ND,ND,1)-
	1               DBB*R(ND)*R(ND)/3.0D0/TCHI(ND)/TCHI(ND)
	ELSE
	  T1= ( MIDF(ND)*AV(ND) - MIDF(ND-1)*Q(ND-1)*AV(ND-1) )
	1           / DTAU(ND-1)/DTAU(ND-1)
	  DO L=1,ND
	    VK(ND,L,1)=VK(ND,L,1)+T1*dTAUdCHI(ND-1,L)
	  END DO
	END IF
C
C 
C
C Convert from total line opacity, and emissivity to line and continuum
C contibutions.
C
	IF(NM .EQ. 4)THEN
	  DO J=1,ND
	    DO I=1,ND
	      VK(I,J,3)=VK(I,J,1)
	      VK(I,J,4)=VK(I,J,2)
	    END DO			!I
	   END DO			!J
	END IF
C
C Multiply line opacity and emissivity variation by line profile.
C
	DO J=1,ND
	  DO I=1,ND
	    VK(I,J,1)=VK(I,J,1)*DEPTHPRO(J)
	    VK(I,J,2)=VK(I,J,2)*DEPTHPRO(J)
	  END DO
	END DO
C
C 
C
C Now we compute the varaition of the equation which updates the
C flux variation. We don't correct for the line profile, since
C we would the require two matrices. Note that HU(I), HL(I) and
C HS(I) depend directly on CHI(I) and CHI(I+1).
C
	DO I=1,ND-1
	  T1=dHUdTAU(I)*AV(I+1)-dHLdTAU(I)*AV(I)
	  DO L=1,ND
	    RHS_dHdCHI(I,L)=RHS_dHdCHI(I,L)+T1*dTAUdCHI(I,L)
	  END DO
	  T1=dHUdCHI(I)*AV(I+1) - dHLdCHI(I)*AV(I)
	1                       + dHSdCHI(I)*CVM1(I)
	  RHS_dHdCHI(I,I)=RHS_dHdCHI(I,I) + T1
	  RHS_dHdCHI(I,I+1)=RHS_dHdCHI(I,I+1) + T1
	END DO
C
	RETURN
	END
