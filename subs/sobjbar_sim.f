C
C Subroutine to compute JBAR and ZNET using the escape probability
C approximation. Six vectors are returned with new values.
C
C Input:
C
C    CHIL         = Total line opacity
C    ETAL         = Total line emissivity
C    CHIL_MAT(,J) = Line opacity for line J
C    ETAL_MAT(,J) = Line emissivity for line J
C    BBCOR(,J)    = Correction to emissivity to allow for frequncy difference.
C
C Output: Note that VB_2 and VC_2 are independent of the line J
C
C    JBAR  = Line mean intensity
C    BETA  = Sobolev Escape Probability.
C    VB_2  : dZnet(,J)/dCHIL(,K)=STOT*VB_2/SRCE(,J)
C    VC_2  : dZnet(,J)/dETAL(,K)=BB_COR(,J)*STOT*VC_2/SRCE(,J)
C
C Matrices (1 vector for each line)
C
C    ZNET(,J)    = Net rate line for line J
C    VB_SIM(,J)  = dZNET(,J)/dCHIL(,J)
C    VC_SIM(,J)  = dZNET(,J)/dETAL(,J)
C    BETAC(,J)   = Sobolev continuum term appearing in net rate.
C
	SUBROUTINE SOBJBAR_SIM(SOURCE,CHI,CHIL,ETAL,V,SIGMA,R,P,AQW,
	1                  JBAR,BETA,
	1                  ZNET,VB,VC,VB_2,VC_2,BETAC,
	1                  CHIL_MAT,ETAL_MAT,BBCOR,
	1                  FL,DIF,DBB,IC,THICK,
	1                  NLF,NC,NP,ND,NSIM,METHOD)
	IMPLICIT NONE
C
C Altered 26-Oct-2014 : We now write ABS(dmuV/dz) so as to handle non-monotonic velocity laws.
C Altered 22-Aug-1997 : Adjusted to handle ND=NC+NP (i.e. NI=1 or 2 is
C                         now handled --- based on FQCOMP_IBC).
C Altered 28-Oct-1996 : Bug Fix: COS converted back to ACOS in TOR expression.
C Altered 26-May-1996 : Scratch blocked removed.
C                       Dynamic allocation for scratch vectors.
C                       IONE pased in call to THOMAS.
C                       Genrical calls to EXP, COS
C
C Altered 15-Jan-1992 VB_2, VC_2 inserted. Bug fixes.
C Created 13-Nov-1992 Allows for the simultaenous treatment of several
C                     lines. Based on SOBJBAR.
C
	INTEGER NLF,NC,NP,ND,NSIM
	REAL*8 SOURCE(ND),CHI(ND),CHIL(ND),ETAL(ND)
	REAL*8 V(ND),SIGMA(ND),R(ND),P(NP),AQW(ND,NP)
C
	REAL*8 JBAR(ND)
	REAL*8 BETA(ND)
C
	REAL*8 ZNET(ND,NSIM)
	REAL*8 VB(ND,NSIM)
	REAL*8 VC(ND,NSIM)
	REAL*8 VB_2(ND)
	REAL*8 VC_2(ND)
	REAL*8 BETAC(ND,NSIM)
C
	REAL*8 CHIL_MAT(ND,NSIM)
	REAL*8 ETAL_MAT(ND,NSIM)
	REAL*8 BBCOR(ND,NSIM)
C
	REAL*8 DBB,IC,FL
	CHARACTER*(*) METHOD
	LOGICAL DIF,THICK
C
	REAL*8 TA(ND),TB(ND),TC(ND),AV(ND),DTAU(ND),Z(ND)
	REAL*8 GAM(ND),GAMH(ND),dCHIdR(ND)
	REAL*8 EX_VEC(ND),dEX_VEC(ND)
	REAL*8 ZNET_BL(ND),BETAC_BL(ND)
	REAL*8 VB_BL(ND),VC_BL(ND)
	REAL*8 ST_ON_SL
C
	REAL*8 EXPONX,d_EXPONX_dX
C
C Local variables.
C
	INTEGER, PARAMETER :: IONE=1
	INTEGER I,J,LS,NI
	REAL*8 T1,T2,DBC,TOR
	REAL*8 IBOUND,E1,E2,E3  		!To handle NI=2
C
C Zero arrays which are incremented as we integrate over angle.
C Evaluate the SOBOLEV optical depth (GAMH) without angle factor.
C
	DO I=1,ND
	  JBAR(I)=0.0D0
	  BETA(I)=0.0D0
	  GAMH(I)=CHIL(I)*3.0D-10*R(I)/V(I)/FL    	!C/dex(15)/dex(5)
	  ZNET_BL(I)=0.0D0
	  VB_BL(I)=0.0D0
	  VC_BL(I)=0.0D0
	  BETAC_BL(I)=0.0D0
	END DO
C
	CALL DERIVCHI(dCHIdR,CHI,R,ND,METHOD)
C
C Enter loop to perform integration along each ray.
C
	DO LS=1,NP
	  NI=ND-(LS-NC-1)
	  IF(LS .LE. NC)NI=ND
C
C Initialize AV vector. GAM is the ANGLE DEPENDENT Sobolev optical depth.
C
	  AV(1:NI)=0.0D0
	  CALL ZALONGP(R,Z,P(LS),NI)
	  DO I=1,NI
	    GAM(I)=GAMH(I)/ABS(1.0D0+Z(I)*Z(I)/R(I)/R(I)*SIGMA(I))
	  END DO
C
C SOURCE(1) is the boundary continuum source function.
C
	  IF(THICK)THEN
	    IF(P(LS) .GT. 0)THEN
	      TOR=CHI(1)*R(1)*R(1)*(1.570796-ACOS(P(LS)/R(1)))/P(LS)
	    ELSE
 	      TOR=CHI(1)*R(1)
	    END IF
	    IBOUND=SOURCE(1)*(1.0D0-EXP(-TOR))
	  END IF
C
C For NI=1 and 2 AV reprents the intensity variable U (the average
C intensity on the ray). For NI>2, AV only represnts u after the call
C to THOMAS.
C
	  IF(NI .EQ. 1)THEN
	    AV(1)=IBOUND
	  ELSE IF(NI .EQ. 2)THEN
	    CALL NORDTAU(DTAU,CHI,Z,R,dCHIdR,NI)
	    E1=EXP(-DTAU(1))
	    E2=1.0D0-(1.0D0-E1)/DTAU(1)
	    E3=(1.0D0-E1)/DTAU(1)-E1
	    IF(DTAU(1) .LT. 1.0D-03)THEN
	      E2=DTAU(1)*0.5+DTAU(1)*DTAU(1)/6.0D0
	      E3=DTAU(1)*0.5-DTAU(1)*DTAU(1)/3.0D0
	    END IF
	    AV(2)=IBOUND*E1+SOURCE(2)*E2+SOURCE(1)*E3
            AV(1)=0.5*(IBOUND+AV(2)*E1+SOURCE(1)*E2+SOURCE(2)*E3)
	  ELSE
	    IF(DIF .AND. LS .LE. NC)THEN
	      DBC=DBB*SQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)/CHI(ND)
	    END IF
	    CALL NORDTAU(DTAU,CHI,Z,R,dCHIdR,NI)
	    CALL TCOMPD(TA,TB,TC,DTAU,DIF,LS,NC,ND,NI)
	    CALL XVECD(DTAU,SOURCE,AV,DIF,DBC,IC,LS,NC,ND,NI)
	    AV(1)=-IBOUND
C
C Solve for the radiation field along ray for this frequency.
C
	    CALL THOMAS(TA,TB,TC,AV,NI,IONE)
	  END IF
C
C The function EXPONX is given by (1.0-EXP(-X))/X.
C The function d_EXPONX_dX is given by d[ (1.0-EXP(-X))/X ]/dX. These
C functions are called to allow for cancellation when X is small. Note
C GAM(I)/CHIL(I) = d_GAM(I)/d_CHIL(I) .
C
	  DO I=1,NI
            EX_VEC(I)=EXPONX(GAM(I))*AQW(I,LS)
            dEX_VEC(I)=d_EXPONX_dX(GAM(I))*AQW(I,LS)
	    BETA(I)=BETA(I)+EX_VEC(I)
	  END DO
C
	  DO I=1,NI
	    T2=AV(I)*CHIL(I)/ETAL(I)
	    T1=1.0D0-T2
C
	    ZNET_BL(I)=ZNET_BL(I)+T1*EX_VEC(I)
	    BETAC_BL(I)=BETAC_BL(I)+EX_VEC(I)*AV(I)
	    VB_BL(I)=VB_BL(I)+( T1*dEX_VEC(I)*GAM(I)/CHIL(I)-
	1                 AV(I)*EX_VEC(I)/ETAL(I))
	    VC_BL(I)=VC_BL(I)+T2*EX_VEC(I)/ETAL(I)
C
	  END DO		!End I loop
C
	END DO			!End LS Loop
C
C Now compute the net rates, and variations of the net rates. At this
C stage it is not clear whether this should be done as below, or done
C in loop above to improve accuracy.
C
	DO J=1,NSIM
	  DO I=1,ND
            ST_ON_SL=(ETAL(I)/ETAL_MAT(I,J))*
	1               (CHIL_MAT(I,J)/CHIL(I))/BBCOR(I,J)
	    ZNET(I,J)=ST_ON_SL*ZNET_BL(I)+(1.0D0-ST_ON_SL)
	    BETAC(I,J)=BETAC_BL(I)*CHIL_MAT(I,J)/ETAL_MAT(I,J)/BBCOR(I,J)
            VB(I,J)=ST_ON_SL*( VB_BL(I)+(ZNET_BL(I)-1.0D0)*
	1                      (1.0D0/CHIL_MAT(I,J)-1.0D0/CHIL(I)) )
            VC(I,J)=ST_ON_SL*( VC_BL(I)+(ZNET_BL(I)-1.0D0)*
	1                      (BBCOR(I,J)/ETAL(I)-1.0D0/ETAL_MAT(I,J)) )
	  END DO
	END DO
C
	DO I=1,ND
          VB_2(I)=VB_BL(I)-(ZNET_BL(I)-1.0D0)/CHIL(I)
          VC_2(I)=VC_BL(I)+(ZNET_BL(I)-1.0D0)/ETAL(I)
	END DO
C
C Compute the mean line intensity. This is not very accurate, and
C fails if CHIL(I,1) .EQ. 0
C
	DO I=1,ND
	  JBAR(I)=(1.0D0-ZNET_BL(I))*ETAL(I)/CHIL(I)
	END DO
C
	RETURN
	END
