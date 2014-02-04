!
! This routine solves for the mean intensity as a function of depth using the
! Feautrier Technique for a spherical gray atmosphere. A diffusion approaximation
! is used for the lower boundary condition.
!
! This routine includes the zerth order correction to the transfer equations because
! of the velocity. It arises from the V dI/dv term in the transfer equation, which
! is the only v term usually included for modelling stellar winds. Other terms, which 
! were negelected, are also important to first order in V in the grey transfer equation.
!
	SUBROUTINE JGREY_WITH_FVT(RJ,RSQ_HFLUX,CHI,R,VEL,SIGMA,
	1                  P,JQW,HQW,KQW,NQW,
	1                  LUMINOSITY,METHOD,DIFF_APPROX,IC,
	1                  ACCURACY,ND,NC,NP)
	IMPLICIT NONE
!
!   Altered 06-Nov-2007: Aditional and improved diagnostics inserted.
! Finalized 14-Jan-2003: Gives ame answer as JGREY in limit of VEL=0
!
	INTEGER NC
	INTEGER ND
	INTEGER NP
!
	REAL*8 RJ(ND)			!Mean intensity (computed and returned)
	REAL*8 RSQ_HFLUX(ND)            !r^2 . Flux
	REAL*8 CHI(ND)			!Opacity
	REAL*8 R(ND)			!Radius grid (in units of 10^10 cm)
	REAL*8 VEL(ND)			!Velocity (in km/s)
	REAL*8 SIGMA(ND)		!dlnv/dlnr -1
!
	REAL*8 P(NP)			!Impact parameters
	REAL*8 JQW(ND,NP)		!Quadrature weight for J (on grid)
	REAL*8 KQW(ND,NP) 		!Quadrature weight for K (on grid)
	REAL*8 HQW(ND-1,NP)		!Quadrature weight for H (at midpoints)
	REAL*8 NQW(ND-1,NP)		!Quadrature weight for N (at midpoints)
!
	REAL*8 LUMINOSITY               !Luminosity at inner bounary in Lsun.
	REAL*8 IC
	REAL*8 ACCURACY			!Convergence accuracy for computing f.
	CHARACTER*6 METHOD
	LOGICAL DIFF_APPROX		!Use a diffusion approximation (as opposed to a Schuster core)
!
! Local vectors & arrays
!
	REAL*8 TA(ND),TB(ND),TC(ND)
	REAL*8 HU(ND),HL(ND)
	REAL*8 XM(ND)
	REAL*8 Z(ND)
	REAL*8 CHI_MOD(ND)
	REAL*8 dCHIdR(ND)
	REAL*8 dCHI_MODdR(ND)
	REAL*8 Q(ND)
	REAL*8 F(ND)
	REAL*8 JFAC(ND)
	REAL*8 N_ON_H(ND)
	REAL*8 DTAU(ND)
	REAL*8 BETA(ND)
	REAL*8 VU(ND)
	REAL*8 CV(ND)
	REAL*8 AVE_DTAU(ND)
!
! FS indcates the following quanties (J, H, K & N) have been computed using the
! formal soulution.
!
	REAL*8 FS_RSQJ(ND)
	REAL*8 FS_RSQH(ND)
	REAL*8 FS_RSQK(ND)
	REAL*8 FS_RSQN(ND)
!
	REAL*8 HBC		!Eddington factor for H at outer boundary.
	REAL*8 IN_HBC		!Eddington factor for H at inner boundary.
	REAL*8 HMOD
	REAL*8 IBOUND
	REAL*8 DBB
	REAL*8 DBC
	REAL*8 C_KMS
!
	INTEGER, PARAMETER :: IONE=1
	REAL*8 PI
	REAL*8 T1,T2
	REAL*8 E1,E2,E3
	INTEGER I,NI,LS
	INTEGER LU_DIAG
	REAL*8 SPEED_OF_LIGHT
	EXTERNAL SPEED_OF_LIGHT
	LOGICAL VERBOSE
!
! Set initial values.
!
	VERBOSE=.FALSE.
	LU_DIAG=7
	C_KMS=1.0D-05*SPEED_OF_LIGHT()
	PI=ACOS(-1.0D0)
	DO I=1,ND
	  F(I)=0.33333D0
	  N_ON_H(I)=0.5D0
	  BETA(I)=VEL(I)/C_KMS
	END DO
	HBC=1.0D0
	IN_HBC=0.2D0
	IBOUND=0.0D0
!
! Loop to converge Eddington factors.
!
1000	CONTINUE
!
! Form the sphericity factor Q from F
!
	CALL QFROMF(F,Q,R,TA,TB,ND)	!TA,TB work vectors
!
	DO I=1,ND
	  CHI_MOD(I)=CHI(I)+BETA(I)/R(I)*(1.0D0+SIGMA(I)*N_ON_H(I))
	END DO
!
! Form "SPHERICAL" optical depth scale.
!
	DO I=1,ND
	  TA(I)=Q(I)*CHI_MOD(I)
	END DO
	CALL DERIVCHI(dCHIdR,TA,R,ND,METHOD)
	CALL NORDTAU(DTAU,TA,R,R,dCHIdR,ND)
!
! Compute the optical depth step on the nodes.
!
        DO I=2,ND-1
	  AVE_DTAU(I)=0.5D0*(DTAU(I)+DTAU(I-1))
        END DO
!
! HUL,HL are used to compute r^2.h
!
	DO I=1,ND-1
	  HU(I)=F(I+1)*Q(I+1)/DTAU(I)
	  HL(I)=F(I)*Q(I)/DTAU(I)
	  JFAC(I)=BETA(I)*(1.0D0+SIGMA(I)*F(I))/R(I)/Q(I)/CHI_MOD(I)
	END DO
!
! Compute the TRIDIAGONAL operators, and the RHS source vector. These
! vetors are corrupted in the solution.
!
        DO I=2,ND-1
          TA(I)=HL(I-1)
          TC(I)=HU(I)
          TB(I)=-AVE_DTAU(I)*JFAC(I) - HL(I) - HU(I-1)
          XM(I)=0.0D0
        END DO
!
! Evaluate TA,TB,TC for boudary conditions
!
! Outer boundary
!
        TC(1)=-F(2)*Q(2)/DTAU(1)
        TB(1)= F(1)*Q(1)/DTAU(1) + HBC
        XM(1)=0.0D0
        TA(1)=0.0D0
!
! NB: HMOD is actually R^2 . H
!
	HMOD=3.826D+13*LUMINOSITY/16.0D0/PI/PI
        TA(ND)=-Q(ND-1)*F(ND-1)/DTAU(ND-1)
        IF(DIFF_APPROX)THEN
          TB(ND)=F(ND)/DTAU(ND-1)
          XM(ND)=HMOD
        ELSE
          TB(ND)=F(ND)/DTAU(ND-1)+IN_HBC
          XM(ND)=R(ND)*R(ND)*IC*(0.25D0+0.5D0*IN_HBC)
        END IF
        TC(ND)=0.0D0
!
	CALL THOMAS(TA,TB,TC,XM,ND,1)
	RJ(1:ND)=XM(1:ND)/R(1:ND)/R(1:ND)
!
	DO I=1,ND-1
	  RSQ_HFLUX(I)=HU(I)*XM(I+1)-HL(I)*XM(I)
	END DO
!
	FS_RSQJ(:)=0.0D0
	FS_RSQH(:)=0.0D0
	FS_RSQK(:)=0.0D0
	FS_RSQN(:)=0.0D0
	HBC=0.0D0
!
! DBB =3L/16(piR)**2 and is used for the lower boundary diffusion approximation.
!
	T1=3.826D+13*LUMINOSITY/16.0D0/PI/PI
	DBB=3.0D0*T1/R(ND)/R(ND)
!
	CALL DERIVCHI(dCHIdr,CHI,R,ND,METHOD)
!
! ENTER LOOP FOR EACH IMPACT PARAMETER P
!
	DO 2000 LS=1,NP
	  NI=ND-(LS-NC-1)
	  IF(LS .LE. NC)THEN
	    NI=ND
	    DBC=DBB*DSQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)
	  END IF
!
! Compute Z for this imapct parameter
!
	  IF(NI .GT. 1)THEN
	    DO I=1,NI
	      Z(I)=DSQRT( (R(I)-P(LS))*(R(I)+P(LS)) )
	      T1=Z(I)/R(I)
	      TA(I)=BETA(I)*(1.0D0+SIGMA(I)*T1*T1)/R(I)
	      CHI_MOD(I)=CHI(I)+TA(I)
	    END DO
	    IF(NI .EQ. 2)THEN
	      TB(1)=0.0D0         !Cance so values unimportant
	      TB(2)=TB(1)
	    ELSE
	      CALL DERIVCHI(TB,TA,R,NI,METHOD)
	    END IF
	    DO I=1,NI
	      dCHI_MODdR(I)=dCHIdR(I)+TB(I)
	    END DO
	  END IF
!
! Compute Z for this imapct parameter
!
	  IF(NI .GT. 2)THEN
	    DO I=1,NI-1
	      DTAU(I)=0.5D0*(Z(I)-Z(I+1))*(CHI_MOD(I)+CHI_MOD(I+1)+(Z(I)-Z(I+1))
	1       *(dCHI_MODdR(I+1)*Z(I+1)/R(I+1)-dCHI_MODdR(I)*Z(I)/R(I))/6.0D0)
	    END DO
	  END IF
!
! Compute TA (tridiagonal matric) and XM vector. Since it is a grey
! atmosphere, the source function is simply RJ*CHI(I)/CHI_MOD(I).
!
	  IF(NI .GT. 2)THEN
!
	    DO I=1,NI-1
	      VU(I)=1.0D0/DTAU(I)
	    END DO
! 
	    XM(1)=-IBOUND
	    TA(1)=0.0D0
	    TC(1)=1.0D0/DTAU(1)
	    TB(1)=-1.0D0-TC(1)
	    DO I=2,NI-1
	      TA(I)=VU(I-1)
	      TC(I)=VU(I)
	      TB(I)=-0.5D0*(DTAU(I-1)+DTAU(I))-TA(I)-TC(I)
	      XM(I)=-RJ(I)*CHI(I)*(DTAU(I-1)+DTAU(I))*0.5D0/CHI_MOD(I)
	    END DO
!
	    IF(LS .LE. NC .AND. DIFF_APPROX)THEN
	      TB(NI)=VU(NI-1)
	      TA(NI)=-VU(NI-1)
	      XM(NI)=DBC
	    ELSE IF(LS .GT. NC)THEN
	      TA(NI)=-TC(NI-1)
	      TB(NI)=-TA(NI)+DTAU(NI-1)/2.0D0
	      XM(NI)=0.5D0*DTAU(NI-1)*RJ(NI)*CHI(NI)/CHI_MOD(NI)
	    ELSE
	      TA(NI)=-VU(NI-1)
	      TB(NI)=1.0D0+VU(NI-1)
	      XM(NI)=IC
	    END IF
	    TC(NI)=0.0D0
!
! Solve the tridiagonal system of equations.
!
	    CALL THOMAS(TA,TB,TC,XM,NI,IONE)
C
	  ELSE IF(NI .EQ. 1)THEN
	    XM(1)=0.0D0
	  ELSE IF(NI .EQ. 2)THEN
	    Z(1)=SQRT(R(1)*R(1)-P(LS)*P(LS))
	    DTAU(1)=0.5D0*Z(1)*(CHI_MOD(1)+CHI_MOD(2))		!Z(2)=0.0
	    E1=DEXP(-DTAU(1))
	    E2=1.0D0-(1.0D0-E1)/DTAU(1)
	    E3=(1.0D0-E1)/DTAU(1)-E1
	    IF(DTAU(1) .LT. 1.0D-03)THEN
	      E2=DTAU(1)*0.5D0+DTAU(1)*DTAU(1)/6.0D0
	      E3=DTAU(1)*0.5D0-DTAU(1)*DTAU(1)/3.0D0
	    END IF
C
	    XM(2)=TA(2)*E2+TA(1)*E3
            XM(1)=0.5D0*(XM(2)*E1+TA(1)*E2+TA(2)*E3)
	  END IF
!
! Update J, H, K, & N for this angle.
!
	  DO I=1,NI
	    FS_RSQJ(I)=FS_RSQJ(I)+JQW(I,LS)*XM(I)
	    FS_RSQK(I)=FS_RSQK(I)+KQW(I,LS)*XM(I)
	  END DO
!
! CV is the flux (v in usual Mihalas notation).
!
	  DO I=1,NI-1
	    CV(I)=VU(I)*(XM(I+1)-XM(I))
	    FS_RSQH(I)=FS_RSQH(I)+HQW(I,LS)*CV(I)
	    FS_RSQN(I)=FS_RSQN(I)+NQW(I,LS)*CV(I)
	  END DO
!
	  HBC=HBC+JQW(1,LS)*(XM(1)-IBOUND)*Z(1)/R(1)
!  
2000	CONTINUE
!
! Compute the new Feautrier factors. These are stored in FS_RSQK so as not
! to destroy the old factors.
!
	DO I=1,ND
	  FS_RSQK(I)=FS_RSQK(I)/FS_RSQJ(I)
	END DO
	DO I=1,ND-1
	  N_ON_H(I)=FS_RSQN(I)/FS_RSQH(I)
	END DO
!
! Compute the factor for the outer boundary condition.
!
	HBC=HBC/FS_RSQJ(1)
!
	T1=MAXVAL( ABS(F(1:ND)-FS_RSQK(1:ND)) )
	WRITE(6,*)'Current accuracy is', T1
	IF(T1 .GT. ACCURACY)THEN
	  F(1:ND)=FS_RSQK(1:ND)
	  GOTO 1000
	END IF
!
! 
!
! The following allows checking whether difference equatins were evaluated successfully.
! It evaluations should be identical to that above --- done again becuase terms were corrupted.
!
	IF(VERBOSE)THEN
	  DO I=1,ND
	    CHI_MOD(I)=CHI(I)+BETA(I)/R(I)*(1.0D0+SIGMA(I)*N_ON_H(I))
	    TA(I)=Q(I)*CHI_MOD(I)
	  END DO
	  CALL DERIVCHI(dCHIdR,TA,R,ND,METHOD)
	  CALL NORDTAU(DTAU,TA,R,R,dCHIdR,ND)
!
! Compute the TRIDIAGONAL operators, and the RHS source vector. These
! vetors are corrupted in the solution.
!
          DO I=2,ND-1
            TA(I)=HL(I-1)
            TC(I)=HU(I)
            TB(I)=-AVE_DTAU(I)*JFAC(I) - HL(I) - HU(I-1)
            XM(I)=0.0D0
          END DO
!
! Evaluate TA,TB,TC for boudary conditions
!
! Outer boundary
!
          TC(1)=-F(2)*Q(2)/DTAU(1)
          TB(1)= F(1)*Q(1)/DTAU(1) + HBC
          XM(1)=0.0D0
          TA(1)=0.0D0
!
! NB: HMOD is actually R^2 . H
!
	  HMOD=3.826D+13*LUMINOSITY/16.0D0/PI/PI
          TA(ND)=-Q(ND-1)*F(ND-1)/DTAU(ND-1)
          IF(DIFF_APPROX)THEN
            TB(ND)=F(ND)/DTAU(ND-1)
            XM(ND)=HMOD
          ELSE
            TB(ND)=F(ND)/DTAU(ND-1)+IN_HBC
            XM(ND)=R(ND)*R(ND)*IC*(0.25D0+0.5D0*IN_HBC)
          END IF
          TC(ND)=0.0D0
!
          OPEN(UNIT=LU_DIAG,STATUS='UNKNOWN',FILE='NEW_GREY_CHK')
	    Z(1:ND)=RJ(1:ND)*R(1:ND)*R(1:ND)
	    WRITE(LU_DIAG,'(I5,ES12.5,3ES12.3,2ES14.5)')1,R(1),0.0D0,TB(1)*Z(1),TC(1)*Z(2),
	1                     TB(1)*Z(1)+TC(1)*Z(2),XM(1)
	    DO I=2,ND-1
	      WRITE(LU_DIAG,'(I5,ES12.5,3ES12.3,2ES14.5)')I,R(I),TA(I)*Z(I-1),TB(I)*Z(I),TC(I)*Z(I+1),
	1                     TA(I)*Z(I-1)+TB(I)*Z(I)+TC(I)*Z(I+1),XM(I)
	    END DO
	    WRITE(LU_DIAG,'(I5,ES12.5,3ES12.3,2ES14.5)')ND,R(ND),TA(ND)*Z(ND-1),TB(ND)*Z(ND),0.0D0,
	1                     TA(ND)*Z(ND-1)+TB(ND)*Z(ND),XM(ND)
	  CLOSE(UNIT=LU_DIAG)
!
	END IF		!VERBOSE
!
! Check that "Radiative Equilibrium is being achieved. We first regrid H from the mid-points onto
! the nodes.
!
! Allowing for the negeleted terms, we should have
!
!       H(d=1) = H(d) - Integral[r,rmax] r. V/c . J
!
! Because H(d) >> H(d=1) at depth, may not be accurately satisfied at depth
! (Integration rule has errors etc).
!
! Alternatively:
!       H(d) = H(ND) + Integral[rmin,r] r. V/c . J
!
! This later formulation does not easily allow the conservation equation to be checked
! in the outer regions --- the most important region.
!
        T1=HBC*RJ(1)
        T2=HMOD/R(ND)/R(ND)
        CALL REGRID_H(TB,R,RSQ_HFLUX,T1,T2,ND,TC)
	DO I=1,ND
	  TA(I)=R(I)*VEL(I)*RJ(I)*(1.0D0+SIGMA(I)*F(I))/C_KMS
	END DO
	TC(1:ND)=TA(1:ND)
	CALL LUM_FROM_ETA(TA,R,ND)
!
! Now calculate the integrated correction.
!
	Q(1)=0.0D0
	DO I=2,ND
	  Q(I)=Q(I-1)+TA(I-1)
	END DO
        OPEN(UNIT=LU_DIAG,STATUS='UNKNOWN',FILE='GREY_CHK')
          WRITE(LU_DIAG,'(A)')' '
          WRITE(LU_DIAG,'(A)')' Check on GREY solution.'
          WRITE(LU_DIAG,'(A)')' Nonrelativistic case but with velocity terms.'
          WRITE(LU_DIAG,'(A)')' Routine used is JGREY_WITH_FVT'
          WRITE(LU_DIAG,'(A)')' '
	  WRITE(LU_DIAG,*)'HBC=',HBC
          WRITE(LU_DIAG,'(A)')' '
!
          WRITE(LU_DIAG,'(4X,A,2X,10(3X,A))')'I','        R','  V(km/s)','      Chi','    rsq.J','    rsq.H',
	1                                                    '    rsq.H','    Jterm','     dInt','      Int',
	1                                                    ' Cons. Flux'
	  RSQ_HFLUX(ND)=0.0D0
	  DO I=1,ND
	    T1=RJ(I)*R(I)*R(I)
	    WRITE(LU_DIAG,'(I5,2X,ES12.5,8ES12.3,ES14.5)')I,R(I),VEL(I),CHI(I),T1,RSQ_HFLUX(I),
	1                                TB(I),TC(I),TA(I),Q(I),TB(I)-Q(I)
	  END DO
!
! Compare solutions of the zero moment equation, allowing for the different averaging of CHI.
!
	  WRITE(LU_DIAG,'(A)')' '
	  WRITE(LU_DIAG,'(A)')' Comparison with differencing in R'
	  WRITE(LU_DIAG,'(A)')' In ideal world, last 2 terms should be identical'
	  WRITE(LU_DIAG,'(A)')' In practice, not identical because of chi averaging'
	  WRITE(LU_DIAG,'(A)')' '
	  WRITE(LU_DIAG,'(4X,A,4(7X,A),2(6X,A))')'I','    R','    V','  Chi','rsq.J',' drsqHdr','r.Beta.J'
	  DO I=2,ND-1
	    T1=RJ(I)*R(I)*R(I)
	    T2=2.0D0*(RSQ_HFLUX(I-1)-RSQ_HFLUX(I))/(R(I-1)-R(I+1))
	    WRITE(LU_DIAG,'(I5,2X,ES12.5,3ES12.3,2ES14.5)')I,R(I),VEL(I),CHI(I),T1,T2,BETA(I)*RJ(I)*R(I)
	  END DO
!
	  WRITE(LU_DIAG,'(A)')' '
	  WRITE(LU_DIAG,'(A)')' Comparison with differencing in TAU'
	  WRITE(LU_DIAG,'(A)')' Last 2 terms should be identical'
	  WRITE(LU_DIAG,'(A)')' '
	  WRITE(LU_DIAG,'(4X,A,4(7X,A),2(6X,A))')'I','    R','    V','  Chi','rsq.J',' drsqHdr','rsq.Jfac'
	  DO I=2,ND-1
	    T1=RJ(I)*R(I)*R(I)
	    T2=(RSQ_HFLUX(I-1)-RSQ_HFLUX(I))/AVE_DTAU(I)
	    WRITE(LU_DIAG,'(I5,2X,ES12.5,3ES12.3,2ES14.5)')I,R(I),VEL(I),CHI(I),T1,T2,JFAC(I)*T1
	  END DO
!
! We use VU for RJ in the call to JGREY.
!
	WRITE(LU_DIAG,'(A)')' '
	WRITE(LU_DIAG,'(A)')' '
	WRITE(LU_DIAG,'(A)')' Comparison with static solution'
	WRITE(LU_DIAG,'(A)')' NB: This only provides a check when V(ND) is effectively zero.'
	WRITE(LU_DIAG,'(A)')' '
	WRITE(LU_DIAG,'(4X,A,6(2X,A))')'I','         R','   V(km/s)','    J(cmf)','   J(stat)',
	1                                  ' rsqH(cmf)','rsqH(stat)'
	WRITE(LU_DIAG,'(A)')' '
	DO I=1,10
	  CALL JGREY(TA,TB,TC,XM,DTAU,R,Z,P,VU,FS_RSQJ,FS_RSQK,Q,F,
	1         CHI,dCHIdr,JQW,KQW,LUMINOSITY,HBC,T2,NC,ND,NP,METHOD)
	  HBC=T2
	END DO
	T1=3.826D+13*LUMINOSITY/16.0D0/PI/PI		!E-20 as R in units of 10^10 cm
	DO I=1,ND
	  WRITE(LU_DIAG,'(X,I4,ES12.5,5ES12.3)')I,R(I),VEL(I),RJ(I),VU(I),RSQ_HFLUX(I),T1
	END DO
!
	CLOSE(UNIT=LU_DIAG)
!
	RETURN
	END
