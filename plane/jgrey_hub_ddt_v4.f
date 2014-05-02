!
! This routine solves for the mean intensity as a function of depth using the
! Feautrier Technique for a spherical gray atmosphere. A diffusion approximation
! is used for the lower boundary condition. Radiative equilibrium is assumed.
!
! The time dependent radiative transfer equation is solved including all terms 
! to first order in v/c, and assumes a Hubble flow. Thus sigma=dlnv/dlnr-1 =0. 
! A Langrangian description is used to handle the time dependence. The frequency 
! integrated moments J and H (and the associated boundary conditions) at the previous
! time step need to be available, and are read in from the file JH_AT_PREV_TIME.
! If desired, the full D/Dt term can be neglected, in which case the code should give 
! identical answers to JGREY_WITH_FVT.
!
! Partially based on the routine ../subs/jgrey_with_fvt.f
!
	SUBROUTINE JGREY_HUB_DDT_V4(RJ,RSQ_HFLUX,CHI,CHI_PLANCK,R,VEL,SIGMA,POPS,
	1              P,JQW,HQW,KQW,LUMINOSITY,METHOD,DIFF_APPROX,IC,
	1              ACCURACY,DO_TIME_VAR,TIME_SEQ_NO,COMPUTED,ND,NC,NP,NT)
	IMPLICIT NONE
!
! Altered 14-Nov-2014: Changed to V4 and installed COMPUTED in call.
! Altered 04-Nov-2012; WORK was not being set before evaluating f.
! Altered 28-Nov-2011: REDUCTION_FATOR reduced to 10, and BAD_J_COUNTER now only updated
!                          per iteration, not per bad depth. WORK is now set to zero
!                          before jumping to compute "f" on the first iteraton step. The
!                          advanced computation of "f" was added, since a initial estimate for 
!                          "f" was causing some of the ocnvergnce difficuties.
! Altered 13-Nov-2011: Better handling of error conditions - RJ adjusted before T udated.
!                          Comments added 23-Nov-2011.
! Altered 06-Jun-2010: Assume zero-flux option when not diffusion option.
! Altered 05-Aug-2008: Fixed bug with check SUM written to SN_GREY_CHK.
!                          Replaced T1+T2+T3-E3 with T1+T2+T3-E3 (.ie., -T3 --> +T3)
!                          NB: E3 refers to the work term, T3 to the radictive decay term.
!                          Earlier fixed from July 08 was updated.
!                          Included BAD_J_COUNTER and F_LOOP_COUNTER
! Created 22-July-2006
!
	INTEGER NC
	INTEGER ND
	INTEGER NP
	INTEGER NT
	INTEGER TIME_SEQ_NO
!
	REAL*8 RJ(ND)			!Mean intensity (computed and returned)
	REAL*8 RSQ_HFLUX(ND)            !r^2 . Flux
	REAL*8 R(ND)			!Radius grid (in units of 10^10 cm)
	REAL*8 CHI_PLANCK(ND)
	REAL*8 CHI(ND)			!Opacity
	REAL*8 VEL(ND)			!Velocity (in km/s)
	REAL*8 SIGMA(ND)		!dlnV/dlnr-1
	REAL*8 POPS(NT,ND)
!
	REAL*8 P(NP)			!Impact parameters
	REAL*8 JQW(ND,NP)		!Quadrature weight for J (on grid)
	REAL*8 KQW(ND,NP) 		!Quadrature weight for K (on grid)
	REAL*8 HQW(ND,NP)		!Quadrature weight for H (at midpoints)
!
	REAL*8 LUMINOSITY               !Luminosity at inner boundary in Lsun.
	REAL*8 IC                       !Not used
	REAL*8 ACCURACY			!Convergence accuracy for computing f.
	LOGICAL COMPUTED		!Indicates whether grey solution was successfully computed
	LOGICAL DIFF_APPROX		!Use a diffusion approximation (as opposed to a Schuster core)
	CHARACTER(LEN=6) METHOD
!
! Local vectors & arrays
!
	REAL*8 TA(ND),TB(ND),TC(ND)
	REAL*8 XM(ND)
	REAL*8 HU(ND)
	REAL*8 HL(ND)
	REAL*8 HT(ND)
	REAL*8 JFAC(ND)
	REAL*8 JT(ND)
	REAL*8 WMID(ND)
	REAL*8 Z(ND)
	REAL*8 CHI_MOD(ND)
	REAL*8 dCHIdR(ND)
	REAL*8 dCHI_MODdR(ND)
	REAL*8 Q(ND)
	REAL*8 F(ND)
	REAL*8 DTAU(ND)
	REAL*8 BETA(ND)
	REAL*8 VU(ND)
	REAL*8 CV(ND)
	REAL*8 AVE_DTAU_ON_Q(ND)
	REAL*8 RAD_DECAY_TERM(ND)
	REAL*8 E_RAD_DECAY(ND)
	REAL*8 WORK(ND)
	REAL*8 OLD_T(ND)
	REAL*8 T_FROM_J(ND)
	REAL*8 TSTORE(ND)
!
	REAL*8 RSQ_J_OLDT(ND)
	REAL*8 RSQ_H_OLDT(ND)
!
! FS indicates the following quantities (J, H, K & N) have been computed using the
! formal solution.
!
	REAL*8 FS_RSQJ(ND)
	REAL*8 FS_RSQH(ND)
	REAL*8 FS_RSQK(ND)
!
	REAL*8 H_OUTBC		!Eddington factor for H at outer boundary.
	REAL*8 H_INBC
	REAL*8 H_OUTBC_OLDT	!Eddington factor for H at outer boundary (prev. time step).
	REAL*8 H_INBC_OLDT
	REAL*8 DTAU_INB_SAVE
	REAL*8 IBOUND
	REAL*8 DBB
	REAL*8 DBC
	REAL*8 C_KMS
	REAL*8 RECIP_CDELTAT
	REAL*8 ROLD_ON_R
	REAL*8 DELTA_TIME_SECS
!
	INTEGER, PARAMETER :: IONE=1
	REAL*8 PI
	REAL*8 T1,T2,T3
	REAL*8 E1,E2,E3
	REAL*8 REDUCTION_FACTOR
	INTEGER I,J,NI,LS
	REAL*8 SPEED_OF_LIGHT
	EXTERNAL SPEED_OF_LIGHT
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	INTEGER LUER,ERROR_LU,LU
	INTEGER LUWARN,WARNING_LU
	INTEGER BAD_J_COUNTER
	INTEGER F_LOOP_COUNTER
	EXTERNAL ERROR_LU,WARNING_LU
	LOGICAL DO_TIME_VAR
	LOGICAL BAD_J
	LOGICAL VERBOSE_OUTPUT
!
! Set initial values.
!
	C_KMS=1.0D-05*SPEED_OF_LIGHT()
	PI=ACOS(-1.0D0)
	LUER=ERROR_LU()
	LUWARN=WARNING_LU()
        CALL GET_VERBOSE_INFO(VERBOSE_OUTPUT)
	DO I=1,ND
	  F(I)=0.33333D0
	  BETA(I)=VEL(I)/C_KMS
	END DO
	H_OUTBC=1.0D0
	H_INBC=0.2D0
	IBOUND=0.0D0
	COMPUTED=.TRUE.
!
	CALL GET_RAD_DECAY_ENERGY(E_RAD_DECAY,ND)
	T1=1.0D+10/4.0D0/PI
	E_RAD_DECAY=T1*E_RAD_DECAY
!
! When DO_TIME_VAR=.FALSE., this routine should give the same answers as 
! JGREY_WITH_FVT.
!
! NB: The factor of 10^10 occurs because c. /\t is a length, and R in
!     cmfgen is in units of 10^10 cm. NB: In the differenced equations
!     we always have terms like 1/(c . /\t . chi)
!
!     The factor of 10^{-5} in ROLD_ON_R occurs because V is in km/s, and
!      R is in units of 10^10cm..
!
	IF(DO_TIME_VAR)THEN
	  T1=0.0D0
	  CALL GET_JH_AT_PREV_TIME_STEP(RSQ_J_OLDt,RSQ_H_OLDt,
	1        H_INBC_OLDT,H_OUTBC_OLDT,DELTA_TIME_SECS,
	1        T1,R,VEL,ND,L_TRUE,'GREY')
	  RECIP_CDELTAT=1.0D+10/SPEED_OF_LIGHT()/DELTA_TIME_SECS
	  ROLD_ON_R=1.0D0-1.0D-05*VEL(ND)*DELTA_TIME_SECS/R(ND)
	ELSE
	  RECIP_CDELTAT=0.0D0
	  ROLD_ON_R=1.0D0
	  RSQ_J_OLDt=0.0D0
	  RSQ_H_OLDt=0.0D0
	  H_INBC_OLDT=0.0D0
	  H_OUTBC_OLDT=0.0D0
	  DELTA_TIME_SECS=0.0D0
	END IF
!
! Loop to converge Eddington factors.
!
	T_FROM_J(1:ND)=0.0D0
	OLD_T(1:ND)=0.0D0
	TSTORE(1:ND)=0.0D0
	REDUCTION_FACTOR=10.0D0
        BAD_J_COUNTER=0
        F_LOOP_COUNTER=0
!
	F(1:ND)=0.0D0
	RJ(1:ND)=RSQ_J_OLDt(1:ND)/R(1:ND)/R(1:ND)
!
! We need WORK to evaluate f.
!
	IF(DO_TIME_VAR)THEN
	  CALL DDT_WORK(WORK,POPS,T_FROM_J,OLD_T,TIME_SEQ_NO,ND,NT)
	ELSE
	  OLD_T=0.0D0
	  T_FROM_J=0.0D0
	  WORK=0.0D0
	END IF
!
! Loop to evaluate the Eddington factor f. A poor initial estimate was causing 
! convergence issues in some models.
!
	GOTO 5000
!
1000	CONTINUE
!
	IF(DO_TIME_VAR)THEN
	  CALL DDT_WORK(WORK,POPS,T_FROM_J,OLD_T,TIME_SEQ_NO,ND,NT)
	ELSE
	  OLD_T=0.0D0
	  T_FROM_J=0.0D0
	  WORK=0.0D0
	END IF
	WORK=WORK/REDUCTION_FACTOR
!
! NB: Some of the vectors below could be computed outside of the
!     loop since they do not depend on J. However, most of the
!     time is spent in the impact parameter loop.
!
! Form the sphericity factor Q from F
!
	CALL QFROMF(F,Q,R,TA,TB,ND)	!TA,TB work vectors
!
! Form "SPHERICAL" optical depth scale.
!
	TA(1:ND)=CHI(1:ND)*Q(1:ND)
	CALL DERIVCHI(dCHIdR,TA,R,ND,METHOD)
	CALL NORDTAU(DTAU,TA,R,R,dCHIdR,ND)
	DTAU_INB_SAVE=DTAU(ND-1)
!
! Compute the optical depth step on the nodes.
!
	AVE_DTAU_ON_Q(1)=0.0D0; AVE_DTAU_ON_Q(ND)=0.0D0
        DO I=2,ND-1
          AVE_DTAU_ON_Q(I)=0.5D0*(DTAU(I)+DTAU(I-1))/Q(I)
        END DO
!
! HU,HL and HT are used to compute r^2.h. The equation for r^2.H is
!
! r^2.H(K+0.5) = HU(K).J(K+1) - HL(K).J(K) + HT(K)
!
	DO I=1,ND-1
	  T1=2.0D0*RECIP_CDELTAT/(CHI(I)+CHI(I+1))
	  T2=(BETA(I)+BETA(I+1))*(SIGMA(I)+SIGMA(I+1))/(R(I)+R(I+1))/(CHI(I)+CHI(I+1))
	  WMID(I)=1.0D0+T1+2.0D0*T2
	  HU(I)=F(I+1)*Q(I+1)/DTAU(I)/WMID(I)
	  HL(I)=F(I)*Q(I)/DTAU(I)/WMID(I)
	  HT(I)=T1*ROLD_ON_R*ROLD_ON_R*RSQ_H_OLDT(I)/WMID(I)
	END DO
!
	DO I=1,ND
	  JFAC(I)=AVE_DTAU_ON_Q(I)*(BETA(I)*SIGMA(I)*(1.0D0+F(I))/R(I)+RECIP_CDELTAT)/CHI(I)
	  JT(I)=AVE_DTAU_ON_Q(I)*ROLD_ON_R*ROLD_ON_R*RECIP_CDELTAT/CHI(I)
	  RAD_DECAY_TERM(I)=AVE_DTAU_ON_Q(I)*R(I)*R(I)*(E_RAD_DECAY(I)-WORK(I))/CHI(I)
	END DO
!
! Compute the TRIDIAGONAL operators, and the RHS source vector. These
! vectors are corrupted in the solution.
!
        DO I=2,ND-1
          TA(I)=HL(I-1)
          TC(I)=HU(I)
          TB(I)=-JFAC(I) - HL(I) - HU(I-1)
          XM(I)=-JT(I)*RSQ_J_OLDt(I)-(HT(I)-HT(I-1)) - RAD_DECAY_TERM(I)
        END DO
!
! Evaluate TA,TB,TC & XM for boundary conditions
!
! Outer boundary. NB. Since T1 (=W in earlier notation) multiplies H,
! we don't need it elsewhere in the equations.
!
	T1=1.0D0+RECIP_CDELTAT/CHI(1)
        TC(1)=-F(2)*Q(2)/DTAU(1)
        TB(1)= F(1)*Q(1)/DTAU(1) + H_OUTBC*T1
        XM(1)=RECIP_CDELTAT*ROLD_ON_R*ROLD_ON_R*H_OUTBC_OLDT*RSQ_J_OLDT(1)/CHI(1)
        TA(1)=0.0D0
!
! NB: H_INBC is actually R^2 . H
!
        IF(DIFF_APPROX)THEN
	  H_INBC=3.826D+13*LUMINOSITY/16.0D0/PI/PI
	  T1=1.0D0+RECIP_CDELTAT/CHI(ND)
          TA(ND)=-Q(ND-1)*F(ND-1)/DTAU(ND-1)
          TB(ND)=F(ND)/DTAU(ND-1)
          XM(ND)=H_INBC*T1 - ROLD_ON_R*ROLD_ON_R*RECIP_CDELTAT*H_INBC_OLDT/CHI(ND)
        ELSE
	  H_INBC=0.0D0
	  T1=1.0D0+RECIP_CDELTAT/CHI(ND)
          TA(ND)=-Q(ND-1)*F(ND-1)/DTAU(ND-1)
          TB(ND)=F(ND)/DTAU(ND-1)
          XM(ND)=0.0D0
        END IF
        TC(ND)=0.0D0
!
	IF(.TRUE.)THEN !VERBOSE_OUTPUT)THEN
	  WRITE(121,'(2ES16.6)')RECIP_CDELTAT,ROLD_ON_R
	  WRITE(121,'(2ES16.6)')H_INBC,H_INBC_OLDT
	  WRITE(121,'(3ES16.6)')H_OUTBC,H_OUTBC_OLDT,RSQ_J_OLDt(1)
	  DO I=1,ND
	    WRITE(121,'(I4,7ES16.6)')I,TA(I),TB(I),TC(I),JFAC(I),JT(I),XM(I),F(I)
	  END DO
	END IF
	CALL THOMAS(TA,TB,TC,XM,ND,1)
!
	DO I=1,ND-1
	  RSQ_HFLUX(I)=(HU(I)*XM(I+1)-HL(I)*XM(I)) + HT(I)
	END DO
	RJ(1:ND)=XM(1:ND)/R(1:ND)/R(1:ND)
!
	LU=122
	OPEN(UNIT=LU,FILE='SN_GREY_CHK',STATUS='UNKNOWN')
	  WRITE(LU,'(A,ES12.4)')'  RECIP_CDELTAT=',RECIP_CDELTAT
	  WRITE(LU,'(A,ES12.4)')'      R_ON_ROLD=',1.0D0/ROLD_ON_R
	  WRITE(LU,'(4X,A1,11(A12))')'I','       R','     CHI','    RSQJ','RSQJ_OLD','    RSQH','    DJDT',
	1                      '   HTERM','  MHTERM','   RSQdE','    WORK','     SUM'
	  WRITE(LU,'(I5,11ES12.4)')1,R(1),CHI(1),XM(1),ROLD_ON_R*ROLD_ON_R*RSQ_J_OLDt(1),RSQ_HFLUX(1)
	  DO I=2,ND-1
	    T1=RECIP_CDELTAT*(XM(I)-ROLD_ON_R*ROLD_ON_R*RSQ_J_OLDT(I))
	    T2=-2.0D0*Q(I)*CHI(I)*(RSQ_HFLUX(I)-RSQ_HFLUX(I-1))/(DTAU(I-1)+DTAU(I))
	    E3=R(I)*R(I)*E_RAD_DECAY(I) 
	    T3=R(I)*R(I)*WORK(I) 
	    WRITE(LU,'(I5,11ES12.4)')I,R(I),CHI(I),XM(I),ROLD_ON_R*ROLD_ON_R*RSQ_J_OLDt(I),RSQ_HFLUX(I),
	1         T1,T2,2.0D0*(RSQ_HFLUX(I-1)-RSQ_HFLUX(I))/(R(I-1)-R(I+1)),E3,T3,(T1+T2+T3-E3)
	  END DO
	  I=ND
	  WRITE(LU,'(I5,11ES12.4)')I,R(I),CHI(I),XM(I),ROLD_ON_R*ROLD_ON_R*RSQ_J_OLDt(I),RSQ_HFLUX(I)
	  WRITE(LU,'(4X,A1,11(A12))')'I','       R','     CHI','    RSQJ','RSQJ_OLD','    RSQH','    DJDT',
	1                      '   HTERM','  MHTERM','   RSQdE','    WORK','     SUM'
	  WRITE(LU,'(A)')' '
	  WRITE(LU,'(A)')' NB: RSQJ_OLD has been scaled by ROLD_ON_R^2'
	  WRITE(LU,'(A)')'     HTERM is the term in the transfer equation (dTAU)'
	  WRITE(LU,'(A)')'     MHTERM is the term in the transfer equation (CHIdR)'
	  WRITE(LU,'(A)')'     WORK is the term due to the change in gas enthalpy'
	  WRITE(LU,'(A)')'     RSQdE is the energy term due to radioactive decay'
	  WRITE(LU,'(A)')' '
	CLOSE(LU)
!
	BAD_J=.FALSE.
	DO I=1,ND
	  IF(RJ(I) .LE. 0.0D0)THEN
	    IF(RJ(I+1) .LE. 0.0D0 .AND. I .NE. 1)THEN
	      T1=R(I-1)/R(I)
	      RJ(I)=T1*T1*RJ(I-1)
	    ELSE IF(RJ(I+1) .LE. 0.0D0)THEN
	      DO J=I+2,ND
	        IF(RJ(J) .GT. 0)THEN
	          T1=R(J)/R(I)
	          RJ(I)=T1*T1*RJ(J)
	          EXIT
	        END IF
	      END DO
	    ELSE
	      RJ(I)=SQRT(RJ(I+1)*RJ(I-1))
	    END IF 
	    BAD_J=.TRUE.
	  END IF
	  TA(I)=RJ(I)+(E_RAD_DECAY(I)-WORK(I))/MAX(0.1D0*CHI(I),CHI_PLANCK(I))
	END DO
        IF(BAD_J)BAD_J_COUNTER=BAD_J_COUNTER+1
!
	IF(BAD_J_COUNTER .GT. 300)THEN
	  WRITE(6,*)'Exceeded 300 iterations for BAD_J_COUNTER'
	  WRITE(6,*)'Check starting conditions'
	  WRITE(6,*)'Check DJDT_GREY_ERROR file for more information'
	  OPEN(UNIT=LU,FILE='DJDT_GREY_ERRROR',STATUS='UNKNOWN')
	    WRITE(LU,*)'Reduction factor is',REDUCTION_FACTOR
	    WRITE(LU,'(A,9A14)')'  J','T_FROM_J','    OLD_T','E(Decay)','    WORK',
	1                        '     CHI','    CHIP','      RJ','    MOD_J'
	    DO J=1,ND
	      WRITE(LU,'(I3,8ES14.4)')J,T_FROM_J(J),OLD_T(J),E_RAD_DECAY(J),WORK(J),
	1                                CHI(J),CHI_PLANCK(J),RJ(J),TA(J)
	    END DO
	  CLOSE(UNIT=LU)
	  COMPUTED=.FALSE.
	  RETURN
	END IF
!
! Compute the temperature distribution, and the Rosseland optical depth scale.
! NB sigma=5.67E-05 and the factor of 1.0E-04 is to convert T from units of 
! K to units of 10^4 K. 
!
	DO I=1,ND
	  IF(TA(I) .LE. 0.0D0)THEN
	    IF(WORK(I) .GT. 0)THEN
	      T_FROM_J(I)=0.5D0*(OLD_T(I)*ROLD_ON_R+T_FROM_J(I))
	    ELSE
	      T_FROM_J(I)=0.5D0*(OLD_T(I)*ROLD_ON_R+T_FROM_J(I))
	    END IF
	  ELSE
            T_FROM_J(I)=0.9D0*T_FROM_J(I)+0.1D0*((PI/5.67D-05*TA(I))**0.25D0)*1.0D-04
	  END IF
        END DO
	IF(BAD_J)GOTO 1000
!
	WRITE(146,'(5ES14.4)')T_FROM_J(ND-1),WORK(ND-1),T_FROM_J(ND),WORK(ND),REDUCTION_FACTOR
	BAD_J=.FALSE.
	DO I=1,ND
	  IF(ABS(T_FROM_J(I)-TSTORE(I)) .GT. 0.01D0)BAD_J=.TRUE.
	END DO
	TSTORE(1:ND)=T_FROM_J(1:ND)
	IF(BAD_J)GOTO 1000
	IF(REDUCTION_FACTOR .GT. 1.0D0)THEN
!	  REDUCTION_FACTOR=MAX(1.0D0,REDUCTION_FACTOR-1.0D0)
	  REDUCTION_FACTOR=MAX(1.0D0,REDUCTION_FACTOR/1.2D0)
	  GOTO 1000
	END IF	
	BAD_J_COUNTER=0
!
! Output gray fluxes and boundary Eddington factors for next model in time sequence.
!
! T1 is used as a dummy argument for NU.
!  I is used as a dummy argument for NCF. Neither value is used when the
!  'GREY' option is set. The logical variable L_TRUE is also not accessed.
!
	CALL OUT_JH(XM,RSQ_HFLUX,H_INBC,H_OUTBC,T1,I,R,VEL,ND,L_TRUE,'GREY')
!
! 
!
5000	CONTINUE
!
! Solve for the Eddington factors using a ray by ray solution. The D/Dt terms
! are presently ignored.
!
	FS_RSQJ(:)=0.0D0
	FS_RSQH(:)=0.0D0
	FS_RSQK(:)=0.0D0
	H_OUTBC=0.0D0
!
! DBB =3L/16(piR)**2 and is used for the lower boundary diffusion approximation.
!
	T1=3.826D+13*LUMINOSITY/16.0D0/PI/PI
	DBB=3.0D0*T1/R(ND)/R(ND)
!
	CALL DERIVCHI(dCHIdr,CHI,R,ND,METHOD)
!
! Enter loop for each impact parameter P
!
	DO 2000 LS=1,NP
	  NI=ND-(LS-NC-1)
	  IF(LS .LE. NC)THEN
	    NI=ND
	    DBC=DBB*DSQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)
	  END IF
!
! Compute Z for this impact parameter
!
	  IF(NI .GT. 1)THEN
	    DO I=1,NI
	      Z(I)=DSQRT( (R(I)-P(LS))*(R(I)+P(LS)) )
	      T1=Z(I)/R(I)
	      TA(I)=BETA(I)/R(I)
	      CHI_MOD(I)=CHI(I)+TA(I)
	    END DO
	    IF(NI .EQ. 2)THEN
	      TB(1)=0.0D0         !Cancel so values unimportant
	      TB(2)=TB(1)
	    ELSE
	      CALL DERIVCHI(TB,TA,R,NI,METHOD)
	    END IF
	    DO I=1,NI
	      dCHI_MODdR(I)=dCHIdR(I)+TB(I)
	    END DO
	  END IF
!
! Compute Z for this impact parameter
!
	  IF(NI .GT. 2)THEN
	    DO I=1,NI-1
	      DTAU(I)=0.5D0*(Z(I)-Z(I+1))*(CHI_MOD(I)+CHI_MOD(I+1)+(Z(I)-Z(I+1))
	1       *(dCHI_MODdR(I+1)*Z(I+1)/R(I+1)-dCHI_MODdR(I)*Z(I)/R(I))/6.0D0)
	    END DO
	  END IF
!
! Compute TA (tridiagonal matrix) and XM vector. Since it is a gray
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
	      XM(I)=-(RJ(I)*CHI(I)+E_RAD_DECAY(I)-WORK(I))*(DTAU(I-1)+DTAU(I))*0.5D0/CHI_MOD(I)
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
!
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
	  END DO
!
	  H_OUTBC=H_OUTBC+JQW(1,LS)*(XM(1)-IBOUND)*Z(1)/R(1)
!  
2000	CONTINUE
!
! Compute the new Feautrier factors. These are stored in FS_RSQK so as not
! to destroy the old factors.
!
	DO I=1,ND
	  FS_RSQK(I)=FS_RSQK(I)/FS_RSQJ(I)
	END DO
!
! Compute the factor for the outer boundary condition.
!
	H_OUTBC=H_OUTBC/FS_RSQJ(1)
!
	T1=MAXVAL( ABS(F(1:ND)-FS_RSQK(1:ND)) )
	WRITE(LUWARN,*)'Current accuracy in JGREY_HUB_DDT_V3 is', T1
	F_LOOP_COUNTER=F_LOOP_COUNTER+1
	IF(F_LOOP_COUNTER .GT. 200)THEN
	   WRITE(6,*)'Exceeded 200 iterations or F in JGREY_HUB_DDT_V3'
	   WRITE(6,*)'Check starting conditions'
	   COMPUTED=.FALSE.
	   RETURN
	ELSE IF(T1 .GT. ACCURACY)THEN
	  F(1:ND)=FS_RSQK(1:ND)
	  GOTO 1000
	END IF
!
	OPEN(UNIT=LU,FILE='SN_GREY_CHK',STATUS='UNKNOWN', POSITION='APPEND')
	WRITE(LU,'(A,13X,A,11X,A,9X,A,13(A))')'    I','R','V','CHI','  CHI_PLANCK','        V/cR',
	1           '   RSQE(rad)','        RSQW','      RSQJ  ','RSQ(W+E)/CHI','   RSQJ(old)',
	1           '        RSQH','   RSQH(old)',
	1           '        TMOD','        T(B)',
	1           '          TR','     TR(old)'
	T1=1.0D-04*(3.14159265D0/5.67D-05)**0.25D0
	DO I=1,ND
	  T2=R(I)-1.0D-05*DELTA_TIME_SECS*VEL(I)
	  T2=RSQ_J_OLDt(I)/T2/T2
	  T3=T1*( RJ(I)+(E_RAD_DECAY(I)-WORK(I))/MAX(0.1D0*CHI(I),CHI_PLANCK(I)) )**0.25D0
	  WRITE(LU,'(I5,ES14.5,15ES12.4)')I,R(I),VEL(I),CHI(I),CHI_PLANCK(I),BETA(I)/R(I),
	1       R(I)*R(I)*E_RAD_DECAY(I),R(I)*R(I)*WORK(I),R(I)*R(I)*RJ(I),
	1       R(I)*R(I)*(WORK(I)+E_RAD_DECAY(I))/CHI(I),RSQ_J_OLDt(I),
	1       RSQ_HFLUX(I),RSQ_H_OLDt(I),
	1       POPS(NT,I),T3,T1*(RJ(I)**0.25D0),T1*(T2 **0.25D0) 
	END DO
	WRITE(LU,'(A,13X,A,11X,A,9X,A,13(A))')'    I','R','V','CHI','  CHI_PLANCK','        V/cR',
	1           '   RSQE(rad)','        RSQW','      RSQJ  ','RSQ(W+E)/CHI','   RSQJ(old)',
	1           '        RSQH','   RSQH(old)',
	1           '        TMOD','        T(B)',
	1           '          TR','      TR(old)'
!
        IF(DIFF_APPROX)THEN
	  T1=1.0D0+RECIP_CDELTAT/CHI(ND)
          TA(ND)=-Q(ND-1)*F(ND-1)/DTAU_INB_SAVE
          TB(ND)=F(ND)/DTAU_INB_SAVE
          XM(ND)=H_INBC*T1 - ROLD_ON_R*ROLD_ON_R*RECIP_CDELTAT*H_INBC_OLDT/CHI(ND)
	  I=ND
	  T1=TA(I)*RJ(I-1)*R(I-1)*R(I-1)+TB(I)*RJ(I)*R(I)*R(I)
	  WRITE(LU,'(A,13X,A,12X,A,10X,A,9(1X,A))')'    I','R','V','CHI','        V/cR',
	1           '   RSQE(rad)','        RSQJ','      H_INBC','  H_INBC(old)',
	1           '         LHS','         RHS','         SUM'
	  WRITE(LU,'(I5,ES14.5,11ES13.4)')I,R(I),VEL(I),CHI(I),BETA(I)/R(I),
	1       R(I)*R(I)*E_RAD_DECAY(I),R(I)*R(I)*RJ(I),H_INBC,H_INBC_OLDT,
	1       T1,XM(ND),T1-XM(ND)
	END IF
	CLOSE(UNIT=LU)
!
! We really want to return B which is J+E(radio decay)/CHI
!
	DO I=1,ND
	  RJ(I)=RJ(I)+(E_RAD_DECAY(I)-WORK(I))/MAX(0.1D0*CHI(I),CHI_PLANCK(I))
	END DO
!
	RETURN
	END
