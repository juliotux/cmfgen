C
C Subroutine to copute the mean intensity JBAR and then net radiative
C bracket ZNET. The F and G gaunt factors must be supplied. Q is
C computed internally. Two additional vectors  are returned. These
C contain flux information, and information for the constraint of
C radiative equilibrium. It is assumed that the source function  is
C frequency independent.
C
C The particular choice of boundary adopted is irrelevant for this
C routine. Such information is incorporated by the outer boundary
C Eddington factors HBC, and NBC.
C
C The program computes JNU and r^2 HNU. Initially the JNU vector
C also contains r^2 . JNU but this is corrected before leaving the
C program.
C
	SUBROUTINE MOMJBAR(ETA,CHI,ESEC,THETA,JCONT,CHIL,ETAL,
	1                  V,SIGMA,R,JBAR,ZNET,
	1                  JNU,RSQHNU,F,G,HBC,IN_HBC,NBC,JBLANK,HBLANK,
	1                  PF,PROF,LFQW,FL,DIF,DBB,IC,METHOD,
	1                  EW,CONT_INT,LINE_BL,FULL_ES,
	1                  NLF,NC,NP,ND)
	IMPLICIT NONE
C
C Altered 05-Jun-1996 - LFQW(ML)*PROF(ML) set to T1 and removed from inner loop
C                         when evaluating JBAR(I) because of prolem with CRAY
C                         F90 compiler.
C Altered 24-May-1996 - SCRTEMP blk removed.
C                       Dynamic memory allocation now used (NV deleted).
C                       Calls to DP_ZERO removed.
C
C Altered  6-Jun-1989 - TCHIPREV removed as experience suggests that
C                       it makes little difference to the solution accuracy.
C                       Its inclusion would complicate the linearization.
C                       Cleaning. DTAUONQ vector included. For ML=1,
C                       PROF is now assumed to be zero when evaluating the
C                       opacity and souce function.
C Altered 6-Jun-1989 - Now divide diffusion correction by TCHI (was
C                      errantly dividing by CHI).
C Altered 16-May-1989 - FULL_ES option installed. This option implies
C                       that photons scattered by an electron are not
C                       absorbed by the line.
C Altered 12-May-1989 - Cleaned. Small bug fixes. HBLANK and JBLANK
C                       are now returned with the continuum subtracted out.
C Finalized and Tested 28-Apr-1989.
C Created 10-JAN-1989
C
	INTEGER NLF,NC,NP,ND
	REAL*8 ETA(ND),CHI(ND),THETA(ND),ESEC(ND)
	REAL*8 CHIL(ND),ETAL(ND),JCONT(ND)
	REAL*8 V(ND),SIGMA(ND),R(ND)
C
C Radiation field variables.
C
	REAL*8 JNU(ND,NLF+1),RSQHNU(ND,NLF+1)
	REAL*8 F(ND,NLF+1),G(ND,NLF+1)
	REAL*8 HBC(3,NLF+1),NBC(3,NLF+1),IN_HBC(NLF+1)
	REAL*8 JBAR(ND),ZNET(ND),HBLANK(ND),JBLANK(ND)
C
C Profile information
C
	REAL*8 PF(NLF),PROF(NLF),LFQW(NLF)
C
	REAL*8 DBB,IC,FL,EW,CONT_INT,SRCEBND
	CHARACTER*6 METHOD
	LOGICAL DIF,LINE_BL,FULL_ES
C
	REAL*8 TA(ND),TB(ND),TC(ND),DTAU(ND),TCHI(ND),DTAUONQ(ND)
	REAL*8 Q(ND),XM(ND),SOURCE(ND),JEX_SCAT(ND)
	REAL*8 VB(ND),VC(ND),HU(ND),HL(ND),HS(ND)
	REAL*8 GAM(ND),GAMH(ND),W(ND),WPREV(ND)
	REAL*8 PSI(ND),PSIPREV(ND)
C
C Local variables.
C
	INTEGER, PARAMETER :: IONE=1
	INTEGER I,ML
	REAL*8 T1,T2,DNU,SCALE
C
C Zero vecters. These were originally zerroed in a common block. There is
C currently 21 vectors.
C
	TA(:)=0.0D0;        TB(:)=0.0D0;      TC(:)=0.0D0
	DTAU(:)=0.0D0;      TCHI(:)=0.0D0;    DTAUONQ(:)=0.0D0
	Q(:)=0.0D0;         XM(:)=0.0D0;      SOURCE(:)=0.0D0
        JEX_SCAT(:)=0.0D0;  VB(:)=0.0D0;      VC(:)=0.0D0
        HU(:)=0.0D0;        HL(:)=0.0D0;      HS(:)=0.0D0
	GAM(:)=0.0D0;       GAMH(:)=0.0D0;    W(:)=0.0D0;
	WPREV(:)=0.0D0;     PSI(:)=0.0D0;     PSIPREV(:)=0.0D0
C
C 
C
C Zero relevant vectors and matrices.
C
	JBAR(:)=0.0D0
	HBLANK(:)=0.0D0
	JBLANK(:)=0.0D0
	JNU(:,:)=0.0D0			!ND*(NLF+1)
	RSQHNU(:,:)=0.0D0		!ND*(NLF+1) )
C
C*****************************************************************************
C
	DO ML=1,NLF
C
C Compute the total opacity. Store opacity from previous frequency.
C Compute the total source function.
C
C                 *************************
C                 *************************
C It is currently assumed that ETA is the total continuum emssivity
C and hence already contains the continuum scattering term. This may
C need to be altered.
C                 *************************
C                 *************************
C
	  IF(ML .EQ. 1)THEN
	    DO I=1,ND
	      TCHI(I)=CHI(I)
	      SOURCE(I)=ETA(I)/CHI(I)
	    END DO
	  ELSE
	    DO I=1,ND
	      TCHI(I)=CHI(I)+CHIL(I)*PROF(ML)
	      SOURCE(I)=(ETA(I)+ETAL(I)*PROF(ML))/TCHI(I)
	    END DO
	  END IF
C
C Compute the Q factors from F. Then compute optical depth scale.
C
	  CALL QFROMF(F(1,ML),Q,R,TA,TB,ND)	!TA work vector
	  DO I=1,ND
	    TA(I)=TCHI(I)*Q(I)
	  END DO
	  CALL DERIVCHI(TB,TA,R,ND,METHOD)
	  CALL NORDTAU(DTAU,TA,R,R,TB,ND)
C
C Assume (1)	SIGMAd+1/2 = 0.5*( SIGMAd+1+SIGMAd )
C 	 (2)	Vd+1/2=0.5*( Vd + Vd+1 )
C Note that V is in km/s and SIGMA=(dlnV/dlnR-1.0)
C
	  IF(ML .NE. 1)THEN
	    DNU=PF(ML-1)-PF(ML)
C
C NB - By definition, G is defined at the mesh midpoints.
C
	    DO I=1,ND-1
	      GAMH(I)=2.0D0*3.33564D-06*(V(I)+V(I+1))/(R(I)+R(I+1))
	1         /DNU/( TCHI(I)+TCHI(I+1) )
	      W(I)=GAMH(I)*( 1.0D0+0.5D0*(SIGMA(I)+SIGMA(I+1))*G(I,ML) )
	      WPREV(I)=GAMH(I)*( 1.0D0+
	1                0.5D0*(SIGMA(I)+SIGMA(I+1))*G(I,ML-1) )
	    END DO
C
	    DO I=1,ND
	      GAM(I)=3.33564D-06*V(I)/R(I)/TCHI(I)/DNU
	    END DO
C
C PSIPREV is equivalent to the U vector of FORMSOL.
C
	    PSI(1)=GAM(1)*(HBC(1,ML)+NBC(1,ML)*SIGMA(1))
	    PSIPREV(1)=GAM(1)*(HBC(1,ML-1)+NBC(1,ML-1)*SIGMA(1))
	    SRCEBND=GAM(1)*( HBC(2,ML)+NBC(2,ML)*SIGMA(1)
	1                   -HBC(2,ML-1)-NBC(2,ML-1)*SIGMA(1) )
	    DO I=2,ND
	      DTAUONQ(I)=0.5D0*(DTAU(I)+DTAU(I-1))/Q(I)
	      PSI(I)=DTAUONQ(I)*GAM(I)*( 1.0D0+SIGMA(I)*F(I,ML) )
	      PSIPREV(I)=DTAUONQ(I)*GAM(I)*(  1.0D0+SIGMA(I)*F(I,ML-1) )
	    END DO
C
	  ELSE
	    DO I=2,ND
	      DTAUONQ(I)=0.5D0*(DTAU(I)+DTAU(I-1))/Q(I)
	    END DO
	    DO I=1,ND
	      GAMH(I)=0.0D0
	      W(I)=0.0D0
	      WPREV(I)=0.0D0
	      GAM(I)=0.0D0
	      PSI(I)=0.0D0
	      PSIPREV(I)=0.0D0
	    END DO
	    SRCEBND=0.0D0
	  END IF
C 
C
C Compute vectors used to compute the flux vector H.
C
	  DO I=1,ND-1
	    HU(I)=F(I+1,ML)*Q(I+1)/(1.0D0+W(I))/DTAU(I)
	    HL(I)=F(I,ML)*Q(I)/(1.0D0+W(I))/DTAU(I)
	    HS(I)=WPREV(I)/(1.0D0+W(I))
	  END DO
C
C Compute the TRIDIAGONAL operators, and the RHS source vector.
C
	  DO I=2,ND-1
	    TA(I)=-HL(I-1)
	    TC(I)=-HU(I)
	    TB(I)=DTAUONQ(I) + PSI(I) + HU(I-1) + HL(I)
	    VB(I)=-HS(I-1)
	    VC(I)=HS(I)
	    XM(I)=DTAUONQ(I)*SOURCE(I)*R(I)*R(I)
	  END DO
C
C Evaluate TA,TB,TC for boudary conditions
C
	  TC(1)=-F(2,ML)*Q(2)/DTAU(1)
	  TB(1)=F(1,ML)*Q(1)/DTAU(1) + HBC(1,ML) + PSI(1)
	  XM(1)=R(1)*R(1)*(HBC(2,ML)+SRCEBND)
	  TA(1)=0.0D0
	  VB(1)=0.0D0
	  VC(1)=0.0D0
C
	  TA(ND)=-F(ND-1,ML)*Q(ND-1)/DTAU(ND-1)
	  IF(DIF)THEN
	    TB(ND)=F(ND,ML)/DTAU(ND-1)
	    XM(ND)=DBB*R(ND)*R(ND)/3.0D0/TCHI(ND)
	  ELSE
	    TB(ND)=F(ND,ML)/DTAU(ND-1)+IN_HBC(ML)
	    XM(ND)=R(ND)*R(ND)*IC*(0.25D0+0.5D0*IN_HBC(ML))
	  END IF
	  TC(ND)=0.0D0
	  VB(ND)=0.0D0
	  VC(ND)=0.0D0
	  PSIPREV(ND)=0.0D0
C
	  IF(ML .NE. 1)THEN
	    XM(1)=XM(1) + PSIPREV(1)*JNU(1,ML-1)
	    DO I=2,ND-1
	      XM(I)=XM(I) + VB(I)*RSQHNU(I-1,ML-1) + VC(I)*RSQHNU(I,ML-1)
	1          + PSIPREV(I)*JNU(I,ML-1)
	    END DO
	     XM(ND)=XM(ND)
	  END IF
C
C Solve for the radiation field along ray for this frequency.
C
	  CALL THOMAS(TA,TB,TC,XM,ND,IONE)
C
	  DO I=1,ND
	    JNU(I,ML)=XM(I)
	  END DO
C
	  IF(ML .EQ. 1)THEN
	    DO I=1,ND-1
	      RSQHNU(I,ML)=HU(I)*XM(I+1)-HL(I)*XM(I)
	    END DO
	  ELSE
	    DO I=1,ND-1
	      RSQHNU(I,ML)=HU(I)*XM(I+1)-HL(I)*XM(I)+HS(I)*RSQHNU(I,ML-1)
	    END DO
	  END IF
	END DO
C
C 
C
C*****************************************************************************
C
C Compute JBAR.
C
	DO ML=1,NLF
	  T1=LFQW(ML)*PROF(ML)
	  DO I=1,ND
	    JBAR(I)=JBAR(I)+T1*JNU(I,ML)
	  END DO
	END DO
C
	DO ML=1,NLF
	  DO I=1,ND
	    HBLANK(I)=HBLANK(I)+RSQHNU(I,ML)*LFQW(ML)
	    JBLANK(I)=JBLANK(I)+JNU(I,ML)*LFQW(ML)
	  END DO
	END DO
C
C Scale JBLANK and HBLANK to allow for the fact that the
C LFQW dont necessarily add to v*(PF(1)-PF(NLF)) (since LFQW is
C normalized so that integral over the line profile is unity).
C
	SCALE=0.0D0
	DO ML=1,NLF
	  SCALE=SCALE+LFQW(ML)
	END DO
	SCALE=1.0D+15*FL*(PF(1)-PF(NLF))/SCALE
C
C Note that the matrice JNU and RSQHNU contain r^2.J and r^2.H respectively.
C Also put both the HBLANK and JBLANK in the correct units - need to multiply
C by 10^15 . FL because of the frequency integration.
C Correct JNU so that mean intensity is returned.
C
	DO I=1,ND
	  T1=R(I)*R(I)
	  JBAR(I)=JBAR(I)/T1
	  JBLANK(I)=JBLANK(I)*(SCALE/T1)
	  HBLANK(I)=HBLANK(I)*SCALE
	  ZNET(I)=1.0D0-JBAR(I)*CHIL(I)/ETAL(I)
	  DO ML=1,NLF
	    JNU(I,ML)=JNU(I,ML)/T1
	  END DO
	END DO
C
	IF(FULL_ES)THEN
	  T1=1.0D-15/FL
	  T2=PF(1)-PF(NLF)
	  DO I=1,ND
	    JEX_SCAT(I)=JBLANK(I)*T1-JCONT(I)*T2
	  END DO
	ELSE
	  JEX_SCAT(:)=0.0D0
	END IF
C
C 
C
C This section of the routine provides a means of estimating the
C line blanketing for selected lines. Once where on the red side of the
C line profile the opacity is constant, hence we can integrate
C the transfer equation with respect to frequency. This leaves
C us with a simple "continuum like" transfer problem. To
C do the frequency integration it is also necessary to assume that
C the source function is frequency independent. We have thus assumed
C that the line photons are scattered coherently. This should be an
C excellent approximation since the continuum transfer of the line
C photons should not be strongly influenced by the photons frequency.
C
C
	IF(LINE_BL)THEN
C
C DNU is the frequency bandwidth over which we integrate to obtain
C JINT and is determined by the maximum expansion velocity of the
C atmosphere. We define DNU with out the factor FL. The bandwidth
C must be identical to that in FGCOMP.
C
	  DNU=3.33564D-06*V(1)*2.0D0
C
C Compute the total opacity. Store opacity from previous frequency.
C Compute the total source function. We will assume that for the
C the propopgation of the line photons that the electron scattering
C is coherent. Since ETA contains a continuum scattering term, we
C first need to subtract this out since it will automatically
C be included in the computations. NB - JEX_SCAT is non-zero only
C if FULL_ES has been specified.
C
	  DO I=1,ND
	    SOURCE(I)=( (ETA(I)-ESEC(I)*JCONT(I))*DNU +
	1                   JEX_SCAT(I)*ESEC(I) )/CHI(I)
	  END DO
C
C Compute the Q factors from F. Then compute optical depth scale.
C
	  CALL QFROMF(F(1,NLF+1),Q,R,TA,TB,ND)	!TA work vector
	  DO I=1,ND
	    TA(I)=CHI(I)*Q(I)
	  END DO
	  CALL DERIVCHI(TB,TA,R,ND,METHOD)
	  CALL NORDTAU(DTAU,TA,R,R,TB,ND)
C
C Assume (1)	SIGMAd+1/2 = 0.5*( SIGMAd+1+SIGMAd )
C 	 (2)	Vd+1/2=0.5*( Vd + Vd+1 )
C Note that V is in km/s and SIGMA=(dlnV/dlnR-1.0)
C
	  DO I=1,ND-1
	    GAMH(I)=2.0D0*3.33564D-06*(V(I)+V(I+1))/(R(I)+R(I+1))
	1         /( CHI(I)+CHI(I+1) )
	    W(I)=1.0D0+0.5D0*( SIGMA(I)+SIGMA(I+1) )*G(I,NLF)
	    WPREV(I)=1.0D0+0.5D0*( SIGMA(I)+SIGMA(I+1) )*G(I,1)
	  END DO
C
	  DO I=1,ND
	    GAM(I)=3.33564D-06*V(I)/R(I)/CHI(I)
	  END DO
C
C Compute RHS of tridiagonal system of equations. Note that the
C SOURCE vector has been already corrected for the integration over
C DNU. NB - JNU at this point contains only J, but RSQHNU contains
C r^2H.
C
	  XM(1)=GAM(1)*R(1)*R(1)*( (HBC(1,1)+SIGMA(1)*NBC(1,1))*JNU(1,1)
	1                         -(HBC(2,1)+SIGMA(1)*NBC(2,1))
	1                - (HBC(1,NLF)+SIGMA(1)*NBC(1,NLF))*JNU(1,NLF)
	1                         +(HBC(2,NLF)+SIGMA(1)*NBC(2,NLF)) )
	  DO I=2,ND-1
	    XM(I)=R(I)*R(I)*(  SOURCE(I)
	1        +GAM(I)*( (1.0D0+SIGMA(I)*F(I,NLF))*JNU(I,NLF)
	1        -(1.0D0+SIGMA(I)*F(I,1))*JNU(I,1) )  )/Q(I)
	1        +(  GAMH(I)*( W(I)*RSQHNU(I,NLF)-WPREV(I)*RSQHNU(I,1) )
	1        -GAMH(I-1)*( W(I-1)*RSQHNU(I-1,NLF)
	1        -WPREV(I-1)*RSQHNU(I-1,1) )  )
	1        *2.0D0/(DTAU(I)+DTAU(I-1))
	  END DO
C
C Note well - DBB =dB/dR (and Q(ND)=1.0 by definition)
C
	  IF(DIF)THEN
	    XM(ND)=DNU*R(ND)*R(ND)*DBB/CHI(ND)/3.0D0
	  ELSE
	    XM(ND)=DNU*R(ND)*R(ND)*IC*( 0.25D0+0.5D0*IN_HBC(NLF+1) )
	  END IF
C
C Compute T ( a tridiagonal matrix) and store it as three vectors
C TA,TB and TC .
C
	  CALL TFEAU(TA,TB,TC,R,Q,F(1,NLF+1),THETA,DTAU
	1             ,HBC(1,NLF+1),IN_HBC(NLF+1),DIF,ND)
C
C Find the solution.
C
C After soultion, XM is integral [ jnu ] (no r^2 since used TFEAU routine.
C
	  CALL THOMAS(TA,TB,TC,XM,ND,IONE)
	  DO I=1,ND
	    JNU(I,NLF+1)=XM(I)
	  END DO
C
C Compute the H vector (=r^2 H).
C
	  DO I=1,ND-1
	    RSQHNU(I,NLF+1)=( R(I+1)*R(I+1)*F(I+1,NLF+1)*Q(I+1)*XM(I+1)
	1      -R(I)*R(I)*F(I,NLF+1)*Q(I)*XM(I) )/DTAU(I)
	1      +GAMH(I)*( W(I)*RSQHNU(I,NLF)-WPREV(I)*RSQHNU(I,1) )
	  END DO
C 
C
C Evaluate JBLANK which is deined by Int{Jv dv} over the WHOLE line.
C HBLANK is similarly defined. Note that JBLANK and HBLANK have already
C been scaled. The frequency factor of 10^15 is included in the JBLANK
C (and HBLANK) definition.
C
	  T1=1.0D+15*FL
	  DO I=1,ND
	    JBLANK(I)=JBLANK(I)+JNU(I,NLF+1)*T1
	    HBLANK(I)=HBLANK(I)+RSQHNU(I,NLF+1)*T1
	  END DO
C
C Evaluate the line EW. The units are Angstroms. Also evaluate
C the continuum intensity ( Jys/kpc/kpc ). Note that H is
C defined midway between R(1) and R(2).
C
	  T1=( (PF(1)-PF(NLF))+DNU )*FL*1.0D+15
	  EW=2.99794D-12*( HBLANK(1)-RSQHNU(1,1)*T1 )/
	1                 RSQHNU(1,1)/FL/FL
	  CONT_INT=13.19868D0*RSQHNU(1,1)
C
C Change HBLANK and JBLANK to be the int{across line -Jc} (i.e
C integral of J or H above the continuum). Thus, for a weak line,
C HBLANK and JBLANK should be zero.
C
	  DO I=1,ND
	    HBLANK(I)=HBLANK(I)-RSQHNU(I,1)*T1
	    JBLANK(I)=JBLANK(I)-JNU(I,1)*T1
	  END DO
C
C JNU(I,NLF+1) is only used by FG_HAM to determine the electron scattering
C source function. Thus if FULL_ES is specified, we can include the extra
C intensity directly into the definition of JNU( , NLF+1).
C
	  DO I=1,ND
	    JNU(I,NLF+1)=JNU(I,NLF+1)+JEX_SCAT(I)
	  END DO
	END IF
C
	RETURN
	END
