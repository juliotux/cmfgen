C
C Created 08-May-1989 - Subroutine to copute the mean intensity JBAR
C                       and then net radiative bracket ZNET. Second
C                       Order differencing in frequency is used. The F
C                       and G gaunt factors must be supplied. Q is
C                       computed internally. Two additional vectors
C                       are returned. These contain flux information,
C                       and information for the constraint of radiative
C                       equilibrium. It is assumed that the source function
C                       is frequency independent.
C
C The particular choice of the outer boundary condition adopted is irrelevant
C for this routine. Such information is incorporated by the outer boundary
C Eddington factors HBC, and NBC.
C
C The program computes JNU and r^2 HNU. Initially the JNU vector
C also contains r^2 . JNU but this is corrected before leaving the
C program.
C
C Created 05-May-1989
C Altered 12-May-1989 - Several small bug fixes. JBLANK and HBLANK
C                       now have continuum subtracted out.
C Altered 15-May-1989 - Bug fix. HS was incorrect.
C Altered 16-May-1989 - FULL_ES option installed. This option implies
C                       that photons scattered by an electron are not
C                       absorbed by the line.
C
	SUBROUTINE MOMHAM(ETA,CHI,ESEC,THETA,JCONT,CHIL,ETAL,
	1                  V,SIGMA,R,JBAR,ZNET,
	1                  JNU,RSQHNU,F,G,HBC,IN_HBC,NBC,JBLANK,HBLANK,
	1                  PF,PROF,LFQW,FL,DIF,DBB,IC,METHOD,
	1                  EW,CONT_INT,LINE_BL,FULL_ES,
	1                  NLF,NC,NP,ND)
	IMPLICIT NONE
C
	INTEGER NLF,NC,NP,ND,NV
	PARAMETER (NV=100)
	REAL*8 ETA(ND),CHI(ND),THETA(ND),ESEC(ND)
	REAL*8 CHIL(ND),ETAL(ND),JCONT(ND)
	REAL*8 V(ND),SIGMA(ND),R(ND)
C
C Radiation field variables.
C
	REAL*8 JNU(ND,NLF+1),RSQHNU(ND,NLF+1)
	REAL*8 F(ND,NLF+1),G(ND,NLF+1)
	REAL*8 HBC(NLF+1),NBC(NLF+1),IN_HBC(NLF+1)
	REAL*8 JBAR(ND),ZNET(ND),HBLANK(ND),JBLANK(ND)
C
C Profile information
C
	REAL*8 PF(NLF),PROF(NLF),LFQW(NLF)
C
	REAL*8 DBB,IC,FL,EW,CONT_INT
	CHARACTER*6 METHOD
	LOGICAL DIF,LINE_BL,FULL_ES
C
	COMMON /SCRATCH/  PROGDESC,TA,TB,TC,DTAU,TCHI,TCHIPREV,
	1                   XM,SOURCE,MIDF,Q,JEX_SCAT,
	1                   UA,UB,UC,VB,VC,
	1                   HU,HL,HS,GAM,GAMH,W,WPREV,PSI,PSIPREV
	REAL*8 TA(NV),TB(NV),TC(NV),DTAU(NV),TCHI(NV),TCHIPREV(NV)
	REAL*8 XM(NV),SOURCE(NV),MIDF(NV),Q(NV),JEX_SCAT(NV)
	REAL*8 UA(NV),UB(NV),UC(NV)
	REAL*8 VB(NV),VC(NV),HU(NV),HL(NV),HS(NV)
	REAL*8 GAM(NV),GAMH(NV),W(NV),WPREV(NV)
	REAL*8 PSI(NV),PSIPREV(NV)
C
	REAL*8 PROGDESC	
C
C Local variables.
C
	INTEGER I,ML
	REAL*8 T1,T2,T3,DNU,SCALE,MID_PRO
C
C PROGDESC is a variable use to confirm that the scratch block is not
C being used by some other routine.
C
	PROGDESC=236742249.0D0		!Must be unique
C
	IF(ND .GT. NV)THEN
	  WRITE(6,*)'Error in MOMHAM - NV smaller than ND'
	  WRITE(6,*)'ND=',ND,'NV',NV
	  STOP
	END IF
C
C Zero common block. There are currently 25 vectors in the common block.
C TA must be the first vector, and PSIPREV the last.
C
	PSIPREV(NV-1)=1.0D0
	PSIPREV(NV)=1.0D0
	I=(NV*25)-1
	CALL DP_ZERO(TA,I)
	IF(PSIPREV(NV-1) .NE. 0 .AND. PSIPREV(NV) .NE. 1)THEN
	  WRITE(6,*)'Error in zeroing SCRATCH block in MOMJBAR'
	  STOP
	ELSE
	  PSIPREV(NV)=0.0D0
	END IF
C
C
C
C Zero relevant vectors and matrices.
C
	DO I=1,ND
	  JBAR(I)=0.0D0
	  HBLANK(I)=0.0D0
	  JBLANK(I)=0.0D0
	END DO
	CALL DP_ZERO(JNU, ND*(NLF+1) )
	CALL DP_ZERO(RSQHNU, ND*(NLF+1) )
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
	    MID_PRO=0.0D0
	  ELSE
	    MID_PRO=0.50D0*( PROF(ML)+PROF(ML-1) )
	  END IF
	  DO I=1,ND
	    TCHI(I)=CHI(I)+CHIL(I)*MID_PRO
	    SOURCE(I)=(ETA(I)+ETAL(I)*MID_PRO)/TCHI(I)
	  END DO
C
C Compute the Q factors from F. Then compute optical depth scale.
C
	  IF(ML .EQ. 1)THEN
	    DO I=1,ND
	      MIDF(I)=F(I,1)
	    END DO
	  ELSE
	    DO I=1,ND
	      MIDF(I)=(F(I,ML-1)+F(I,ML))*0.5D0
	    END DO
	  END IF
	  CALL QFROMF(MIDF,Q,R,TA,TB,ND)	!TA work vector
	  DO I=1,ND
	    TA(I)=TCHI(I)*Q(I)
	  END DO
	  IF(METHOD .EQ. 'ZERO')THEN
	    CALL TAU(DTAU,TA,R,ND)
	  ELSE
	    CALL DERIVCHI(TB,TA,R,ND,METHOD)
	    CALL NORDTAU(DTAU,TA,R,R,TB,ND)
	  END IF
C
C
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
	      W(I)=2.0D0*GAMH(I)
	1               *( 1.0D0+0.5D0*(SIGMA(I)+SIGMA(I+1))*G(I,ML) )
	      WPREV(I)=2.0D0*GAMH(I)
	1               *( 1.0D0+0.5D0*(SIGMA(I)+SIGMA(I+1))*G(I,ML-1) )
	    END DO
C
	    DO I=1,ND
	      GAM(I)=3.33564D-06*V(I)/R(I)/TCHI(I)/DNU
	    END DO
	    DO I=2,ND-1
	      T1=GAM(I)*(DTAU(I-1)+DTAU(I))/Q(I)
	      PSI(I)=T1*( 1.0D0+SIGMA(I)*F(I,ML) )
	      PSIPREV(I)=T1*(  1.0D0+SIGMA(I)*F(I,ML-1) )
	    END DO
	  ELSE
	    DO I=1,ND
	      GAMH(I)=0.0D0
	      W(I)=0.0D0
	      WPREV(I)=0.0D0
	      GAM(I)=0.0D0
	      PSI(I)=0.0D0
	      PSIPREV(I)=0.0D0
	    END DO
	  END IF
C
C
C Compute vectors used to compute the flux vector H.
C
	  DO I=1,ND-1
	    HU(I)=MIDF(I+1)*Q(I+1)/(1.0D0+W(I))/DTAU(I)
	    HL(I)=MIDF(I)*Q(I)/(1.0D0+W(I))/DTAU(I)
	    HS(I)=(WPREV(I)-1.0D0)/(1.0D0+W(I))
	  END DO
C
C Compute the TRIDIAGONAL operators, and the RHS source vector.
C UA, UB, UC, VB AND VC are identically zero when ML=0 (i.e. for
C continuum calculation), and are not accessed.
C
	  IF(ML .EQ. 1)THEN
	    DO I=2,ND-1
	      T3=-0.5D0*( DTAU(I)+DTAU(I-1) )/Q(I)
	      TA(I)=HL(I-1)
	      TC(I)=HU(I)
	      TB(I)=T3-(HL(I)+HU(I-1))
	      XM(I)=SOURCE(I)*R(I)*R(I)*T3
	    END DO
C
C Evaluate TA,TB,TC for boudary conditions
C
	    TC(1)=MIDF(2)*Q(2)/DTAU(1)
	    TB(1)=-MIDF(1)*Q(1)/DTAU(1)-HBC(ML)
	    TA(1)=0.0D0
	    XM(1)=0.0D0
C
	    TC(ND)=0.0
	    TA(ND)=-MIDF(ND-1)*Q(ND-1)/DTAU(ND-1)
	    IF(DIF)THEN
	      TB(ND)=MIDF(ND)/DTAU(ND-1)
	      XM(ND)=DBB*R(ND)*R(ND)/3.0D0/CHI(ND)
	    ELSE
	      TB(ND)=MIDF(ND)/DTAU(ND-1)+IN_HBC(ML)
	      XM(ND)=R(ND)*R(ND)*IC*( 0.25D0+0.5D0*IN_HBC(ML) )
	    END IF
	  ELSE
	    DO I=2,ND-1
	      T3=0.5D0*(DTAU(I)+DTAU(I-1))/Q(I)
	      TA(I)=HL(I-1)
	      TC(I)=HU(I)
	      TB(I)=-T3-PSI(I)-(HL(I)+HU(I-1))
	      UA(I)=-HL(I-1)
	      UC(I)=-HU(I)
	      UB(I)=(HL(I)+HU(I-1))+T3-PSIPREV(I)
	      VB(I)=HS(I-1)+1.0D0
	      VC(I)=-HS(I)-1.0D0
	      XM(I)=-2.0D0*SOURCE(I)*R(I)*R(I)*T3
	    END DO
C
C Evaluate TA,TB,TC for boudary conditions
C
	    TC(1)=0.5D0*MIDF(2)*Q(2)/DTAU(1)
	    TB(1)=-0.5D0*MIDF(1)*Q(1)/DTAU(1) - 0.5D0*HBC(ML)
	1                - GAM(1)*(HBC(ML)+NBC(ML)*SIGMA(1))
	    UC(1)=-TC(1)
	    UB(1)=0.5D0*MIDF(1)*Q(1)/DTAU(1) + 0.5D0*HBC(ML-1)
	1                - GAM(1)*(HBC(ML-1)+NBC(ML-1)*SIGMA(1))
	    XM(1)=0.0D0
	    TA(1)=0.0D0
	    UA(1)=0.0D0
	    VB(1)=0.0D0
	    VC(1)=0.0D0
C
	    TC(ND)=0.0
	    UA(ND)=0.5D0*MIDF(ND-1)*Q(ND-1)/DTAU(ND-1)
	    TA(ND)=-UA(ND)
	    IF(DIF)THEN
	      TB(ND)=0.5D0*MIDF(ND)/DTAU(ND-1)
	      UB(ND)=-TB(ND)
	      XM(ND)=DBB*R(ND)*R(ND)/3.0D0/TCHI(ND)
	    ELSE
	      TB(ND)=0.5D0*MIDF(ND)/DTAU(ND-1)+0.5D0*IN_HBC(ML)
	      UB(ND)=-TB(ND)
	      XM(ND)=R(ND)*R(ND)*IC*( 0.25D0
	1                         + 0.25D0*(IN_HBC(ML)+IN_HBC(ML-1)) )
	    END IF
	    VB(ND)=0.0D0
	    VC(ND)=0.0D0
C
	    XM(1)=XM(1) + UB(1)*JNU(1,ML-1) + UC(1)*JNU(2,ML-1)
	    DO I=2,ND-1
	      XM(I)=XM(I) + UB(I)*JNU(I,ML-1)
	1          + UA(I)*JNU(I-1,ML-1) + UC(I)*JNU(I+1,ML-1)
	1          + VB(I)*RSQHNU(I-1,ML-1) + VC(I)*RSQHNU(I,ML-1)
	    END DO
	    XM(ND)=XM(ND) + UB(ND)*JNU(ND,ML-1) + UA(ND)*JNU(ND-1,ML-1)
	  END IF
C
C
C
C Solve for the radiation field along ray for this frequency.
C
	  CALL THOMAS(TA,TB,TC,XM,ND,1)
C
	  IF(ML .EQ. 1)THEN
	    DO I=1,ND-1
	      JNU(I,ML)=XM(I)
	      RSQHNU(I,ML)=HU(I)*XM(I+1)-HL(I)*XM(I)
	    END DO
	    JNU(ND,ML)=XM(ND)
	  ELSE
	    DO I=1,ND-1
	      JNU(I,ML)=XM(I)
	      RSQHNU(I,ML)=HU(I)*(XM(I+1)+JNU(I+1,ML-1))
	1                  -HL(I)*(XM(I)+JNU(I,ML-1))
	1                   +HS(I)*RSQHNU(I,ML-1)
	    END DO
	    JNU(ND,ML)=XM(ND)
	  END IF
	END DO
C
C
C
C*****************************************************************************
C*****************************************************************************
C
C Compute JBAR.
C
	DO ML=1,NLF
	  DO I=1,ND
	    JBAR(I)=JBAR(I)+LFQW(ML)*PROF(ML)*JNU(I,ML)
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
	  CALL DP_ZERO(JEX_SCAT,ND)
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
	  XM(1)=GAM(1)*R(1)*R(1)*( (HBC(1)+SIGMA(1)*NBC(1))*JNU(1,1)
	1                - (HBC(NLF)+SIGMA(1)*NBC(NLF))*JNU(1,NLF) )
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
	1             ,HBC(NLF+1),IN_HBC(NLF+1),DIF,ND)
C
C Find the solution.
C
C After soultion, XM is integral [ jnu ] (no r^2 since used TFEAU routine.
C
	  CALL THOMAS(TA,TB,TC,XM,ND,1)
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
	  CONT_INT=13.19868*RSQHNU(1,1)
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
C intensity to be scattered directly into the definition of JNU( , NLF+1).
C
	  DO I=1,ND
	    JNU(I,NLF+1)=JNU(I,NLF+1)+JEX_SCAT(I)
	  END DO
	END IF
C
	IF(PROGDESC .NE. 236742249.0D0)THEN
	  WRITE(6,*)'Error - SCRATCH block corrupted in MOMJBAR'
	  STOP
	END IF
C
	RETURN
	END
