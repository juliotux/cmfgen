C
C Routine to compute the mean intensity J at a single frequency in the
C Comoving-Frame. The computed intensity thus depends on the intensity
C computed for the previous (bluer) frequency.
C
C The F, G, and RSQN_ON_RSQJ Eddingto factors must be supplied.
C
C NB:
C	F = K / J
C	G=  N / H
C	RSQN_ON_RSQJ(I) = RSQ_N(I)/( RSQ_J(I)+ RQS_J(I+1))
C
C where
C	RSQ_N=0.25*N(I)*(R(I)+R(I+1))**2
C
C NB: Only one of G and RSQN_ON_RSQJ is defined at a given depth. This
C     avoids having to test what mode I am using for the Eddington factors.
C
C     IF N_TYPE='G_ONLY' (in FG_J_CMF_V4) G is defined at all depths, and
C       RSQN_ON_RSQJ=0 at all depths.
C     IF N_TYPE='N_ON_J' (in FG_J_CMF_V4) RSQN_ON_RSQJ is defined at all
C       depths, and G=0 at all depths.
C     IF N_TYPE='MIXED' (in FG_J_CMF_V4) one of G or RSQN_ON_RSQJ is
C       non-zero, and is the value to be used in MOM_J_CMF
C
	SUBROUTINE MOM_J_CMF_V5(ETA,CHI,ESEC,V,SIGMA,R,
	1		   F,G,RSQN_ON_RSQJ,
	1                  F_PREV,G_PREV,RSQN_ON_RSQJ_PREV,
	1                  JNU,RSQHNU,JNU_PREV,RSQHNU_PREV,
	1                  HBC,IN_HBC,NBC,
	1                  HBC_PREV,IN_HBC_PREV,NBC_PREV,
	1                  FREQ,dLOG_NU,DIF,DBB,IC,METHOD,COHERENT,
	1                  INIT,NEW_FREQ,NC,NP,ND)
	IMPLICIT NONE
C
C Altered:  02-Jul-1998  LUER and ERROR_LU installed.
C Altered:  05-Dec-1996  PROGDESC set to number not character. Keep as REAL*8
C                          to avoid possible alignment probelms.
C Altered:  03-May-1996  Treatemnet of negative mean intensities adjusted.
C Altered:  25-Jan-1996  HBC, and HBC_PREV are now scalers. HBC(2) and
C                           HBC(3) are no longer used with the EXTENSION
C                           thick boundary condition.
C                           Changed to _V5
C Altered:   09-Mar-1995  RSQ_N_ONJ installed in effort to solve problem caused
C                           by zero H, and hence undefined G values.
C
C Finalized: 11-Nov-1994
C Created:   08-Sep-1994
C
	INTEGER NC,NP,ND,NV
	PARAMETER (NV=300)
	REAL*8 ETA(ND),CHI(ND),ESEC(ND)
	REAL*8 V(ND),SIGMA(ND),R(ND)
C
C Radiation field variables. F, G, JNU_PREV, and RSQHNU_PREV must be supplied.
C JNU and RSQHNU recomputed.
C
	REAL*8 F(ND),G(ND),RSQN_ON_RSQJ(ND)
	REAL*8 F_PREV(ND),G_PREV(ND),RSQN_ON_RSQJ_PREV(ND)
	REAL*8 JNU(ND),RSQHNU(ND)
	REAL*8 JNU_PREV(ND),RSQHNU_PREV(ND)
C
	INTEGER N_ERR_MAX,MOM_ERR_CNT
	PARAMETER (N_ERR_MAX=1000)
	REAL*8 MOM_ERR_ON_FREQ
	COMMON /MOM_J_CMF_ERR/MOM_ERR_ON_FREQ(N_ERR_MAX),MOM_ERR_CNT
	LOGICAL RECORDED_ERROR
C
C Boundary conditions.
C
	REAL*8 HBC,NBC,IN_HBC
	REAL*8 HBC_PREV,NBC_PREV,IN_HBC_PREV
C
	REAL*8 DBB,IC,FREQ,dLOG_NU
	CHARACTER*6 METHOD
C
	REAL*8, PARAMETER :: PROG_ID=1.46310043D+08  !Must be unique (MOM_J_CM)
C
C INIT is used to indicate that there is no coupling to the previous frequency.
C We are thus solving the normal continuum transfer equation (i.e. the absence
C of velocity fields)
C
C NEW_FREQ is used to indicatae that we are computing J for a new frequency.
C If we were iterating between computing J and the Eddington factors, NEW_FREQ
C would be set to false.
C
C COHERENT indicates whether the scattering is coherent. If it is, we
C explicitly take it into account. If COHERENT is FALSE, any electron
C scattering term should be included directly in the ETA that is passed
C to the routine.
C
	LOGICAL DIF,INIT,COHERENT,NEW_FREQ
C
	REAL*8 CON_GAM(NV),CON_GAMH(NV),AV_SIGMA(NV)
	SAVE CON_GAM,CON_GAMH,AV_SIGMA
C
	COMMON /SCRATCH/  PROGDESC,TA,TB,TC,DTAU,DTAUONQ,
	1                   Q,XM,SOURCE,VB,VC,HU,HL,HS,COH_VEC,GAM,GAMH,
	1                   EPS,EPS_PREV,RSQJNU_PREV,
	1                   W,WPREV,PSI,PSIPREV
	REAL*8 TA(NV),TB(NV),TC(NV),DTAU(NV),DTAUONQ(NV)
	REAL*8 Q(NV),XM(NV),SOURCE(NV)
	REAL*8 VB(NV),VC(NV),HU(NV),HL(NV),HS(NV),COH_VEC(NV)
	REAL*8 GAM(NV),GAMH(NV),W(NV),WPREV(NV)
	REAL*8 PSI(NV),PSIPREV(NV)
	REAL*8 EPS(NV),EPS_PREV(NV)
	REAL*8 RSQJNU_PREV(NV)
C
	REAL*8 PROGDESC	
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
C Local variables.
C
	INTEGER I,J
C
C PROGDESC is a variable use to confirm that the scratch block is not
C being used by some other routine.
C
	PROGDESC=PROG_ID
	IF(ND .GT. NV)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in MOM_SING_JBAR - NV smaller than ND'
	  WRITE(LUER,*)'ND=',ND,'NV',NV
	  STOP
	END IF
C
C Zero common block. There is currently 23 vectors in the common block.
C TA must be the first vector, and PSIPREV the last.
C
	PSIPREV(NV-1)=1.0D0
	PSIPREV(NV)=1.0D0
	I=(NV*23)-1
	CALL DP_ZERO(TA,I)
	IF(PSIPREV(NV-1) .NE. 0 .AND. PSIPREV(NV) .NE. 1)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in zeroing SCRATCH block in MOMJBAR'
	  STOP
	ELSE
	  PSIPREV(NV)=0.0D0
	END IF
C
C 
C
C Zero relevant vectors and matrices.
C
	CALL DP_ZERO(JNU, ND )
	CALL DP_ZERO(RSQHNU, ND )
C
	IF(INIT)THEN
	  DO I=1,N_ERR_MAX
	    MOM_ERR_ON_FREQ(I)=0.0D0
	  END DO
	  MOM_ERR_CNT=0
	END IF
C
C*****************************************************************************
C
	DO I=1,ND
	  SOURCE(I)=ETA(I)/CHI(I)
	END DO
	IF(COHERENT)THEN
	  DO I=1,ND
	    COH_VEC(I)=ESEC(I)/CHI(I)
	  END DO
	ELSE
	  DO I=1,ND
	    COH_VEC(I)=0.0D0
	  END DO
	END IF
!
	DO I=1,ND
	  IF(G(I) .GT. 1.0D0)G(I)=1.0D0
	END DO
C
C NB: We actually solve for r^2 J, not J.
C
C Compute the Q factors from F. Then compute optical depth scale.
C
	CALL QFROMF(F,Q,R,TA,TB,ND)	!TA work vector
	DO I=1,ND
	  TA(I)=CHI(I)*Q(I)
	END DO
	CALL DERIVCHI(TB,TA,R,ND,METHOD)
	CALL NORDTAU(DTAU,TA,R,R,TB,ND)
C
	IF(INIT)THEN
	  DO I=1,ND
	    GAMH(I)=0.0D0
	    GAM(I)=0.0D0
	    W(I)=0.0D0
	    WPREV(I)=0.0D0
	    PSI(I)=0.0D0
	    PSIPREV(I)=0.0D0
	    JNU_PREV(I)=0.0D0
	    RSQHNU_PREV(I)=0.0D0
	    EPS(I)=0.0D0
	    EPS_PREV(I)=0.0D0
	  END DO
C
C Assume (1)	SIGMAd+1/2 = 0.5*( SIGMAd+1+SIGMAd )
C 	 (2)	Vd+1/2=0.5*( Vd + Vd+1 )
C Note that V is in km/s and SIGMA=(dlnV/dlnR-1.0)
C
	  DO I=1,ND-1
	    CON_GAMH(I)=2.0D0*3.33564D-06*(V(I)+V(I+1))/(R(I)+R(I+1))
	    AV_SIGMA(I)=0.5D0*(SIGMA(I)+SIGMA(I+1))
	    CON_GAM(I)=3.33564D-06*V(I)/R(I)
	  END DO
	  CON_GAM(ND)=3.33564D-06*V(ND)/R(ND)
	ELSE
C
C Since we are intgerating from blue to red, FL_PREV is always larger than
C FL. dLOG_NU is define as vd / dv which is the same as d / d ln v.
C
C EPS is used if we define N in terms of J rather than H, This is sometimes
C useful as H can approach zero, and hence N/H is undefined.
C
	  DO I=1,ND-1
	    GAMH(I)=CON_GAMH(I)/dLOG_NU/( CHI(I)+CHI(I+1) )
	    W(I)=GAMH(I)*( 1.0D0+AV_SIGMA(I)*G(I) )
	    WPREV(I)=GAMH(I)*( 1.0D0+AV_SIGMA(I)*G_PREV(I) )
	    EPS(I)=GAMH(I)*AV_SIGMA(I)*RSQN_ON_RSQJ(I)/(1.0D0+W(I))
	    EPS_PREV(I)=GAMH(I)*AV_SIGMA(I)*RSQN_ON_RSQJ_PREV(I)/(1.0D0+W(I))
	  END DO
C
	  DO I=1,ND
	    GAM(I)=CON_GAM(I)/CHI(I)/dLOG_NU
	  END DO
C
C PSIPREV is equivalent to the U vector of FORMSOL.
C
	  PSI(1)=GAM(1)*(HBC+NBC*SIGMA(1))
	  PSIPREV(1)=GAM(1)*( HBC_PREV+NBC_PREV*SIGMA(1) )
	END IF
C
	DO I=2,ND
	  DTAUONQ(I)=0.5D0*(DTAU(I)+DTAU(I-1))/Q(I)
	  PSI(I)=DTAUONQ(I)*GAM(I)*( 1.0D0+SIGMA(I)*F(I) )
	  PSIPREV(I)=DTAUONQ(I)*GAM(I)*(  1.0D0+SIGMA(I)*F_PREV(I) )
	END DO
C
C NB: We are initiually computing  R^2 J. We need to multiply the
C     original JNU_PREV by R^2, since is was divided by R^2 before
C     it was stored.
C
	DO I=1,ND
	  RSQJNU_PREV(I)=R(I)*R(I)*JNU_PREV(I)
	END DO
C
C 
C
C Compute vectors used to compute the flux vector H.
C
	DO I=1,ND-1
	  HU(I)=F(I+1)*Q(I+1)/(1.0D0+W(I))/DTAU(I)
	  HL(I)=F(I)*Q(I)/(1.0D0+W(I))/DTAU(I)
	  HS(I)=WPREV(I)/(1.0D0+W(I))
	END DO
C
C Compute the TRIDIAGONAL operators, and the RHS source vector.
C
	DO I=2,ND-1
	  TA(I)=-HL(I-1)-EPS(I-1)
	  TC(I)=-HU(I)+EPS(I)
	  TB(I)=DTAUONQ(I)*(1.0D0-COH_VEC(I)) + PSI(I) + HL(I) + HU(I-1)
	1             -EPS(I-1)+EPS(I)
	  VB(I)=-HS(I-1)
	  VC(I)=HS(I)
	  XM(I)=DTAUONQ(I)*SOURCE(I)*R(I)*R(I)
	END DO
C
C Evaluate TA,TB,TC for boudary conditions
C
	TC(1)=-F(2)*Q(2)/DTAU(1)
	TB(1)=F(1)*Q(1)/DTAU(1) + PSI(1) + HBC
	XM(1)=0.0D0
	TA(1)=0.0D0
	VB(1)=0.0D0
	VC(1)=0.0D0
C
	TA(ND)=-F(ND-1)*Q(ND-1)/DTAU(ND-1)
	IF(DIF)THEN
	  TB(ND)=F(ND)/DTAU(ND-1)
	  XM(ND)=DBB*R(ND)*R(ND)/3.0D0/CHI(ND)
	ELSE
	  TB(ND)=F(ND)/DTAU(ND-1)+IN_HBC
	  XM(ND)=R(ND)*R(ND)*IC*(0.25D0+0.5D0*IN_HBC)
	END IF
	TC(ND)=0.0D0
	VB(ND)=0.0D0
	VC(ND)=0.0D0
	PSIPREV(ND)=0.0D0
C
C Note that often EPS and EPS_PREV will be zero hence we could speed this
C section up. However, MOM_J_XMF is much faster than FG_J_CMF anyway, hence
C not much tome would be gained.
C
	XM(1)=XM(1) + PSIPREV(1)*RSQJNU_PREV(1)
	DO I=2,ND-1
	  XM(I)=XM(I) + VB(I)*RSQHNU_PREV(I-1) + VC(I)*RSQHNU_PREV(I)
	1          + PSIPREV(I)*RSQJNU_PREV(I)
	1          - EPS_PREV(I-1)*(RSQJNU_PREV(I-1)+RSQJNU_PREV(I))
	1          + EPS_PREV(I)*(RSQJNU_PREV(I)+RSQJNU_PREV(I+1))
	END DO
	XM(ND)=XM(ND)
C
C Solve for the radiation field along ray for this frequency.
C
	CALL THOMAS(TA,TB,TC,XM,ND,1)
C
C Check that no negative mean intensities have been computed.
C
	IF(MINVAL(XM(1:ND)) .LE. 0)THEN
	   WRITE(47,*)'Freq=',FREQ
	   TA(1:ND)=XM(1:ND)/R(1:ND)/R(1:ND)
	   CALL WRITV(TA,ND,'XM Vec',47)
	   CALL WRITV(F,ND,'F Vec',47)
	   CALL WRITV(G,ND,'G Vec',47)
	   CALL WRITV(ETA,ND,'ETA Vec',47)
	   CALL WRITV(ESEC,ND,'ESEC Vec',47)
	   CALL WRITV(CHI,ND,'CHI Vec',47)
	END IF
!
	DO I=1,ND
	  IF(XM(I) .LT. 0)THEN
	    XM(I)=ABS(XM(I))/10.0D0
	    RECORDED_ERROR=.FALSE.
	    J=1
	    DO WHILE (J .LE. MOM_ERR_CNT .AND. .NOT. RECORDED_ERROR)
	      IF(MOM_ERR_ON_FREQ(J) .EQ. FREQ)RECORDED_ERROR=.TRUE.
	      J=J+1
	    END DO
	    IF(.NOT. RECORDED_ERROR .AND. MOM_ERR_CNT .LT. N_ERR_MAX)THEN
	      MOM_ERR_CNT=MOM_ERR_CNT+1
	      MOM_ERR_ON_FREQ(MOM_ERR_CNT)=FREQ
	    END IF	
	  END IF
	END DO
C
C Store J, correcting for the fact that we actually compute r^2 J
C
	DO I=1,ND
	  JNU(I)=XM(I)/R(I)/R(I)
	END DO
C
	DO I=1,ND-1
	  RSQHNU(I)=HU(I)*XM(I+1)-HL(I)*XM(I)+HS(I)*RSQHNU_PREV(I)+
	1              ( EPS_PREV(I)*(RSQJNU_PREV(I)+RSQJNU_PREV(I+1)) -
	1                  EPS(I)*(XM(I)+XM(I+1)) )
	END DO
C
	IF(PROGDESC .NE. PROG_ID)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error - SCRATCH block corrupted in MOM_J_CMF'
	  STOP
	END IF
C
	RETURN
	END
