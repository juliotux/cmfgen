
! Data module for PP_MOM_CMF_V1. Data placed in this module is automatically
! saved between subroutine calls..
!
	MODULE PP_MOM_CMF_MOD
	IMPLICIT NONE
!
! To be dimenensioned ND_SM where ND_SM is the size of the R grid
! as passed to MOM_J_CMF.
!
	REAL*8, ALLOCATABLE :: LOG_R_SM(:)
!
! Dimensioned ND_SM,4
!
	REAL*8, ALLOCATABLE :: V_COEF(:,:)
	REAL*8, ALLOCATABLE :: SIGMA_COEF(:,:)
	REAL*8, ALLOCATABLE :: ESEC_COEF(:,:)
	REAL*8, ALLOCATABLE :: CHI_COEF(:,:)
	REAL*8, ALLOCATABLE :: ETA_COEF(:,:)
	REAL*8, ALLOCATABLE :: F_COEF(:,:)
	REAL*8, ALLOCATABLE :: G_COEF(:,:)
	REAL*8, ALLOCATABLE :: N_ON_J_COEF(:,:)
!
! All the following rays have dimension ND, where ND >= ND_SM.
! Some of the data in the arrays is need in subsequent calls.
! 
	REAL*8, ALLOCATABLE :: R(:)
	REAL*8, ALLOCATABLE :: V(:)
	REAL*8, ALLOCATABLE :: SIGMA(:)
	REAL*8, ALLOCATABLE :: ESEC(:)
	REAL*8, ALLOCATABLE :: CHI(:)
	REAL*8, ALLOCATABLE :: ETA(:)
!
	REAL*8, ALLOCATABLE :: JNU(:)
	REAL*8, ALLOCATABLE :: HNU(:)
	REAL*8, ALLOCATABLE :: F(:)
	REAL*8, ALLOCATABLE :: G(:)
	REAL*8, ALLOCATABLE :: N_ON_J(:)
	REAL*8, ALLOCATABLE :: F_SAV(:)
	REAL*8, ALLOCATABLE :: G_SAV(:)
	REAL*8, ALLOCATABLE :: N_ON_J_SAV(:)
!
	REAL*8, ALLOCATABLE :: JNU_PREV(:)
	REAL*8, ALLOCATABLE :: HNU_PREV(:)
	REAL*8, ALLOCATABLE :: F_PREV(:)
	REAL*8, ALLOCATABLE :: G_PREV(:)
	REAL*8, ALLOCATABLE :: N_ON_J_PREV(:)
!
	REAL*8, ALLOCATABLE :: CON_GAM(:)
	REAL*8, ALLOCATABLE :: CON_GAMH(:)
	REAL*8, ALLOCATABLE :: AV_SIGMA(:)
!
	REAL*8, ALLOCATABLE :: DTAU(:)
	REAL*8, ALLOCATABLE :: DTAU_MID(:)
	REAL*8, ALLOCATABLE :: TA(:)
	REAL*8, ALLOCATABLE :: TB(:)
	REAL*8, ALLOCATABLE :: TC(:)
	REAL*8, ALLOCATABLE :: XM(:)
	REAL*8, ALLOCATABLE :: SOURCE(:)
	REAL*8, ALLOCATABLE :: VB(:)
	REAL*8, ALLOCATABLE :: VC(:)
	REAL*8, ALLOCATABLE :: HU(:)
	REAL*8, ALLOCATABLE :: HL(:)
	REAL*8, ALLOCATABLE :: HS(:)
	REAL*8, ALLOCATABLE :: COH_VEC(:)
	REAL*8, ALLOCATABLE :: GAM(:)
	REAL*8, ALLOCATABLE :: GAMH(:)
	REAL*8, ALLOCATABLE :: W(:)
	REAL*8, ALLOCATABLE :: WPREV(:)
	REAL*8, ALLOCATABLE :: PSI(:)
	REAL*8, ALLOCATABLE :: PSIPREV(:)
	REAL*8, ALLOCATABLE :: EPS(:)
	REAL*8, ALLOCATABLE :: EPS_PREV(:)
!
! Boundary conditions
!
        REAL*8 HBC_PREV
	REAL*8 IN_HBC_PREV
	REAL*8 NBC_PREV
	REAL*8 NBC_INCID_PREV
!
	REAL*8 HBC_SAV
	REAL*8 IN_HBC_SAV
	REAL*8 NBC_SAV
	REAL*8 NBC_INCID_SAV
	REAL*8 VDOP_FRAC_SAV
!
	INTEGER ND
!
! R_PNT(K) defines the interpolation for the variable at depth K.
!
	INTEGER, ALLOCATABLE :: R_PNT(:)
!
! ?_INDX are used to indicate the location of J and H on the small grid
! in the larger array.
!
	INTEGER, ALLOCATABLE :: J_INDX(:)
	INTEGER, ALLOCATABLE :: H_INDX(:)
!
	LOGICAL FIRST_TIME
	DATA FIRST_TIME/.TRUE./
	DATA VDOP_FRAC_SAV/-10001.1D0/    !Absurd value
!
	END MODULE PP_MOM_CMF_MOD
!
!
!
!
! Routine to compute the mean intensity J at a single frequency in the
! Comoving-Frame. The computed intensity thus depends on the intensity
! computed for the previous (bluer) frequency. Routine is designed for
! a plane-parallel atomoshere.
!
! The F, G, and N_ON_J Eddingto factors must be supplied.
!
! NB:
!	F = K / J
!	G=  N / H
!	N_ON_J(I) = N(I)/(J(I)+J(I+1))
!
! NB: Only one of G and N_ON_J is defined at a given depth. This
!     avoids having to test what mode I am using for the Eddington factors.
!
!     IF N_TYPE='G_ONLY' (in FG_J_CMF_V4) G is defined at all depths, and
!       N_ON_J=0 at all depths.
!     IF N_TYPE='N_ON_J' (in FG_J_CMF_V4) N_ON_J is defined at all
!       depths, and G=0 at all depths.
!     IF N_TYPE='MIXED' (in FG_J_CMF_V4) one of G or N_ON_J is
!       non-zero, and is the value to be used in MOM_J_CMF
!
	SUBROUTINE PP_MOM_CMF_V1(ETA_SM,CHI_SM,ESEC_SM,
	1                  V_SM,SIGMA_SM,R_SM,
	1		   F_SM,G_SM,N_ON_J_SM,
	1                  JNU_SM,HNU_SM,
	1                  VDOP_VEC,VDOP_FRAC,
	1                  IN_HBC,HBC,HBC_INCID,NBC,NBC_INCID,
	1                  FREQ,dLOG_NU,DIF,DBB,IC,
	1                  N_TYPE,METHOD,COHERENT,
	1                  INIT,NEW_FREQ,ND_SM)
	USE PP_MOM_CMF_MOD
	IMPLICIT NONE
!
! Altered:   17-Feb-2007 : dNBC_INCID was not zeroed on when INIT set.
! Created:   02-Mar-2006 : Based on MOM_J_CMF_V6
!
	INTEGER ND_SM
	REAL*8 ETA_SM(ND_SM)
	REAL*8 CHI_SM(ND_SM)
	REAL*8 ESEC_SM(ND_SM)
	REAL*8 V_SM(ND_SM)
	REAL*8 SIGMA_SM(ND_SM)
	REAL*8 R_SM(ND_SM)
!
! Radiation field variables. F, G, JNU_PREV, and HNU_PREV must be supplied.
! JNU and HNU recomputed.
!
	REAL*8 F_SM(ND_SM)
	REAL*8 G_SM(ND_SM)
	REAL*8 N_ON_J_SM(ND_SM)
	REAL*8 JNU_SM(ND_SM)
	REAL*8 HNU_SM(ND_SM)
	REAL*8 VDOP_VEC(ND_SM)
	REAL*8 VDOP_FRAC
!
	INTEGER N_ERR_MAX,MOM_ERR_CNT
	PARAMETER (N_ERR_MAX=1000)
	REAL*8 MOM_ERR_ON_FREQ
	COMMON /MOM_J_CMF_ERR/MOM_ERR_ON_FREQ(N_ERR_MAX),MOM_ERR_CNT
	LOGICAL RECORDED_ERROR
!
! Boundary conditions.
!
	REAL*8 HBC,HBC_INCID
	REAL*8 NBC,NBC_INCID
	REAL*8 IN_HBC
!
	REAL*8 DBB,IC,FREQ,dLOG_NU
	CHARACTER*6 METHOD
	CHARACTER*6 N_TYPE 
!
! INIT is used to indicate that there is no coupling to the previous frequency.
! We are thus solving the normal continuum transfer equation (i.e. the absence
! of velocity fields)
!
! NEW_FREQ is used to indicatae that we are computing J for a new frequency.
! If we were iterating between computing J and the Eddington factors, NEW_FREQ
! would be set to false.
!
! COHERENT indicates whether the scattering is coherent. If it is, we
! explicitly take it into account. If COHERENT is FALSE, any electron
! scattering term should be included directly in the ETA that is passed
! to the routine.
!
	LOGICAL DIF,INIT,COHERENT,NEW_FREQ
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
! Local variables.
!
	REAL*8 T1,T2
	REAL*8  DELTA_R
	REAL*8  dNBC_INCID
	INTEGER IT1,I,J,K
	INTEGER IOS
	LOGICAL NEW_R_GRID
!
!
!
	NEW_R_GRID=.FALSE.
	IF(INIT .AND. ALLOCATED(R))THEN
	  DO I=1,ND_SM
	    IF(LOG(R_SM(I)) .NE. LOG_R_SM(I))THEN
	      NEW_R_GRID=.TRUE.
	      WRITE(171,*)'Updating RGRID in MOM_J_CMF_V6'
	      EXIT
	    END IF
	  END DO
          IF(VDOP_FRAC .NE. VDOP_FRAC_SAV)NEW_R_GRID=.TRUE.
	END IF
	IF(FIRST_TIME)THEN
          OPEN(UNIT=47,FILE='MOM_J_ERRORS',STATUS='UNKNOWN')
        ELSE IF(INIT)THEN
	  REWIND(47)
	END IF
!
! Deallocate all arrayes if we have changed VDOP_FRAC. This will only
! be done in testing this routine (e.g., using DISPGEN).
!
	IF(ALLOCATED(R) .AND. NEW_R_GRID)THEN
	  DEALLOCATE ( R, R_PNT, LOG_R_SM )
	  DEALLOCATE ( V_COEF, SIGMA_COEF, ETA_COEF, ESEC_COEF, CHI_COEF )
	  DEALLOCATE ( F_COEF, G_COEF, N_ON_J_COEF )
	  DEALLOCATE ( V, SIGMA, ETA, ESEC, CHI )
	  DEALLOCATE ( JNU, HNU, F, G, N_ON_J )
	  DEALLOCATE ( F_SAV, G_SAV, N_ON_J_SAV )
	  DEALLOCATE ( JNU_PREV, HNU_PREV, F_PREV, G_PREV, N_ON_J_PREV )
	  DEALLOCATE ( CON_GAM, CON_GAMH, AV_SIGMA )
	  DEALLOCATE ( DTAU, DTAU_MID, TA, TB, TC, XM, SOURCE )
	  DEALLOCATE ( VB, VC, HU, HL, HS, COH_VEC )
	  DEALLOCATE ( GAM, GAMH, W, WPREV, PSI, PSIPREV )
	  DEALLOCATE ( EPS, EPS_PREV, J_INDX, H_INDX )
	END IF
	VDOP_FRAC_SAV=VDOP_FRAC
!
! On the very first entry, we define the improved R grid, and allocate all
! data arrays.
!
	IF( FIRST_TIME .OR. .NOT. ALLOCATED(R) )THEN
!
! Determine the number of points for the expanded R grid.
! We always insert an EVEN number of points. This guarentees that
! H_SM (defined at the midpoints of the pass grid) has an exact correspondence
! with H defined on the extended gid.
!
          K=1
          T2=VDOP_FRAC*MINVAL(VDOP_VEC(1:ND_SM))
          DO I=1,ND_SM-1
            IT1=INT( (V_SM(I)-V_SM(I+1))/T2 )
	    IF( MOD(IT1,2) .NE. 0)IT1=IT1+1
            IF(IT1 .GT. 0)K=K+IT1
            K=K+1
          END DO
          ND=K
!
	  ALLOCATE ( R(ND), STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( R_PNT(ND), STAT=IOS )
	  IF(IOS .NE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Unable to allocate R & R_PNT in PP_MOM_CMF_V1'
	    WRITE(LUER,*)'Status=',IOS
	  END IF  
! 
          K=1
	  R(1)=R_SM(1)
          R_PNT(1)=1
	  DO I=1,ND_SM-1
            IT1=INT( (V_SM(I)-V_SM(I+1))/T2 )
	    IF( MOD(IT1,2) .NE. 0)IT1=IT1+1
            IF(IT1 .GT. 0)THEN
              DELTA_R=(R_SM(I+1)-R_SM(I))/(IT1+1)
              DO J=1,IT1
                K=K+1
                R(K)=R(K-1)+DELTA_R
                R_PNT(K)=I
              END DO
            END IF
            K=K+1
            R(K)=R_SM(I+1)
            R_PNT(K)=I
	  END DO
!
	  ALLOCATE ( LOG_R_SM(ND_SM),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (V_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (SIGMA_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ETA_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ESEC_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CHI_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (F_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (G_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (N_ON_J_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	     LUER=ERROR_LU()
	     WRITE(lUER,*)'Unable to allocate COEF memory in PP_MOM_CMF_V1'
	     STOP
	  END IF
!
	  ALLOCATE ( V(ND), SIGMA(ND), ETA(ND), ESEC(ND), CHI(ND), STAT=IOS) 
	  IF(IOS .EQ. 0)ALLOCATE (JNU(ND), HNU(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (F(ND), G(ND), N_ON_J(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (F_SAV(ND), G_SAV(ND), N_ON_J_SAV(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (JNU_PREV(ND), HNU_PREV(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (F_PREV(ND), G_PREV(ND), N_ON_J_PREV(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CON_GAM(ND), CON_GAMH(ND), AV_SIGMA(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (DTAU(ND), DTAU_MID(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (TA(ND),  TB(ND),  TC(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (XM(ND), SOURCE(ND), COH_VEC(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (VB(ND), VC(ND), HU(ND), HL(ND), HS(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (GAM(ND), GAMH(ND), W(ND), WPREV(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (PSI(ND), PSIPREV(ND), EPS(ND), EPS_PREV(ND), STAT=IOS)
	  IF(IOS .NE. 0)THEN
	     LUER=ERROR_LU()
	     WRITE(lUER,*)'Unable to allocate vectors in PP_MOM_CMF_V1'
	     STOP
	  END IF
!
!
	  LOG_R_SM(1:ND_SM)=LOG(R_SM(1:ND_SM))
          CALL MON_INT_FUNS_V2(V_COEF,V_SM,LOG_R_SM,ND_SM)
          CALL MON_INT_FUNS_V2(SIGMA_COEF,SIGMA_SM,LOG_R_SM,ND_SM)
          DO I=1,ND
            K=R_PNT(I)
            T1=LOG(R(I)/R_SM(K))
            V(I)=((V_COEF(K,1)*T1+V_COEF(K,2))*T1+V_COEF(K,3))*T1+V_COEF(K,4)
            SIGMA(I)=((SIGMA_COEF(K,1)*T1+SIGMA_COEF(K,2))*T1+SIGMA_COEF(K,3))*T1+SIGMA_COEF(K,4)
          END DO
!
	  ALLOCATE ( J_INDX(ND_SM),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( H_INDX(ND_SM),STAT=IOS )
	  IF(IOS .NE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Unable to allocate J_INDX & H_INDX in PP_MOM_CMF_V1'
	    WRITE(LUER,*)'Status=',IOS
	  END IF
	  J_INDX(1:ND_SM)=0;    H_INDX(1:ND_SM)=0
!  
	  K=1
	  DO I=1,ND_SM
	    DO WHILE(J_INDX(I) .EQ. 0)
	      IF(R_SM(I) .LE. R(K) .AND. R_SM(I) .GE. R(K+1))THEN
	        IF( (R(K)-R_SM(I)) .LT. (R_SM(I)-R(K+1)) )THEN
	          J_INDX(I)=K
	        ELSE
	          J_INDX(I)=K+1
	        END IF
	      ELSE
	        K=K+1
	      END IF
	    END DO
	  END DO
!
	  K=1
	  DO I=1,ND_SM-1
	    T1=0.5D0*(R_SM(I)+R_SM(I+1))
	    DO WHILE(H_INDX(I) .EQ. 0)
	      IF(T1 .LT. R(K) .AND. T1 .GT. R(K+1))THEN
	        H_INDX(I)=K
	      ELSE
	        K=K+1
	      END IF
	    END DO
	  END DO
!
	  IF(N_TYPE .NE. 'G_ONLY' .AND. N_TYPE .NE. 'N_ON_J')THEN
	    IF(ND .NE. ND_SM)THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Error in MOM_J_CMF_V6'
	      WRITE(LUER,*)'Cannot use N_TYPE=''MIXED'' when inserting extra points'
	    END IF
	  END IF
!
	  FIRST_TIME=.FALSE.
	END IF
!
!
	IF(ND .GT. ND_SM)THEN
!
! Interpolate quantities onto revised grid.
!
	  TA(1:ND_SM)=LOG(CHI_SM(1:ND_SM))
	  CALL MON_INT_FUNS_V2(CHI_COEF,TA,LOG_R_SM,ND_SM)
	  TA(1:ND_SM)=LOG(ESEC_SM(1:ND_SM))
	  CALL MON_INT_FUNS_V2(ESEC_COEF,TA,LOG_R_SM,ND_SM)
	  TA(1:ND_SM)=LOG(ETA_SM(1:ND_SM))
	  CALL MON_INT_FUNS_V2(ETA_COEF,TA,LOG_R_SM,ND_SM)
	  CALL MON_INT_FUNS_V2(F_COEF,F_SM,LOG_R_SM,ND_SM)
	  CALL MON_INT_FUNS_V2(G_COEF,G_SM,LOG_R_SM,ND_SM)
	  CALL MON_INT_FUNS_V2(N_ON_J_COEF,N_ON_J_SM,LOG_R_SM,ND_SM)

	  DO I=1,ND
	    K=R_PNT(I)
	    T1=LOG(R(I)/R_SM(K))
	    T2=((CHI_COEF(K,1)*T1+CHI_COEF(K,2))*T1+CHI_COEF(K,3))*T1+CHI_COEF(K,4)
	    CHI(I)=EXP(T2)
	    T2=((ESEC_COEF(K,1)*T1+ESEC_COEF(K,2))*T1+ESEC_COEF(K,3))*T1+ESEC_COEF(K,4)
	    ESEC(I)=EXP(T2)
	    T2=((ETA_COEF(K,1)*T1+ETA_COEF(K,2))*T1+ETA_COEF(K,3))*T1+ETA_COEF(K,4)
	    ETA(I)=EXP(T2)
	    F(I)=((F_COEF(K,1)*T1+F_COEF(K,2))*T1+F_COEF(K,3))*T1+F_COEF(K,4)
	    G(I)=((G_COEF(K,1)*T1+G_COEF(K,2))*T1+G_COEF(K,3))*T1+G_COEF(K,4)
	    N_ON_J(I)=((N_ON_J_COEF(K,1)*T1+N_ON_J_COEF(K,2))*T1+
	1                               N_ON_J_COEF(K,3))*T1+N_ON_J_COEF(K,4)
	  END DO
!
	ELSE
	  ESEC(1:ND)=ESEC_SM(1:ND)
	  CHI(1:ND)=CHI_SM(1:ND)
	  ETA(1:ND)=ETA_SM(1:ND)
	  F(1:ND)=F_SM(1:ND)
	  G(1:ND)=G_SM(1:ND)
	  N_ON_J(1:ND)=N_ON_J_SM(1:ND)
	END IF
!
! 
!

	IF(INIT)THEN
	  DO I=1,N_ERR_MAX
	    MOM_ERR_ON_FREQ(I)=0.0D0
	  END DO
	  MOM_ERR_CNT=0
          JNU=0.0D0; HNU=0.0D0
	  JNU_PREV=0.0D0; HNU_PREV=0.0D0
	  N_ON_J_SAV=0.0D0; G_SAV=0.0D0; F_SAV=0.0D0
	END IF
!
!*****************************************************************************
!
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
	  IF(G(I) .LT. -1.0D0)G(I)=-1.0D0
	END DO
!
! NB: We solve for J.
!
!
	CALL DERIVCHI(TB,CHI,R,ND,METHOD)
	CALL NORDTAU(DTAU,CHI,R,R,TB,ND)
!
	IF(NEW_FREQ)THEN
	  F_PREV(1:ND)=F_SAV(1:ND)
	  G_PREV(1:ND)=G_SAV(1:ND)
	  N_ON_J_PREV(1:ND)=N_ON_J_SAV(1:ND)
	  JNU_PREV(1:ND)=JNU(1:ND)
	  HNU_PREV(1:ND)=HNU(1:ND)
	  HBC_PREV=HBC_SAV;  IN_HBC_PREV=IN_HBC_SAV
	  NBC_PREV=NBC_SAV
	  NBC_INCID_PREV=NBC_INCID_SAV
	END IF
!
!
	IF(INIT)THEN
	  DO I=1,ND
	    GAMH(I)=0.0D0
	    GAM(I)=0.0D0
	    W(I)=0.0D0
	    WPREV(I)=0.0D0
	    PSI(I)=0.0D0
	    PSIPREV(I)=0.0D0
	    JNU_PREV(I)=0.0D0
	    HNU_PREV(I)=0.0D0
	    EPS(I)=0.0D0
	    EPS_PREV(I)=0.0D0
	  END DO
	  HBC_PREV=0.0D0;  IN_HBC_PREV=0.0D0; NBC_PREV=0.0D0
	  dNBC_INCID=0.0D0
	ELSE
!
! Assume (1)	SIGMAd+1/2 = 0.5*( SIGMAd+1+SIGMAd )
! 	 (2)	Vd+1/2=0.5*( Vd + Vd+1 )
! Note that V is in km/s and SIGMA=(dlnV/dlnR-1.0)
!
	  DO I=1,ND-1
	    CON_GAMH(I)=2.0D0*3.33564D-06*(V(I)+V(I+1))/(R(I)+R(I+1))
	    AV_SIGMA(I)=0.5D0*(SIGMA(I)+SIGMA(I+1))
	    CON_GAM(I)=3.33564D-06*V(I)/R(I)
	  END DO
	  CON_GAM(ND)=3.33564D-06*V(ND)/R(ND)
!
! Since we are intgerating from blue to red, FL_PREV is always larger than
! FL. dLOG_NU is define as vd / dv which is the same as d / d ln v.
!
! EPS is used if we define N in terms of J rather than H, This is sometimes
! useful as H can approach zero, and hence N/H is undefined.
!
	  IF(N_TYPE .EQ. 'G_ONLY')THEN
	    DO I=1,ND-1
	      GAMH(I)=CON_GAMH(I)/dLOG_NU/( CHI(I)+CHI(I+1) )
	      W(I)=GAMH(I)*(1.0D0+AV_SIGMA(I))*G(I)
	      WPREV(I)=GAMH(I)*(1.0D0+AV_SIGMA(I))*G_PREV(I) 
	    END DO
	  ELSE
	    DO I=1,ND-1
	      GAMH(I)=CON_GAMH(I)/dLOG_NU/( CHI(I)+CHI(I+1) )
	      W(I)=GAMH(I)*(1.0D0+AV_SIGMA(I))*G(I)
	      WPREV(I)=GAMH(I)*(1.0D0+AV_SIGMA(I))*G_PREV(I)
	      EPS(I)=GAMH(I)*(1.0D0+AV_SIGMA(I))*N_ON_J(I)/(1.0D0+W(I))
	      EPS_PREV(I)=GAMH(I)*(1.0D0+AV_SIGMA(I))*N_ON_J_PREV(I)/(1.0D0+W(I))
	    END DO
	  END IF
!
	  DO I=1,ND
	    GAM(I)=CON_GAM(I)/CHI(I)/dLOG_NU
	  END DO
!
! PSIPREV is equivalent to the U vector of FORMSOL.
!
	  PSI(1)=GAM(1)*NBC*(1.0D0+SIGMA(1))
	  PSIPREV(1)=GAM(1)*NBC_PREV*(1.0D0+SIGMA(1))
	  dNBC_INCID=GAM(1)*(1.0D0+SIGMA(1))*(NBC_INCID-NBC_INCID_PREV)
	END IF
!
! 
!
	DO I=2,ND-1
	  DTAU_MID(I)=0.5D0*(DTAU(I)+DTAU(I-1))
	  PSI(I)=DTAU_MID(I)*GAM(I)*(1.0D0+SIGMA(I))*F(I)
	  PSIPREV(I)=DTAU_MID(I)*GAM(I)*(1.0D0+SIGMA(I))*F_PREV(I)
	END DO
!
! Compute vectors used to compute the flux vector H.
!
	DO I=1,ND-1
	  HU(I)=F(I+1)/(1.0D0+W(I))/DTAU(I)
	  HL(I)=F(I)/(1.0D0+W(I))/DTAU(I)
	  HS(I)=WPREV(I)/(1.0D0+W(I))
	END DO
!
! Compute the TRIDIAGONAL operators, and the RHS source vector.
!
	IF(N_TYPE .EQ. 'G_ONLY')THEN
	  DO I=2,ND-1
	    TA(I)=-HL(I-1)
	    TC(I)=-HU(I)
	    TB(I)=DTAU_MID(I)*(1.0D0-COH_VEC(I)) + PSI(I) + HL(I) + HU(I-1)
	    VB(I)=-HS(I-1)
	    VC(I)=HS(I)
	    XM(I)=DTAU_MID(I)*SOURCE(I)
	  END DO
	ELSE
	  DO I=2,ND-1
	    TA(I)=-HL(I-1)-EPS(I-1)
	    TC(I)=-HU(I)+EPS(I)
	    TB(I)=DTAU_MID(I)*(1.0D0-COH_VEC(I)) + PSI(I) + HL(I) + HU(I-1)
	1               -EPS(I-1)+EPS(I)
	    VB(I)=-HS(I-1)
	    VC(I)=HS(I)
	    XM(I)=DTAU_MID(I)*SOURCE(I)
	  END DO
	END IF
!
! Evaluate TA,TB,TC for boudary conditions
!
	TC(1)=-F(2)/DTAU(1)
	TB(1)=F(1)/DTAU(1) + PSI(1) + HBC
	XM(1)=HBC_INCID+dNBC_INCID
	TA(1)=0.0D0
	VB(1)=0.0D0
	VC(1)=0.0D0
!
	TA(ND)=-F(ND-1)/DTAU(ND-1)
	IF(DIF)THEN
	  TB(ND)=F(ND)/DTAU(ND-1)
	  XM(ND)=DBB/3.0D0/CHI(ND)
	ELSE
	  TB(ND)=F(ND)/DTAU(ND-1)+IN_HBC
	  XM(ND)=IC*(0.25D0+0.5D0*IN_HBC)
	END IF
	TC(ND)=0.0D0
	VB(ND)=0.0D0
	VC(ND)=0.0D0
	PSIPREV(ND)=0.0D0
!
! Note that EPS and EPS_PREV will be identically zero hence when N_TYPE is
! G_ONLY.
!
	IF(N_TYPE .EQ. 'G_ONLY')THEN
	  XM(1)=XM(1) + PSIPREV(1)*JNU_PREV(1)
	  DO I=2,ND-1
	    XM(I)=XM(I) + VB(I)*HNU_PREV(I-1) + VC(I)*HNU_PREV(I)
	1          + PSIPREV(I)*JNU_PREV(I)
	  END DO
	  XM(ND)=XM(ND)
	ELSE
	  XM(1)=XM(1) + PSIPREV(1)*JNU_PREV(1)
	  DO I=2,ND-1
	    XM(I)=XM(I) + VB(I)*HNU_PREV(I-1) + VC(I)*HNU_PREV(I)
	1          + PSIPREV(I)*JNU_PREV(I)
	1          - EPS_PREV(I-1)*(JNU_PREV(I-1)+JNU_PREV(I))
	1          + EPS_PREV(I)*(JNU_PREV(I)+JNU_PREV(I+1))
	  END DO
	  XM(ND)=XM(ND)
	END IF
!
! Solve for the radiation field along ray for this frequency.
!
	CALL THOMAS(TA,TB,TC,XM,ND,1)
!
! Check that no negative mean intensities have been computed.
!
	IF(MINVAL(XM(1:ND)) .LE. 0.0D0)THEN
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
	RECORDED_ERROR=.FALSE.
	DO I=1,ND
	  IF(XM(I) .LT. 0.0D0)THEN
	    XM(I)=ABS(XM(I))/10.0D0
	  END IF
	  IF(.NOT. RECORDED_ERROR)THEN
	    IF(MOM_ERR_CNT .GT. N_ERR_MAX)THEN
	      MOM_ERR_CNT=MOM_ERR_CNT+1
	    ELSE IF(MOM_ERR_ON_FREQ(MOM_ERR_CNT) .NE. FREQ)THEN
	      MOM_ERR_CNT=MOM_ERR_CNT+1
	      IF(MOM_ERR_CNT .LT. N_ERR_MAX)MOM_ERR_ON_FREQ(MOM_ERR_CNT)=FREQ
	    END IF
	    RECORDED_ERROR=.TRUE.
	  END IF
	END DO
!
! Save R^2 J for next iteration.
!
	JNU(1:ND)=XM(1:ND)
!
	IF(N_TYPE .EQ. 'G_ONLY')THEN
	  DO I=1,ND-1
	    HNU(I)=HU(I)*XM(I+1)-HL(I)*XM(I)+HS(I)*HNU_PREV(I)
	  END DO
	ELSE
	  DO I=1,ND-1
	    HNU(I)=HU(I)*XM(I+1)-HL(I)*XM(I)+HS(I)*HNU_PREV(I)+
	1              ( EPS_PREV(I)*(JNU_PREV(I)+JNU_PREV(I+1)) -
	1                  EPS(I)*(XM(I)+XM(I+1)) )
	  END DO
	END IF
!
! Regrid derived J and H values onto small grid. We devide J by R^2 so that
! we return J.
!
	DO I=1,ND_SM
	  K=J_INDX(I)
	  JNU_SM(I)=JNU(K)
	END DO
!
	DO I=1,ND_SM-1
	  K=H_INDX(I)
	  HNU_SM(I)=HNU(K)
	END DO
!
	F_SAV(1:ND)=F(1:ND)
	G_SAV(1:ND)=G(1:ND)
	N_ON_J_SAV(1:ND)=N_ON_J(1:ND)
	HBC_SAV=HBC
	IN_HBC_SAV=IN_HBC
	NBC_SAV=NBC
	NBC_INCID_SAV=NBC_INCID
!
	RETURN
	END
