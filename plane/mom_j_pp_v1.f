
! Data module for MOM_PP_V1. Data placed in this module is automatically
! saved between subroutine calls.
!
	MODULE MOD_MOM_PP_J_V1
	IMPLICIT NONE
!
! To be dimenensioned ND_SM where ND_SM is the size of the R grid
! as passed to MOM_J_CMF.
!
	REAL*8, ALLOCATABLE :: LOG_R_SM(:)
!
! Dimensioned ND_SM,4
!
	REAL*8, ALLOCATABLE :: ESEC_COEF(:,:)
	REAL*8, ALLOCATABLE :: CHI_COEF(:,:)
	REAL*8, ALLOCATABLE :: ETA_COEF(:,:)
	REAL*8, ALLOCATABLE :: F_COEF(:,:)
!
! All the following rays have dimension ND, where ND >= ND_SM.
! Some of the data in the arrays is need in subsequent calls.
! 
	REAL*8, ALLOCATABLE :: R(:)
	REAL*8, ALLOCATABLE :: ESEC(:)
	REAL*8, ALLOCATABLE :: CHI(:)
	REAL*8, ALLOCATABLE :: ETA(:)
!
	REAL*8, ALLOCATABLE :: JNU(:)
	REAL*8, ALLOCATABLE :: HNU(:)
	REAL*8, ALLOCATABLE :: F(:)
!
	REAL*8, ALLOCATABLE :: DTAU(:)
	REAL*8, ALLOCATABLE :: MID_DTAU(:)
	REAL*8, ALLOCATABLE :: TA(:)
	REAL*8, ALLOCATABLE :: DD(:)		!TB=-DD-TA-TC
	REAL*8, ALLOCATABLE :: TC(:)
	REAL*8, ALLOCATABLE :: XM(:)
	REAL*8, ALLOCATABLE :: SOURCE(:)
	REAL*8, ALLOCATABLE :: HU(:)
	REAL*8, ALLOCATABLE :: HL(:)
	REAL*8, ALLOCATABLE :: COH_VEC(:)
!
	INTEGER ND
	INTEGER, SAVE :: NINS_SAV=0
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
!
	END MODULE MOD_MOM_PP_J_V1
!
!
!
! Routine to compute the mean intensity J at a single frequency for
! a plane-parallel atmosphere in the observer's frame.
!
! The F Eddington factors must be supplied.
!
! NB:
!	F = K / J
!
	SUBROUTINE MOM_J_PP_V1(ETA_SM,CHI_SM,ESEC_SM,
	1                  R_SM,F_SM,JNU_SM,HNU_SM,
	1                  HBC_J,HBC_S,IN_HBC,
	1                  FREQ,DIF,DBB,IC,METHOD,COHERENT,
	1                  PASSED_NINS,INIT,NEW_FREQ,ND_SM)
	USE MOD_MOM_PP_J_V1
	IMPLICIT NONE
!
! Altered 3-Feb-2008 : Changed to allow for a variable R grid,
!
	INTEGER ND_SM
	REAL*8 ETA_SM(ND_SM)
	REAL*8 CHI_SM(ND_SM)
	REAL*8 ESEC_SM(ND_SM)
	REAL*8 R_SM(ND_SM)
!
! Radiation field variables JNU and HNU are recomputed.
!
	REAL*8 F_SM(ND_SM)
	REAL*8 JNU_SM(ND_SM)
	REAL*8 HNU_SM(ND_SM)
!
	INTEGER N_ERR_MAX,MOM_ERR_CNT
	PARAMETER (N_ERR_MAX=1000)
	REAL*8 MOM_ERR_ON_FREQ
	COMMON /MOM_J_CMF_ERR/MOM_ERR_ON_FREQ(N_ERR_MAX),MOM_ERR_CNT
	LOGICAL RECORDED_ERROR
!
	INTEGER PASSED_NINS
	INTEGER NINS
!
! Boundary conditions.
!
	REAL*8 HBC_J,HBC_S,IN_HBC
!
	REAL*8 DBB,IC,FREQ
	CHARACTER*6 METHOD
!
! INIT is used to indicate that this is the first frequency. At present, this
! is only used to initialize the error arrays.
!
! NEW_FREQ is used to indicate that we are computing J for a new frequency.
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
	INTEGER I,J,K
	INTEGER IOS
!
!
!
! Determine the number of points for the expanded R grid. With out
! a velocity field, this is imply related to NINS and ND.
! We always insert an EVEN number of points. This guarentees that
! H_SM (defined at the midpoints of the pass grid) has an exact correspondence
! with H defined on the extended gid.
!
	NINS=PASSED_NINS
	IF(MOD(NINS,2) .NE. 0)NINS=NINS+1
	ND=(ND_SM-1)*(NINS+1)+1
!
! Deallocate all arrayes if we have changed VDOP_FRAC. This will only
! be done in testing this routine (e.g., using DISPGEN).
!
	IF(ALLOCATED(R) .AND. NINS .NE. NINS_SAV)THEN
	  DEALLOCATE ( R )
	  DEALLOCATE ( R_PNT )
	  DEALLOCATE ( LOG_R_SM )
!
	  DEALLOCATE ( ETA_COEF )
	  DEALLOCATE ( ESEC_COEF )
	  DEALLOCATE ( CHI_COEF )
	  DEALLOCATE ( F_COEF )
!
	  DEALLOCATE ( ETA )
	  DEALLOCATE ( ESEC )
	  DEALLOCATE ( CHI )
	  DEALLOCATE ( F )
!
	  DEALLOCATE ( DTAU )
	  DEALLOCATE ( MID_DTAU )
	  DEALLOCATE ( TA )
	  DEALLOCATE ( DD )
	  DEALLOCATE ( TC )
	  DEALLOCATE ( XM )
	  DEALLOCATE ( SOURCE )
	  DEALLOCATE ( HU )
	  DEALLOCATE ( HL )
	  DEALLOCATE ( J_INDX )
	  DEALLOCATE ( H_INDX )
	END IF
	NINS_SAV=NINS
!
! On the very first entry, we define the improved R grid, and allocate all
! data arrays.
!
	IF( FIRST_TIME .OR. .NOT. ALLOCATED(R) )THEN
!
	  ALLOCATE ( R(ND) )
	  ALLOCATE ( R_PNT(ND) )
!
	  ALLOCATE ( LOG_R_SM(ND_SM),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( ETA_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( ESEC_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( CHI_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( F_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	     WRITE(2,*)'Unable to allocate COEF memory in MOM_J_CMF_V6'
	     STOP
	  END IF
!
	  ALLOCATE ( ETA(ND) )
	  ALLOCATE ( ESEC(ND) )
	  ALLOCATE ( CHI(ND) )
!
	  ALLOCATE ( JNU(ND) )             ; JNU(1:ND)=0.0D0
	  ALLOCATE ( HNU(ND) )             ; HNU(1:ND)=0.0D0
	  ALLOCATE ( F(ND) )
!
	  ALLOCATE ( DTAU(ND) )
	  ALLOCATE ( MID_DTAU(ND) )
	  ALLOCATE ( TA(ND) )
	  ALLOCATE ( DD(ND) )
	  ALLOCATE ( TC(ND) )
	  ALLOCATE ( XM(ND) )
	  ALLOCATE ( SOURCE(ND) )
	  ALLOCATE ( HU(ND) )
	  ALLOCATE ( HL(ND) )
	  ALLOCATE ( COH_VEC(ND) )
!
	  ALLOCATE ( J_INDX(ND_SM) )
	  ALLOCATE ( H_INDX(ND_SM) )
!
	  FIRST_TIME=.FALSE.
	END IF
!
!
	IF(INIT)THEN
!
	  K=1
	  R(1)=R_SM(1)
          R_PNT(1)=1
	  DO I=1,ND_SM-1
            DELTA_R=(R_SM(I+1)-R_SM(I))/(NINS+1)
            DO J=1,NINS
              K=K+1
              R(K)=R(K-1)+DELTA_R
              R_PNT(K)=I
            END DO
            K=K+1
            R(K)=R_SM(I+1)
            R_PNT(K)=I
	  END DO
	  LOG_R_SM(1:ND_SM)=LOG(R_SM(1:ND_SM))
!
	  J_INDX(1:ND_SM)=0
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
	  H_INDX(1:ND_SM)=0
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
	END IF
!
!
!
! Interpolate quantities onto revised grid.
!
	IF(ND .GT. ND_SM)THEN
	  TA(1:ND_SM)=LOG(CHI_SM(1:ND_SM))
	  CALL MON_INT_FUNS_V2(CHI_COEF,TA,LOG_R_SM,ND_SM)
	  TA(1:ND_SM)=LOG(ESEC_SM(1:ND_SM))
	  CALL MON_INT_FUNS_V2(ESEC_COEF,TA,LOG_R_SM,ND_SM)
	  TA(1:ND_SM)=LOG(ETA_SM(1:ND_SM))
	  CALL MON_INT_FUNS_V2(ETA_COEF,TA,LOG_R_SM,ND_SM)
	  CALL MON_INT_FUNS_V2(F_COEF,F_SM,LOG_R_SM,ND_SM)

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
	  END DO
!
	ELSE
	  ESEC(1:ND)=ESEC_SM(1:ND)
	  CHI(1:ND)=CHI_SM(1:ND)
	  ETA(1:ND)=ETA_SM(1:ND)
	  F(1:ND)=F_SM(1:ND)
	END IF
!
! 
!
	IF(INIT)THEN
	  DO I=1,N_ERR_MAX
	    MOM_ERR_ON_FREQ(I)=0.0D0
	  END DO
	  MOM_ERR_CNT=0
	END IF
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
! NB: We solve J, not for r^2 J as in the spherical case.
!
! Compute optical depth scale.
!
	CALL DERIVCHI(DD,CHI,R,ND,METHOD)
	CALL NORDTAU(DTAU,CHI,R,R,DD,ND)
!
	DO I=2,ND
	  MID_DTAU(I)=0.5D0*(DTAU(I)+DTAU(I-1))
	END DO
!
! Compute vectors used to compute the flux vector H.
!
	DO I=1,ND-1
	  HU(I)=F(I+1)/DTAU(I)
	  HL(I)=F(I)/DTAU(I)
	END DO
!
! Compute the TRIDIAGONAL operators, and the RHS source vector.
!
	DO I=2,ND-1
	  TA(I)=-F(I-1)/DTAU(I-1)
	  TC(I)=-F(I+1)/DTAU(I)
	  DD(I)=-MID_DTAU(I)*(1.0D0-COH_VEC(I)) - (F(I)-F(I-1))/DTAU(I-1) - 
	1                                        (F(I)-F(I+1))/DTAU(I)
	  XM(I)=MID_DTAU(I)*SOURCE(I)
	END DO
!
! Second order boundary conditions.
!
	TC(1)=-F(2)/DTAU(1)
	DD(1)=(F(2)-F(1))/DTAU(1) -0.5D0*DTAU(1)*(1.0D0-COH_VEC(1)) - 
	1           HBC_J + HBC_S*COH_VEC(1)
	XM(1)=0.5D0*DTAU(1)*SOURCE(1)+HBC_S*SOURCE(1)
	TA(1)=0.0D0
!
	TA(ND)=-F(ND-1)/DTAU(ND-1)
	IF(DIF)THEN
	  DD(ND)=-(F(ND)-F(ND-1))/DTAU(ND-1)-0.5D0*DTAU(ND-1)*(1.0D0-COH_VEC(ND))
	  XM(ND)=DBB/3.0D0/CHI(ND)+0.5D0*DTAU(ND-1)*SOURCE(ND)
	ELSE
	  DD(ND)=-(F(ND)-F(ND-1))/DTAU(ND-1)-0.5D0*DTAU(ND-1)*(1.0D0-COH_VEC(ND))-IN_HBC
	  XM(ND)=IC*(0.25D0+0.5D0*IN_HBC)+0.5D0*DTAU(ND-1)*SOURCE(ND)
	END IF
	TC(ND)=0.0D0
!
! Solve for the radiation field along ray for this frequency.
!
	CALL THOMAS_RH(TA,DD,TC,XM,ND,1)
!
! Check that no negative mean intensities have been computed.
!
	IF(MINVAL(XM(1:ND)) .LE. 0)THEN
	   WRITE(47,*)'Freq=',FREQ
	   TA(1:ND)=XM(1:ND)/R(1:ND)/R(1:ND)
	   CALL WRITV(TA,ND,'XM Vec',47)
	   CALL WRITV(F,ND,'F Vec',47)
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
!
! Save J and compute H.
!
	JNU(1:ND)=XM(1:ND)
	DO I=1,ND-1
	  HNU(I)=HU(I)*XM(I+1)-HL(I)*XM(I)
	END DO
!
! Regrid derived J and H values onto small grid.
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
	RETURN
	END
