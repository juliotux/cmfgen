!
! Data module for storing quanties on the FINE R_GRID.
! saved between subroutine calls..
!
	MODULE MOD_FINE_R_GRID
	IMPLICIT NONE
!
! To be dimenensioned ND_SM where ND_SM is the size of the R grid
! as passed to FINE_R_GRID.
!
	REAL*8, ALLOCATABLE :: LOG_R_SM(:)
	REAL*8, ALLOCATABLE :: R_MID_SM(:)
	REAL*8, ALLOCATABLE :: RSQ_GAM_SM(:)
!
! Dimensioned ND_SM,4
!
	REAL*8, ALLOCATABLE :: V_COEF(:,:)
	REAL*8, ALLOCATABLE :: SIGMA_COEF(:,:)
	REAL*8, ALLOCATABLE :: ESEC_COEF(:,:)
	REAL*8, ALLOCATABLE :: CHI_COEF(:,:)
	REAL*8, ALLOCATABLE :: ETA_COEF(:,:)
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
	REAL*8, ALLOCATABLE :: GAM_REL(:)
!
	REAL*8, ALLOCATABLE :: JNU(:)
	REAL*8, ALLOCATABLE :: H_ON_J(:)
	REAL*8, ALLOCATABLE :: K_ON_J(:)
	REAL*8, ALLOCATABLE :: N_ON_J(:)
!
	REAL*8, ALLOCATABLE :: JNU_MID(:)
	REAL*8, ALLOCATABLE :: HNU_MID(:)
	REAL*8, ALLOCATABLE :: KNU_MID(:)
	REAL*8, ALLOCATABLE :: NNU_MID(:)

	REAL*8, ALLOCATABLE :: N_ON_H_MID(:)
	REAL*8, ALLOCATABLE :: RSQN_MID_ON_RSQJ(:)
!
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
	REAL*8 VDOP_FRAC_SAV
	LOGICAL NEW_R_GRID
	DATA VDOP_FRAC_SAV/-1000001.1D0/    !Absurd value
!
	END MODULE MOD_FINE_R_GRID
!
!
!
	SUBROUTINE FINE_R_GRID(V_SM,SIGMA_SM,R_SM,VDOP_VEC,VDOP_FRAC,INIT,ND_SM)
	USE MOD_FINE_R_GRID
	IMPLICIT NONE
!
	INTEGER ND_SM
	REAL*8 R_SM(ND_SM)
	REAL*8 V_SM(ND_SM)
	REAL*8 SIGMA_SM(ND_SM)
	REAL*8 VDOP_VEC(ND_SM)
	REAL*8 VDOP_FRAC
	LOGICAL INIT
!
	REAL*8 SPEED_OF_LIGHT
	REAL*8 C_KMS
	REAL*8 T1,T2
	REAL*8  DELTA_R
	INTEGER IT1,I,J,K
	INTEGER IOS
	EXTERNAL SPEED_OF_LIGHT
!
	NEW_R_GRID=.FALSE.
	C_KMS=1.0D-05*SPEED_OF_LIGHT()
	IF(INIT .AND. ALLOCATED(R))THEN
	  DO I=1,ND_SM
	    IF(LOG(R_SM(I)) .NE. LOG_R_SM(I))THEN
	      NEW_R_GRID=.TRUE.
	      WRITE(171,*)'Updating RGRID in MOM_RGRID'
	      EXIT
	    END IF
	  END DO
          IF(VDOP_FRAC .NE. VDOP_FRAC_SAV)NEW_R_GRID=.TRUE.
	ELSE
	  NEW_R_GRID=.FALSE.
	  RETURN
	END IF
!
	IF(NEW_R_GRID)THEN
	  ALLOCATE (R_MID_SM(ND_SM))
	  ALLOCATE (RSQ_GAM_SM(ND_SM))
	  DO I=1,ND_SM
	   RSQ_GAM_SM(I)=R_SM(I)*R_SM(I)/SQRT(1.0D0-(V_SM(I)/C_KMS)**2)
	  END DO
	  DO I=1,ND_SM-1
	    R_MID_SM(I)=0.5D0*(R_SM(I)+R_SM(I+1))
	  END DO
	END IF
!
! Deallocate all arrayes if we have changed VDOP_FRAC. This will only
! be done in testing this routine (e.g., using DISPGEN).
!
	IF(ALLOCATED(R) .AND. NEW_R_GRID)THEN
	  DEALLOCATE ( R )
	  DEALLOCATE ( R_PNT )
	  DEALLOCATE ( LOG_R_SM )
	  DEALLOCATE ( V )
	  DEALLOCATE ( SIGMA )
	  DEALLOCATE ( J_INDX )
	  DEALLOCATE ( H_INDX )
	END IF
	VDOP_FRAC_SAV=VDOP_FRAC
!
! On the very first entry, we define the improved R grid, and allocate all
! data arrays.
!
	IF( .NOT. ALLOCATED(R) )THEN
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
	  ALLOCATE ( R(ND) )
	  ALLOCATE ( R_PNT(ND) )
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
	  IF(IOS .NE. 0)THEN
	     WRITE(2,*)'Unable to allocate COEF memory in EXT_RGRID'
	     STOP
	  END IF
!
	  ALLOCATE ( V(ND) )
	  ALLOCATE ( V(ND) )
	  ALLOCATE ( SIGMA(ND) )
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
! Define the links between the FINE grid, and the old origunal grid. We need 2 links:
! One for the grid points, and one for the mid points of the old grid.
!
	  ALLOCATE ( J_INDX(ND_SM) );       J_INDX(1:ND_SM)=0
	  ALLOCATE ( H_INDX(ND_SM) );       H_INDX(1:ND_SM)=0
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
	END IF
	DEALLOCATE (V_COEF,SIGMA_COEF)
!
	DO I=1,ND
	  GAM_REL(I)=SQRT(1.0D0/(1.0D0-(V(I)/C_KMS)**2))
	END DO
!
	RETURN
	END
