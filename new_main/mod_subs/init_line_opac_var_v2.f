!
! Routine to:
!         (a) Assign storage locations for the variation with new lines. New lines are
!               indicated by NEW_LINE_STORAGE.
!         (b) Zero storage assoicated with new lines.
! This routine is called in the CONTINUUM section of the code, and is used
! in conjunction with SET_LINE_OPAC.
!
        SUBROUTINE INIT_LINE_OPAC_VAR_V2(LAST_LINE,LUER,ND,TX_OFFSET,MAX_SIM,NM)
	USE MOD_CMFGEN
 	USE OPAC_MOD
	USE CONTROL_VARIABLE_MOD
	USE LINE_MOD
	USE VAR_RAD_MOD
        IMPLICIT NONE
!
! Incorporated: 02-Jan-2104: LINE_QW_SUM is now a 2D array.
! Altered 29-Oct-2012: Changed to V2 but call is the same.
!                      Option to allow each line to have it own distinct stoagre location.
!                      Old method (but reprogammed) can still be used.
!
	INTEGER ND
	INTEGER NM
	INTEGER MAX_SIM
	INTEGER TX_OFFSET
	INTEGER LAST_LINE
	INTEGER LUER
!
	REAL*8 FL
	REAL*8 T1,T2
	REAL*8 NU_DOP
!
	INTEGER NUM_BNDS
	INTEGER NDEXT
	INTEGER I,J,K,L
	INTEGER ID
	INTEGER FREQ_INDX
	LOGICAL STORAGE_LOC_FOUND
!
	IF(ACCURATE)NDEXT=SIZE(TX_EXT,1)
	IF(COMPUTE_BA .AND. .NOT. LAMBDA_ITERATION)NUM_BNDS=SIZE(dZ,2)
!
! This section of the routine initializes arrays and data used in the CONTINUUM section.
!
	DO SIM_INDX=1,MAX_SIM
	  IF(NEW_LINE_STORAGE(SIM_INDX))THEN
!
! Zero arrays which is used to store the net rate (ZNET_SIM) and mean intensity
! for each line (JBAR_SIM).
!
	    DO I=1,ND
	      ZNET_SIM(I,SIM_INDX)=0.0D0
	      JBAR_SIM(I,SIM_INDX)=0.0D0
	    END DO
	    LINE_QW_SUM(1:ND,SIM_INDX)=0.0D0
!
! Decide if line is weak, and hence whether we can iterate on the net rates
! rather than use a full linearization.
!
	    IF(WEAK_WITH_NET)THEN
	      IF(USE_WEAK_TAU_LIM)THEN
!
! Compute Sobolev optical depth. Note: 2.998D+10 = C(km/s) / 1.0D+15
!
	        WEAK_LINE(SIM_INDX)=.TRUE.
	        T1=2.998D-10/FL_SIM(SIM_INDX)
	        DO I=1,ND
	          T2=ABS(CHIL_MAT(I,SIM_INDX))*T1*R(I)/V(I)
	          IF(T2 .GT. WEAK_TAU_LINE_LIMIT)WEAK_LINE(SIM_INDX)=.FALSE.
	        END DO
	      ELSE
!
! Compute opacity at line center.
!
	        T1=1.0D-15/1.77245385095516D0		!1.0D-15/SQRT(PI)
	        NU_DOP=FL_SIM(SIM_INDX)*12.85D0*SQRT( TDOP/AMASS_SIM(SIM_INDX) +
	1                        (VTURB/12.85D0)**2 )/2.998D+05
	        T2=T1/NU_DOP
	        WEAK_LINE(SIM_INDX)=.TRUE.
	        DO I=1,ND
	          IF( ABS(CHIL_MAT(I,SIM_INDX))*T2/ESEC(I) .GT. WEAK_LINE_LIMIT)
	1                              WEAK_LINE(SIM_INDX)=.FALSE.
	        END DO
	      END IF
	      IF(WEAK_LINE(SIM_INDX))NUM_OF_WEAK_LINES=NUM_OF_WEAK_LINES+1
	    ELSE
	      WEAK_LINE(SIM_INDX)=.FALSE.
	    END IF
! 
!
! Now need to determine the storage location for the 2 variation parameters.
! NB: Because of our coding, SIM_NL(SIM_INDX) must be set before we can 
!     do this.
!
! NB: TX_OFFSET refers to the number or arrays used (=5 if chi,eta, chi_old,
! eta_old, and esec) when computing the variation of J.
!
! We could provide two separate storage locations for each line. However,
! because of the use of super-levels, many transitions involve the same levels.
! We can therefore use the same storage area. To maintain consistency in the
! linearization we slpit transitions into two groups --- those where both
! levels are regarded as important, and those where at least one level is
! unimportant. Levels in thest two seaprate classes are kept distinct.
! The variable THIS_TRANS_IMP and the vector IMP_TRANS_VEC are used to
! distinguish betwene the two classes.
!
	    CALL TUNE(IONE,'VLSETUP')
	    IF(COMPUTE_BA .AND. .NOT. LAMBDA_ITERATION .AND. .NOT. WEAK_LINE(SIM_INDX))THEN
	      THIS_TRANS_IMP=IMP_VAR(SIM_NL(SIM_INDX)) .AND. IMP_VAR(SIM_NUP(SIM_INDX))
!
! Do the lower level. We first search to see if there is an existing storage location for level NL.
!
	      STORAGE_LOC_FOUND=.FALSE.
	      IF(.NOT. NEW_VAR_STORAGE_METHOD)THEN
                DO I=TX_OFFSET+1,NM
                  IF( (VAR_LEV_ID(I) .EQ. SIM_NL(SIM_INDX)) .AND. (IMP_TRANS_VEC(I) .EQV. THIS_TRANS_IMP) )THEN
                    LOW_POINTER(SIM_INDX)=I
                    VAR_IN_USE_CNT(I)=VAR_IN_USE_CNT(I)+1
                    STORAGE_LOC_FOUND=.TRUE.
	            EXIT
                  END IF
                END DO
	      END IF
!
	      IF(.NOT. STORAGE_LOC_FOUND)THEN
                I=TX_OFFSET+1
	        DO WHILE(VAR_LEV_ID(I) .NE. 0)
	          I=I+1
	          IF(I .GT. NM)THEN
	            WRITE(LUER,*)'Error in INIT_LINE_OPAC_VAR_V2 --- not enough storage locations'
	            WRITE(LUER,*)'LAST_LINE=',LAST_LINE-1
	            STOP
	          END IF
	        END DO
	        VAR_LEV_ID(I)=SIM_NL(SIM_INDX)
	        IMP_TRANS_VEC(I)=THIS_TRANS_IMP
	        LOW_POINTER(SIM_INDX)=I
	        VAR_IN_USE_CNT(I)=1
	      END IF
!
! Now do the upper level. We first search to see if there is an existing storage location for level NUP.
!
	      STORAGE_LOC_FOUND=.FALSE.
	      IF(.NOT. NEW_VAR_STORAGE_METHOD)THEN
                DO I=TX_OFFSET+1,NM
                  IF( (VAR_LEV_ID(I) .EQ. SIM_NUP(SIM_INDX)) .AND. (IMP_TRANS_VEC(I) .EQV. THIS_TRANS_IMP) )THEN
                    UP_POINTER(SIM_INDX)=I
                    VAR_IN_USE_CNT(I)=VAR_IN_USE_CNT(I)+1
                    STORAGE_LOC_FOUND=.TRUE.
	            EXIT
                  END IF
                END DO
	      END IF
!
	      IF(.NOT. STORAGE_LOC_FOUND)THEN
	        I=TX_OFFSET+1
	        DO WHILE(VAR_LEV_ID(I) .NE. 0 )
	          I=I+1
	          IF(I .GT. NM)THEN
	            WRITE(LUER,*)'Error in INIT_LINE_OPAC_VAR_V2 --- not enough storage locations'
	            WRITE(LUER,*)'LAST_LINE=',LAST_LINE
	            STOP
	          END IF
	        END DO
	        VAR_LEV_ID(I)=SIM_NUP(SIM_INDX)
	        IMP_TRANS_VEC(I)=THIS_TRANS_IMP
	        UP_POINTER(SIM_INDX)=I
	        VAR_IN_USE_CNT(I)=1
	      END IF
!
! 
!
! Zero the appropriate dNL and dNUP matrices in TX and TVX IFF they are
! not already in use by another line.
!
	      J=LOW_POINTER(SIM_INDX)
	      IF(VAR_IN_USE_CNT(J) .EQ. 1)THEN	  !Just this line using store.
!	        TX(:,:,J)=0.0D0
!	        TVX(:,:,J)=0.0D0
	        CALL ZERO_2D_MAT(TX(1,1,J),ND,ND)
	        CALL ZERO_2D_MAT(TVX(1,1,J),ND-1,ND)
	        IF(ACCURATE)CALL ZERO_2D_MAT(TX_EXT(1,1,J),NDEXT,NDEXT)
	        IF(ACCURATE)CALL ZERO_2D_MAT(TVX_EXT(1,1,J),NDEXT-1,NDEXT)
!	        IF(ACCURATE)TX_EXT(:,:,J)=0.0D0
!	        IF(ACCURATE)TVX_EXT(:,:,J)=0.0D0
	      END IF
!
	      J=UP_POINTER(SIM_INDX)
	      IF(VAR_IN_USE_CNT(J) .EQ. 1)THEN	  !Just this line using store.
!	        TX(:,:,J)=0.0D0
!	        TVX(:,:,J)=0.0D0
	        CALL ZERO_2D_MAT(TX(1,1,J),ND,ND)
	        CALL ZERO_2D_MAT(TVX(1,1,J),ND-1,ND)
	        IF(ACCURATE)CALL ZERO_2D_MAT(TX_EXT(1,1,J),NDEXT,NDEXT)
	        IF(ACCURATE)CALL ZERO_2D_MAT(TVX_EXT(1,1,J),NDEXT-1,NDEXT)
!	        IF(ACCURATE)TX_EXT(:,:,J)=0.0D0
!	        IF(ACCURATE)TVX_EXT(:,:,J)=0.0D0
	      END IF
!
! Zero the appropriate dCHIL and dETAL matrices in dZ, since we will no
! longer be including the variation of the deleted line.
!
	      J=LOW_POINTER(SIM_INDX)
	      IF(VAR_IN_USE_CNT(J) .EQ. 1)THEN	  !Just this line using store.
	        CALL ZERO_dZ(dZ,J,NM,NUM_BNDS*ND,MAX_SIM)
!	        dZ(J,:,:,:)=0.0D0		!NM,NUM_BNDS,ND,MAX_SIM
	      END IF
	      J=UP_POINTER(SIM_INDX)
	      IF(VAR_IN_USE_CNT(J) .EQ. 1)THEN	  !Just this line using store.
	        CALL ZERO_dZ(dZ,J,NM,NUM_BNDS*ND,MAX_SIM)
!	        dZ(J,:,:,:)=0.0D0		!NM,NUM_BNDS,ND,MAX_SIM
	      END IF
!
! Ensure dZ for this line is zeroed.
!
	      CALL ZERO_2D_MAT(dZ(1,1,1,SIM_INDX),NM,NUM_BNDS*ND)
!	      dZ(:,:,:,SIM_INDX)=0.0D0	!NM,NUM_BNDS,ND,MAX_SIM	
!
	    END IF			!BA computed and weak line check
	    CALL TUNE(ITWO,'VLSETUP')
	  END IF                        !New line storage requires
	END DO
!
	RETURN
	END
!
	SUBROUTINE ZERO_dZ(dZ,J,NM,NX,M)
	IMPLICIT NONE
	INTEGER J,NM,NX,M
	REAL*8 dZ(NM,NX,M)
!
	INTEGER L,K
!
!$OMP PARALLEL DO
	DO L=1,M
	  DO K=1,NX
	     dZ(J,K,L)=0.0D0
	  END DO
	END DO
!
	RETURN
	END
