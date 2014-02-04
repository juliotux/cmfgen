	MODULE BA_J_DATA_MOD_V2
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
! Module which contains dJdN integrated over a small frequency band over
! which the continuum cross-sections are assumed not to change.
!
! VJ_R = Int ( dJ*EXP(_hv/kT)/v) dv
! VJ_P = Int ( dJ/v)
! VJ_C = Int ( dJ)
! RJ_SUM = Int (J)
!
	TYPE dJ
	  REAL*8, POINTER :: VJ_R(:,:,:)
	  REAL*8, POINTER :: VJ_P(:,:,:)
	  REAL*8, POINTER :: VJ_T(:,:,:)
	END DO
	TYPE dJ((BA_NUM_SPECIES*MAX_IONS_PER_SPECIES)
!
	REAL*8, ALLOCATABLE :: RJ_SUM(:)
!
	END MODULE BA_J_DATA_MOD
!
! Subroutine to increment the variation matrix BA due to the variation of J.
!
! BA_PAR is used for the diagonal variation only. It is updated on each call
! rather than BA to improve numerical stability. BA_PAR should contain terms
! of similar size. BA_PAR will need to be added to BA after every approximately
! every 50 frequencies. In this way the BA and BAION matrices should suffer
! less cancelation effects due to the addition of large positive and negative
! terms.
!
! Utilizing the fact that consecutive frequency terms should be correlated, and
! hence similar in size. While some minor cancellation, should be much less
! then adding to full BA matrix in which terms have arbitrary size.
!
	SUBROUTINE BA_UPDATE_V6(
	1              VJ,VCHI,VETA,
	1              CHI_CONT,ESEC,T,POPS,RJ,
	1              NU,FQW,NEW_CONT,FINAL_FREQ,
	1              dJ_CHK_FAC,NION,
	1              NT,NUM_BNDS,ND,DST,DEND)
	USE BA_J_DATA_MOD
	USE STEQ_DATA_MOD 
	IMPLICIT NONE
!
! Altered 05-Apr-2001 :: Change to V6
!                        Now STEQ_DATA_MOD. STurcture BA data array.
!                        Important variables utilized.
! Altered 02-Mar-2001 :: Changed to V5
!                        NIV and LNK_F_TO_IV inserted in call.
!                        Altered to handle important variables. With
!                        the LNK_F_TO_IV array appropriately set, the routine
!                        will also recover the case with out important 
!                        variableas.
!
! Altered: 17-Sep-1997 :: QFV_R and QFV_P installed so that BA is not updated
!                           for every frequency. Call changed so updated to
!                           V4. If BA is to be updated every frequency, NEW_FREQ
!                           and FINAL_FREQ must both be true. In this case
!                           routine should give ``identical'' results to
!                           V3.
! Altered: 01-Feb-1997 :: dJ_CHK_FAC put in call (changed from V2 to V3)
!                           Replaces parameter RMAX_FAC.
!                           The smaller dJ_CHK_FAC, the more accurate the
!                           computation of the BAION matrix.
!                           dj_CHK_FAC is normally around 1.0D-04.
!                           Larger values give less accuracy for BAION,
!                           but allow faster computation.
!
! Altered: 16-Aug-1996 :: COMP_VEC installed to improve vectorization.
!                           Improvement will depend on how many times
!                           innermost loop is executed.
!
! Altered: 24-May-1996 :: IONE, ITWO, RMAX_FAC inserted
!                         DABS changed to ABS.
! Created: 28-Feb-1995
!
	INTEGER NION
  	INTEGER NT,NUM_BNDS,ND,DST,DEND
	REAL*8 VJ(NT,NUM_BNDS,ND)
	REAL*8 POPS(NT,ND)
	REAL*8 VCHI(NT,ND)
	REAL*8 VETA(NT,ND)
	REAL*8 RJ(ND)
	REAL*8 dJ_CHK_FAC
!
	REAL*8 CHI_CONT(ND)
	REAL*8 ESEC(ND)
	REAL*8 T(ND)
!
	REAL*8 NU
	REAL*8 FQW
!
	LOGICAL NON_ZERO_VJ(NT)
!
! NEW_CONT indicates that this is the first frequency of a new continuum band
! in which the continuum cross-sections are constant. FINAL_FREQ indicates
! that it is the last frequency of a continuum band.
!
	LOGICAL FINAL_FREQ
	LOGICAL NEW_CONT
!
	REAL*8 COMP_VEC(NT)
!
! Constants for opacity etc.
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
!
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
	REAL*8 T1,T2,QFV_T
	INTEGER I,J,K,L,LS,IOS,JJ,ID
	INTEGER DIAG_INDX,BNDST,BNDEND
!
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
!
	CALL TUNE(IONE,'BA_UP')
	DIAG_INDX=(NUM_BNDS+1)/2
!
	NON_ZERO_VJ(:)=.TRUE.
!	DO I=1,NT
!	  IF(VJ(I,DIAG_INDX,1) .EQ. 0 .AND.
!	1    VJ(I,DIAG_INDX,ND/2) .EQ. 0 .AND.
!	1    VJ(I,DIAG_INDX,ND) .EQ. 0)NON_ZERO_VJ(I)=.FALSE.
!	END DO
! 
	IF(.NOT. ALLOCATED(dJ(1)%VJ_R))THEN
	  DO ID=1,NION
	    IF(SE(ID)%XzV_PRES)THEN
	      NX=SE(ID)%N_SE+2
	    ELSE
	      NX=1
	    END IF
	  END DO
	  ALLOCATE (SE(ID)%VJ_R(NX,NUM_BNDS,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (SE(ID)%VJ_P(NX,NUM_BNDS,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (VJ_T(NT,NUM_BNDS,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (RJ_SUM(ND),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    I=ERROR_LU()
	    WRITE(I,*)'Error in BA_UPDATE_V4'

	    WRITE(I,*)'Unable to allocate required dynamic memory'
	    STOP
	  END IF
	END IF
!
! dJ_CHK_FAC 
!
	IF(dJ_CHK_FAC .LT. 1.0D-10 .OR. dJ_CHK_FAC .GT. 0.1)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error in BA_UPDATE_V4'
	  WRITE(I,*)'Invalid value for dJ_CHK_FAC: dJ_CHK_FAC=',dJ_CHK_FAC
	  STOP
	END IF
!
! Perform the frequency integral of dJ over. Procedure depends on whether
! this is a new frequency of part of a band. In oder to minimize computation,
! a band which is a single frequency (i.e. NEW_CONT and FINAL_FREQ both TRUE)
! is treated as a special case.
!
	IF(NEW_CONT .AND. .NOT. FINAL_FREQ)THEN
	  T2=FQW/NU
	  DO L=DST,DEND
	    RJ_SUM(L)=T2*RJ(L)
	  END DO
!
	  T2=FQW/NU
	  DO ID=1,NION
	    IF(SE(ID)%ZxV_PRES)THEN
  	      DO L=DST,DEND
	        T1=T2*EXP(-HDKT*NU/T(L))
	        dJ(ID)%VJ_R(:,:,L)=T1*SE(ID)%VJ(:,:,L)
	        dJ(ID)%VJ_P(:,:,L)=T2*SE(ID)%VJ(:,:,L)
	      END DO
	    END IF
	  END DO
!
	  DO ID=1,NION
	    IF(SE(ID)%ZxV_PRES .AND. .NOT. SE(ID)%IMPURITY_SPECIES)THEN
	      DO L=DST,DEND
	        DO K=1,NUM_BNDS
	          IF(K .EQ. DIAG_INDX)THEN
	            dJ(ID)%VJ_T(:,K,L)=FQW*( RJ(L)*SE(ID)%VCHI(:,L) - SE(ID)%VETA(:,L) + 
	1                       (CHI_CONT(L)-ESEC(L))*VJ(:,K,L) )
	          ELSE
	            dJ%VJ_T(:,K,L)=FQW*(CHI_CONT(L)-ESEC(L))*VJ(:,K,L)
	          END IF
	        END DO
	      END DO
	    END IF
	  END DO
	  RETURN
!
	ELSE IF(NEW_CONT)THEN		!and hence `single frequency' band
	  T2=FQW/NU
	  DO L=1,ND
	    T1=T2*EXP(-HDKT*NU/T(L))
	    DO ID=1,NION
	      IF(SE(ID)%XzV_PRES)THEN
	        SE(ID)%QFV_R(:,L)=SE(ID)%QFV_R(:,L)*T1-SE(ID)%QFV_P(:,L)*T2	!1:SE(ID)%N_IV
	      END IF
	    END DO
	  END DO
	ELSE
	  CALL TUNE(1,'BA_UP_VJ')
	  T2=FQW/NU
	  DO L=DST,DEND
	    T1=T2*EXP(-HDKT*NU/T(L))
	    DO K=1,NUM_BNDS
	      DO I=1,NT
	        IF(NON_ZERO_VJ(I))THEN
	          VJ_R(I,K,L)=VJ_R(I,K,L)+T1*VJ(I,K,L)
	          VJ_P(I,K,L)=VJ_P(I,K,L)+T2*VJ(I,K,L)
	        END IF
	      END DO
	    END DO
	    RJ_SUM(L)=RJ_SUM(L)+T2*RJ(L)
	  END DO
	  CALL TUNE(2,'BA_UP_VJ')
!
	  CALL TUNE(1,'BA_UP_VT')
	  DO L=DST,DEND
	    T1=CHI_CONT(L)-ESEC(L)
	    T2=FQW*T1
	    DO K=1,NUM_BNDS
	      IF(K .EQ. DIAG_INDX)THEN
	        DO I=1,NT
	          IF(NON_ZERO_VJ(I))VJ_T(I,K,L)= VJ_T(I,K,L) + 
	1              FQW*(RJ(L)*VCHI(I,L) - VETA(I,L)+ 
	1              T1*VJ(I,K,L) )
	        END DO
	      ELSE
	        DO I=1,NT
	          IF(NON_ZERO_VJ(I))
	1               VJ_T(I,K,L)= VJ_T(I,K,L) + T2*VJ(I,K,L)
	        END DO
	      END IF
	    END DO
	  END DO
	  CALL TUNE(2,'BA_UP_VT')
	END IF
!
	IF(NEW_CONT .AND. FINAL_FREQ)THEN
!
! Done in this way to ensure better cancellation of large terms. Note that
! we loop over all L. This minimizes paging.
!
	  DO L=DST,DEND					!S.E. equation depth
	    QFV_T=FQW*(CHI_CONT(L)-ESEC(L))
	    DO K=1,NUM_BNDS
	      IF(K .EQ. DIAG_INDX)THEN
   	        DO  J=1,NT	  	  	  	!Variable
	          BA_T_PAR(J,L)=BA_T_PAR(J,L) + (
	1            FQW*(RJ(L)*VCHI(J,L)-VETA(J,L)) + QFV_T*VJ(J,K,L) )
	        END DO
	      ELSE
   	        DO  J=1,NT	  	  	  !Variable
	          BA_T(J,K,L)=BA_T(J,K,L)+QFV_T*VJ(J,K,L)
	        END DO
	      END IF
 	    END DO
C
	    BNDST=MAX( 1+DIAG_INDX-L, 1 )
	    BNDEND=MIN( ND+DIAG_INDX-L, NUM_BNDS )
	    DO K=BNDST,BNDEND	  			!Variable depth.
	      LS=L+K-DIAG_INDX
                    DO JJ=1,SE(VID)%N_SE	  	  	  !Variable
  	              J=SE(ID)%LNKF_TO_IV(SE(VID)%LNK_TO_F(JJ))
   	      DO VID=1,SE(ID)%N_SE
!
   	      DO VID=1,NION
!
	        DO J=1,SE(VID)%N_SE+2
	          JJ=SE(VID)%LNK_IV_TO_F(J)
	          COMP_VEC(J)=ABS(SE(VID)%VJ(J,K,L)*POPS(JJ,LS))-RJ(L)*dJ_CHK_FAC
	        END DO
!
	        DO ID=1,NION
	          IF(SE(ID)%XzV_PRES .AND. K .EQ. DIAG_INDX)THEN
   	          DO VID=1,SE(ID)%N_SE
	            IF(.NOT. SE(ID)$IMPURITY_SPECIES .OR. VID .EQ. ID)THEN
                      DO J=1,SE(VID)%N_SE	  	  	  !Variable
  	                JJ=SE(ID)%LNKF_TO_IV(SE(VID)%LNK_TO_F(JJ))
	                IF( COMP_VEC(J) .GE. 0.0D0 )THEN
	                DO  I=1,SE(ID)%N_SE	  	  !Which S.E.
	                  SE(ID)%BA_PAR(I,JJ,L)=SE(ID)%BA_PAR(I,JJ,L) +
	1                           SE(ID)%QFV_R(I,L)*VJ(J,K,L)
	                END DO
	            END IF
	          END DO
!
	        ELSE IF(SE(ID)%XzV_PRES)THEN
   	          DO JJ=1,SE(ID)%N_IV	  	  	  !Variable
  	            J=SE(ID)%LNK_TO_F(JJ)
	            IF( COMP_VEC(J) .GE. 0.0D0 )THEN
	              DO I=1,SE(ID)%N_SE	  	  !Which S.E.
	                 SE(ID)%BA(I,JJ,K,L)=SE(ID)%BA(I,JJ,K,L) +
	1                          SE(ID)%QFV_R(I,L)*VJ(J,K,L)
	              END DO
	            END IF
	          END DO	!Loop over variable
	        END IF
	      END DO		!Loop over species
	    END DO		!Loop over band
	  END DO					!Do DST to DEND
C
C Update BA matrices for several frequencies at once.
C
	ELSE IF(FINAL_FREQ)THEN
	  CALL TUNE(1,'BA_FF')
	  DO L=DST,DEND					!S.E. equation depth
	    BNDST=MAX( 1+DIAG_INDX-L, 1 )
	    BNDEND=MIN( ND+DIAG_INDX-L, NUM_BNDS )
	    DO K=BNDST,BNDEND	  			!Variable depth.
	      IF(K .EQ. DIAG_INDX)THEN
   	        DO  J=1,NT	  	  	  	!Variable
	          BA_T_PAR(J,L)=BA_T_PAR(J,L) + VJ_T(J,K,L)
	        END DO
	      ELSE
   	        DO  J=1,NT	  	  	  	!Variable
	          BA_T(J,K,L)=BA_T(J,K,L) + VJ_T(J,K,L)
	        END DO
	      END IF
	    END DO
C
C NB: We use VJ_P to compute COMP_VEC as this is defined a Int[ (VJ/v) dv].
C            RJ_SUM is defined in the same way (i..e., Int[ (J/v) dv]
C
	    DO K=BNDST,BNDEND	  			!Variable depth.
	      LS=L+K-DIAG_INDX
   	      DO ID=1,NION
!
	        DO J=1,SE(VID)%N_SE+2
	          JJ=SE(VID)%LNK_IV_TO_F(J)
	          COMP_VEC(J)=ABS(SE(VID)%VJ(J,K,L)*POPS(JJ,LS))-RJ(L)*dJ_CHK_FAC
	        END DO
!
	        ID_BEG=1; ID_MAX=NION
	        IF(SE(VID)$IMPURITY_SPECIES)THEN
	          ID_BEG=VID, ID_MAX=VID
	        END IF
!
	        DO VID=ID_BEG,ID_MAX
	          IF(K .EQ. DIAG_INDX .AND. SE(ID)%XzV_PRES)THEN
   	            DO JJ=1,SE(ID)%N_SE+2	  	  	  !Variable
	              IF( COMP_VEC(JJ) .GE. 0.0D0 )THEN
  	                J=SE(ID)%LNKTO_IV(SE(VID)%LNK_TO_F(JJ))
	                DO  I=1,SE(ID)%N_SE	  	  !Which S.E.
	                  SE(ID)%BA_PAR(I,JJ,L)=SE(ID)%BA_PAR(I,JJ,L)+
	1                   ( SE(ID)%QFV_R(I,L)*SE(VID)%VJ_R(J,K,L) - SE(ID)%QFV_P(I,L)*SE(VID)%VJ_P(J,K,L) )
	                END DO
	              END IF
	            END DO
!
	          ELSE IF(SE(ID)%XzV_PRES)THEN
   	            DO JJ=1,SE(ID)%N_SE+2	  	  	  !Variable
	              IF( COMP_VEC(JJ) .GE. 0.0D0 )THEN
  	                J=SE(ID)%LNKTO_IV(SE(VID)%LNK_TO_F(JJ))
	                DO  I=1,SE(ID)%N_SE	  	  !Which S.E.
	                  SE(ID)%BA(I,JJ,K,L)=SE(ID)%BA(I,JJ,K,L)+
	1                   ( SE(ID)%QFV_R(I,L)*SE(VID)%VJ_R(J,K,L) - SE(ID)%QFV_P(I,L)*SE(VID)%VJ_P(J,K,L) )
	                END DO
	              END IF
	            END DO
!
	        END IF
	      END DO			!Loop over ion
	    END DO			!Loop over K (band index)
	  END DO			!Loop over depth (DST to DEND)
	  CALL TUNE(2,'BA_FF')
	END IF
	CALL TUNE(ITWO,'BA_UP')
C
	RETURN
	END
