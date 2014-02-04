	MODULE BA_J_DATA_MOD_V4
	IMPLICIT NONE
C
C Module which contains dJdN integrated over a small frequency band over
C which the continuum cross-sections are assumed not to change.
C
C VJ_R = Int ( dJ*EXP(_hv/kT)/v) dv
C VJ_P = Int ( dJ/v)
C VJ_C = Int ( dJ)
C RJ_SUM = Int (J)
C
	REAL*8, ALLOCATABLE :: VJ_R(:,:,:)
	REAL*8, ALLOCATABLE :: VJ_P(:,:,:)
	REAL*8, ALLOCATABLE :: VJ_T(:,:,:)
	REAL*8, ALLOCATABLE :: RJ_SUM(:)
C
	END MODULE BA_J_DATA_MOD_V4
C
C Subroutine to increment the variation matrix BA due to the variation of J.
C
C BA_PAR is used for the diagonal variation only. It is updated on each call
C rather than BA to improve numerical stability. BA_PAR should contain terms
C of similar size. BA_PAR will need to be added to BA after every approximately
C every 50 frequencies. In this way the BA and BAION matrices should suffer
C less cancelation effects due to the addition of large positive and negative
C terms.
C
C Utilizing the fact that consecutive frequency terms should be correlated, and
C hence similar in size. While some minor cancellation, should be much less
C then adding to full BA matrix in which terms have arbitrary size.
C
	SUBROUTINE BA_UPDATE_V4(BA,BA_PAR,QFV_R,QFV_P,
	1              VJ,VCHI,VETA,
	1              CHI_CONT,ESEC,T,POPS,RJ,
	1              NU,FQW,NEW_CONT,FINAL_FREQ,
	1              dJ_CHK_FAC,NT,NUM_BNDS,ND,DST,DEND)
	USE BA_J_DATA_MOD_V4
	IMPLICIT NONE
C
C Altered: 17-Sep-1997 :: QFV_R and QFV_P installed so that BA is not updated
C                           for every frequency. Call changed so updated to
C                           V4. If BA is to be updated every frequency, NEW_FREQ
C                           and FINAL_FREQ must both be true. In this case
C                           routine should give ``identical'' results to
C                           V3.
C Altered: 01-Feb-1997 :: dJ_CHK_FAC put in call (changed from V2 to V3)
C                           Replaces parameter RMAX_FAC.
C                           The smaller dJ_CHK_FAC, the more accurate the
C                           computation of the BAION matrix.
C                           dj_CHK_FAC is normally around 1.0D-04.
C                           Larger values give less accuracy for BAION,
C                           but allow faster computation.
C
C Altered: 16-Aug-1996 :: COMP_VEC installed to improve vectorization.
C                           Improvement will depend on how many times
C                           innermost loop is executed.
C
C Altered: 24-May-1996 :: IONE, ITWO, RMAX_FAC inserted
C                         DABS changed to ABS.
C Created: 28-Feb-1995
C
  	INTEGER NT,NUM_BNDS,ND,DST,DEND
	REAL*8 BA(NT,NT,NUM_BNDS,ND)
	REAL*8 BA_PAR(NT,NT,ND)
	REAL*8 QFV_R(NT,ND)
	REAL*8 QFV_P(NT,ND)
	REAL*8 VJ(NT,NUM_BNDS,ND)
	REAL*8 POPS(NT,ND)
	REAL*8 VCHI(NT,ND)
	REAL*8 VETA(NT,ND)
	REAL*8 RJ(ND)
	REAL*8 dJ_CHK_FAC
C
	REAL*8 CHI_CONT(ND)
	REAL*8 ESEC(ND)
	REAL*8 T(ND)
C
	REAL*8 NU
	REAL*8 FQW
C
C NEW_CONT indicates that this is the first frequency of a new continuum band
C in which the continuum cross-sections are constant. FINAL_FREQ indicates
C that it is the last frequency of a continuum band.
C
	LOGICAL FINAL_FREQ
	LOGICAL NEW_CONT
C
	REAL*8 COMP_VEC(NT)
C
C Constants for opacity etc.
C
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
C
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
C
	REAL*8 T1,T2,QFV_T
	INTEGER I,J,K,L,LS,IOS
	INTEGER DIAG_INDX,BNDST,BNDEND
C
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
C
	CALL TUNE(IONE,'BA_UP')
	DIAG_INDX=(NUM_BNDS+1)/2
C
	IF(.NOT. ALLOCATED(VJ_R))THEN
	  ALLOCATE (VJ_R(NT,NUM_BNDS,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (VJ_P(NT,NUM_BNDS,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (VJ_T(NT,NUM_BNDS,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (RJ_SUM(ND),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    I=ERROR_LU()
	    WRITE(I,*)'Error in BA_UPDATE_V4'
	    WRITE(I,*)'Unable to allocate required dynamic memory'
	    STOP
	  END IF
	END IF
C
C dJ_CHK_FAC 
C
	IF(dJ_CHK_FAC .LT. 1.0D-10 .OR. dJ_CHK_FAC .GT. 0.1)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error in BA_UPDATE_V4'
	  WRITE(I,*)'Invalid value for dJ_CHK_FAC: dJ_CHK_FAC=',dJ_CHK_FAC
	  STOP
	END IF
C
C Perform the frequency integral of dJ over. Procedure depends on whether
C this is a new frequency of part of a band. In oder to minimize computation,
C a band which is a single frequency (i.e. NEW_CONT and FINAL_FREQ both TRUE)
c is treated as a special case.
C
	IF(NEW_CONT .AND. .NOT. FINAL_FREQ)THEN
	  T2=FQW/NU
	  DO L=DST,DEND
	    T1=T2*EXP(-HDKT*NU/T(L))
	    VJ_R(:,:,L)=T1*VJ(:,:,L)
	    VJ_P(:,:,L)=T2*VJ(:,:,L)
	    RJ_SUM(L)=T2*RJ(L)
	  END DO
	  DO L=DST,DEND
	    DO K=1,NUM_BNDS
	      IF(K .EQ. DIAG_INDX)THEN
	        VJ_T(:,K,L)=FQW*( RJ(L)*VCHI(:,L) - VETA(:,L) + 
	1                         (CHI_CONT(L)-ESEC(L))*VJ(:,K,L) )
	      ELSE
	        VJ_T(:,K,L)=FQW*(CHI_CONT(L)-ESEC(L))*VJ(:,K,L)
	      END IF
	    END DO
	  END DO
	  RETURN
	ELSE IF(NEW_CONT)THEN		!and hence `single frequency' band
	  T2=FQW/NU
	  DO L=1,ND
	    T1=T2*EXP(-HDKT*NU/T(L))
	    QFV_R(:,L)=QFV_R(:,L)*T1-QFV_P(:,L)*T2		!1:NT
	  END DO
	ELSE
	  T2=FQW/NU
	  DO L=DST,DEND
	    T1=T2*EXP(-HDKT*NU/T(L))
	    VJ_R(:,:,L)=VJ_R(:,:,L)+T1*VJ(:,:,L)
	    VJ_P(:,:,L)=VJ_P(:,:,L)+T2*VJ(:,:,L)
	    RJ_SUM(L)=RJ_SUM(L)+T2*RJ(L)
	  END DO
	  DO L=DST,DEND
	    DO K=1,NUM_BNDS
	      IF(K .EQ. DIAG_INDX)THEN
	        VJ_T(:,K,L)= VJ_T(:,K,L) + FQW*(
	1         RJ(L)*VCHI(:,L) - VETA(:,L)+ 
	1                      (CHI_CONT(L)-ESEC(L))*VJ(:,K,L) )
	      ELSE
	        VJ_T(:,K,L)= VJ_T(:,K,L) + 
	1                       FQW*(CHI_CONT(L)-ESEC(L))*VJ(:,K,L)
	      END IF
	    END DO
	  END DO
	END IF
C
	IF(NEW_CONT .AND. FINAL_FREQ)THEN
C
C Done in this way to ensure better cancellation of large terms. Note that
C we loop over all L. This minimizes paging.
C
	  DO L=DST,DEND					!S.E. equation depth
	    QFV_T=FQW*(CHI_CONT(L)-ESEC(L))
	    DO K=1,NUM_BNDS
	      IF(K .EQ. DIAG_INDX)THEN
   	        DO  J=1,NT	  	  	  	!Variable
	          BA_PAR(NT,J,L)=BA_PAR(NT,J,L) + (
	1            FQW*(RJ(L)*VCHI(J,L)-VETA(J,L)) + QFV_T*VJ(J,K,L) )
	        END DO
	      ELSE
   	        DO  J=1,NT	  	  	  !Variable
	          BA(NT,J,K,L)=BA(NT,J,K,L)+QFV_T*VJ(J,K,L)
	        END DO
	      END IF
 	    END DO
C
	    BNDST=MAX( 1+DIAG_INDX-L, 1 )
	    BNDEND=MIN( ND+DIAG_INDX-L, NUM_BNDS )
	    DO K=BNDST,BNDEND	  			!Variable depth.
	      LS=L+K-DIAG_INDX
	      DO J=1,NT
	        COMP_VEC(J)=ABS(VJ(J,K,L)*POPS(J,LS))-RJ(L)*dJ_CHK_FAC
	      END DO
	      IF(K .EQ. DIAG_INDX)THEN
   	         DO  J=1,NT	  	  	  !Variable
                   IF( COMP_VEC(J) .GE. 0.0D0 )THEN
	             DO  I=1,NT-1	  	  !Which S.E.
	               BA_PAR(I,J,L)=BA_PAR(I,J,L) +
	1                               QFV_R(I,L)*VJ(J,K,L)
	             END DO
	           END IF
	        END DO
	      ELSE
   	        DO  J=1,NT	  	  	  !Variable
                  IF( COMP_VEC(J) .GE. 0.0D0)THEN
	            DO  I=1,NT-1	  	  !Which S.E.
	              BA(I,J,K,L)=BA(I,J,K,L)+QFV_R(I,L)*VJ(J,K,L)
	            END DO
	          END IF
	        END DO
	      END IF
	    END DO
	  END DO					!Do DST to DEND
C
C Update BA matrices for several frequencies at once.
C
	ELSE IF(FINAL_FREQ)THEN
	  DO L=DST,DEND					!S.E. equation depth
	    BNDST=MAX( 1+DIAG_INDX-L, 1 )
	    BNDEND=MIN( ND+DIAG_INDX-L, NUM_BNDS )
	    DO K=BNDST,BNDEND	  			!Variable depth.
	      IF(K .EQ. DIAG_INDX)THEN
   	        DO  J=1,NT	  	  	  	!Variable
	          BA_PAR(NT,J,L)=BA_PAR(NT,J,L) + VJ_T(J,K,L)
	        END DO
	      ELSE
   	        DO  J=1,NT	  	  	  	!Variable
	          BA(NT,J,K,L)=BA(NT,J,K,L) + VJ_T(J,K,L)
	        END DO
	      END IF
	    END DO
C
C NB: We use VJ_P to compute COMP_VEC as this is defined a Int[ (VJ/v) dv].
C            RJ_SUM is defined in the same way (i..e., Int[ (J/v) dv]
C
	    DO K=BNDST,BNDEND	  			!Variable depth.
	      LS=L+K-DIAG_INDX
	      DO J=1,NT
	        COMP_VEC(J)=ABS(VJ_P(J,K,L)*POPS(J,LS))-RJ_SUM(L)*dJ_CHK_FAC
	      END DO
	      IF(K .EQ. DIAG_INDX)THEN
   	         DO  J=1,NT	  	  	  !Variable
                   IF( COMP_VEC(J) .GE. 0.0D0 )THEN
	             DO  I=1,NT-1	  	  !Which S.E.
	               BA_PAR(I,J,L)=BA_PAR(I,J,L) +
	1               ( QFV_R(I,L)*VJ_R(J,K,L) - QFV_P(I,L)*VJ_P(J,K,L) )
	             END DO
	           END IF
	        END DO
	      ELSE
   	        DO  J=1,NT	  	  	  !Variable
                  IF( COMP_VEC(J) .GE. 0.0D0)THEN
	            DO  I=1,NT-1	  	  !Which S.E.
	              BA(I,J,K,L)=BA(I,J,K,L)+
	1               ( QFV_R(I,L)*VJ_R(J,K,L) - QFV_P(I,L)*VJ_P(J,K,L) )
	            END DO
	         END IF
	        END DO
	      END IF
	    END DO
	  END DO					!Do DST to DEND
	END IF
	CALL TUNE(ITWO,'BA_UP')
C
	RETURN
	END
