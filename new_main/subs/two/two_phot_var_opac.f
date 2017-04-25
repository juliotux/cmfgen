!
! Routine to compute the variation of CHI and ETA due to 2-photon
! processes. We ignore any temperature dependence (which may arise
! from the use of super levels).
!
	SUBROUTINE TWO_PHOT_VAR_OPAC(VETA,VCHI,POPS,T,FREQ,ND,NT)
	USE TWO_PHOT_MOD
	IMPLICIT NONE
!
! Altered 01-Oct-2015 : Added TWO_METHOD option (passed by module). i
! Created 26-Jun-1998
!
	INTEGER NT,ND
	REAL*8 VETA(NT,ND)
	REAL*8 VCHI(NT,ND)
	REAL*8 POPS(NT,ND)
	REAL*8 T(ND)
!
	REAL*8 FREQ
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Local constants.
!
	REAL*8 h
	REAL*8 PI
	REAL*8 CONST
	REAL*8 ETA_CONST
	REAL*8 CHI_CONST
	REAL*8 T1
	REAL*8 FREQ_B
!
! See TWO_PHOT_OPAC for definitions
!
	REAL*8 AY,Y,U,FU
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	INTEGER J,L
	INTEGER NL,NUP
	
	h=6.626D-27			!cgs units
	PI=4.0D0*ATAN(1.0D0)
	CONST=1.0D+10*H/4.0D0/PI
!
	DO J=1,N_TWO
	  IF(TWO_PHOT_AVAILABLE(J) .AND. FREQ .LT. FREQ_TWO(J))THEN
!
	    FREQ_B=FREQ_TWO(J)-FREQ
	    NL=LOW_LEV_TWO(J)
	    NUP=UP_LEV_TWO(J)
	    IF(TYPE_TWO(J) .EQ. 1)THEN
	      Y=FREQ/FREQ_TWO(J)
	      U=Y*(1.0D0-Y)
	      FU=4.0D0*U
	      AY=24.56D0*COEF_TWO(J,1)*( U*(1.0D0-FU**0.8D0) +
	1                 0.88D0*(U**1.53D0)*(FU**0.8D0) )
	    ELSE
	      WRITE(6,'(/,1X,A)')'Error in TWO_PHOT_VAR_OPAC -- unrecognized type for TWO_PHOTON transition'
	      STOP
	    END IF
!
	    ETA_CONST=CONST*FREQ/FREQ_TWO(J)
	    CHI_CONST=G_UP_TWO(J)*ETA_CONST/TWOHCSQ/FREQ**3
	    IF(TWO_METHOD .EQ. 'OLD_DEFAULT')THEN
	      DO L=1,ND
	        VETA(NUP,L)=VETA(NUP,L) + ETA_CONST*AY*FS_RAT_UP(L,J)
	        T1=EXP(-HDKT*FREQ_B/T(L))      
	        VCHI(NL,L)=VCHI(NL,L)   + CHI_CONST*AY*FS_RAT_LOW(L,J)*T1/G_LOW_TWO(J)
	        VCHI(NUP,L)=VCHI(NUP,L) - CHI_CONST*AY*FS_RAT_UP(L,J)/G_UP_TWO(J)
	      END DO
	    ELSE
	      DO L=1,ND
	        VETA(NUP,L)=VETA(NUP,L) + ETA_CONST*AY*FS_RAT_UP(L,J)*(1.0D0+PHOT_OC_TWO(L,J))
	        VCHI(NL,L)=VCHI(NL,L)   + CHI_CONST*AY*FS_RAT_LOW(L,J)*PHOT_OC_TWO(L,J)/G_LOW_TWO(J)
	        VCHI(NUP,L)=VCHI(NUP,L) - CHI_CONST*AY*FS_RAT_UP(L,J)/G_UP_TWO(J)*(1.0D0+PHOT_OC_TWO(L,J))
	      END DO
	    END IF
	  END IF
	END DO
!
	RETURN
	END
