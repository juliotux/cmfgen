!
! Routine to compute the Upward and Downward rates for all
! 2-photon transitions. The rates only depend on the radiation
! field, and are NOT scaled by the actual atomic populations
!
	SUBROUTINE TWO_PHOT_RATE(T,RJ,FREQ,FQW,ND,NT)
	USE TWO_PHOT_MOD
	IMPLICIT NONE
!
! Altered 19-Aug-2015: Additional options installed to hand 2-photon opacity (cur_hmi, 8-Jul-2015).
! Created 26-Jun-1998
!
	INTEGER NT,ND
	REAL*8 T(ND)
	REAL*8 RJ(ND)
!
	REAL*8 FREQ		!In units of 10^15Hz
	REAL*8 FQW		!In Hz
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Local constants.
!
	REAL*8 DOWN_CONST
	REAL*8 UP_CONST
	REAL*8 T1
	REAL*8 ALPHA_A		!2 h v^3 / c^2
	REAL*8 FREQ_B		!Frequency of other photon
!
! See TWO_PHOT_OPAC for definitions of AY,Y, U, and FU.
!
	REAL*8 AY,Y,U,FU
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	INTEGER J,L
	INTEGER NL,NUP
!
	DO J=1,N_TWO
	  IF(TWO_PHOT_AVAILABLE(J) .AND. FREQ .LT. FREQ_TWO(J))THEN
	    FREQ_B=FREQ_TWO(J)-FREQ
	    NL=LOW_LEV_TWO(J)
	    NUP=UP_LEV_TWO(J)
	    ALPHA_A=TWOHCSQ*(FREQ**3)
	    IF(TYPE_TWO(J) .EQ. 1)THEN
	      Y=FREQ/FREQ_TWO(J)
	      U=Y*(1.0D0-Y)
	      FU=4.0D0*U
	      AY=24.56D0*COEF_TWO(J,1)*( U*(1.0D0-FU**0.8D0) +
	1                 0.88D0*(U**1.53D0)*(FU**0.8D0) )
	    ELSE
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Error in TWO_PHOT_OPAC'
	      WRITE(LUER,*)'Unrecognized two photon transition'
	      STOP
	    END IF
!
! The constant = 0.5/1.0E+15. The factor of 0.5 comes about since
! 2 photons are emitted per electron transition. The factor of 10^15
! arrises since FREQ_TWO is in units of 10^15Hz, while FQW is in units
! of Hz.
!
	    DOWN_CONST=0.5D-15*AY/FREQ_TWO(J)
	    UP_CONST=DOWN_CONST*G_UP_TWO(J)/G_LOW_TWO(J)
	    IF(TWO_METHOD .EQ. 'OLD_DEFAULT')THEN
	      UP_CONST=UP_CONST/ALPHA_A
	      DO L=1,ND
	        DOWN_RATE_TWO(L,J)=DOWN_RATE_TWO(L,J) +
	1          DOWN_CONST*FS_RAT_UP(L,J)*(1.0D0+RJ(L)/ALPHA_A)*FQW
	        T1=EXP(-HDKT*FREQ_B/T(L))
	        UP_RATE_TWO(L,J)=UP_RATE_TWO(L,J) +
	1          UP_CONST*FS_RAT_LOW(L,J)*T1*RJ(L)*FQW
	      END DO
	    ELSE
	      DO L=1,ND
	        T1=RJ(L)/ALPHA_A
	        DOWN_RATE_TWO(L,J)=DOWN_RATE_TWO(L,J) +
	1          DOWN_CONST*FS_RAT_UP(L,J)*(1.0D0+T1+PHOT_OC_TWO(L,J)+T1*PHOT_OC_TWO(L,J))*FQW
	        UP_RATE_TWO(L,J)=UP_RATE_TWO(L,J) +
	1          UP_CONST*FS_RAT_LOW(L,J)*T1*PHOT_OC_TWO(L,J)*FQW
	      END DO
	    END IF
	  END IF
	END DO
!
	RETURN
	END
