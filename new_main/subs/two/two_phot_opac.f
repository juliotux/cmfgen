!
! Routine to increment the Emissivity and Opacity for
! 2-photon processes.
!
	SUBROUTINE TWO_PHOT_OPAC(ETA,CHI,POPS,T,FREQ,ND,NT)
	USE TWO_PHOT_MOD
	IMPLICIT NONE
!
! Created 26-Jun-1998
!
	INTEGER NT,ND
	REAL*8 ETA(ND)		!Emisivity 
	REAL*8 CHI(ND)		!Opacity [in (10^10 cm)^-1]
	REAL*8 POPS(NT,ND)	!Atomic populations
	REAL*8 T(ND)		!in 10^4 K
!
	REAL*8 FREQ		!Current frequency (10^15 Hz)
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Local constants.
!
	REAL*8 PLANKS_CONSTANT	!cgs units
	REAL*8 PI
	REAL*8 CONST
	REAL*8 ETA_CONST	!Used to evaluate ETA
	REAL*8 CHI_CONST	!Used to evaluate CHI
	REAL*8 FREQ_B		!Frequency of other photon
	REAL*8 T1
!
! The 2-photon distribution functions, AY, are usually in wrtten in terms
! of the variable y=FREQ/MAX_FREQ, and which extends from 0 to 1.
! As the distribution, AY, is symmetric about y/2 we have defined the 
! variables  U=Y(1-Y) and FU=4*U
!
	REAL*8 AY,Y,U,FU
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	INTEGER J,L
	INTEGER NL,NUP
	
	PLANKS_CONSTANT=6.626D-27			!cgs units
	PI=4.0D0*ATAN(1.0D0)
!
! The factor 10^10 arises since R is units of 10^10cm, and we
! scale CHI (and ETA) so that R.CHI is dimensionless, and
! ETA/CHI has the dimensions of the Planck Function.
! We don't have to worry about the frequency units, as it multilplied
! by a ratio of 2 frequencies.
!
	CONST=1.0D+10*PLANKS_CONSTANT/4.0D0/PI
!
	DO J=1,N_TWO
	  IF(TWO_PHOT_AVAILABLE(J) .AND. FREQ .LT. FREQ_TWO(J))THEN
!
	    FREQ_B=FREQ_TWO(J)-FREQ
	    NL=LOW_LEV_TWO(J)
	    NUP=UP_LEV_TWO(J)
!
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
! The expressions for CHI and ETA are approximate. Their ratio gives the
! Planck function at depth, and the expression for ETA is exact at low
! densities and for small dilution factors.
!
	    ETA_CONST=CONST*FREQ/FREQ_TWO(J)
	    CHI_CONST=G_UP_TWO(J)*ETA_CONST/TWOHCSQ/FREQ**3
	    DO L=1,ND
	      ETA(L)=ETA(L) + ETA_CONST*AY*POPS(NUP,L)*FS_RAT_UP(L,J)
	      T1=EXP(-HDKT*FREQ_B/T(L))      
	      CHI(L)=CHI(L) + CHI_CONST*AY*( 
	1               POPS(NL,L)*FS_RAT_LOW(L,J)*T1/G_LOW_TWO(J)-
	1               POPS(NUP,L)*FS_RAT_UP(L,J)/G_UP_TWO(J) )
	    END DO
	  END IF
	END DO
!
	RETURN
	END
