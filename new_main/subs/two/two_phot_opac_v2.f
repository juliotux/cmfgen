!
! Routine to increment the Emissivity and Opacity for 2-photon processes.
!
	SUBROUTINE TWO_PHOT_OPAC_V2(ETA,CHI,POPS,T,NU,FREQ,TWO_PHOTON_METHOD,
	1              NCF,ND,NT,LUEDD,IREC_ST)
	USE TWO_PHOT_MOD
	IMPLICIT NONE
!
! Altered 14-Jul-2015: Added two options to improve two photon absorption.
!                        USE_TWO  - Uses actual radiation field from EDDFACTOR file.
!                        USE_W    - Uses dilution factor.
!	
! Created 26-Jun-1998
!
	INTEGER LUEDD
	INTEGER IREC_ST
	INTEGER NCF
	INTEGER NT,ND
	REAL*8 NU(NCF)
	REAL*8 ETA(ND)		!Emisivity 
	REAL*8 CHI(ND)		!Opacity [in (10^10 cm)^-1]
	REAL*8 POPS(NT,ND)	!Atomic populations
	REAL*8 T(ND)		!in 10^4 K
	REAL*8 FREQ		!Current frequency (10^15 Hz)
!
	CHARACTER(LEN=*) TWO_PHOTON_METHOD
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Local vectors and constants.
!
	REAL*8 TA(ND)
	REAL*8 TB(ND)
	REAL*8 PLANKS_CONSTANT	!cgs units
	REAL*8 PI
	REAL*8 CONST
	REAL*8 ETA_CONST	!Used to evaluate ETA
	REAL*8 CHI_CONST	!Used to evaluate CHI
	REAL*8 FREQ_B		!Frequency of other photon
	REAL*8 T1,T2
!
! The 2-photon distribution functions, AY, are usually in wrtten in terms
! of the variable y=FREQ/MAX_FREQ, and which extends from 0 to 1.
! As the distribution, AY, is symmetric about y/2 we have defined the 
! variables  U=Y(1-Y) and FU=4*U
!
	REAL*8 AY,Y,U,FU
!
	INTEGER LUER,ERROR_LU,GET_TWO_INDX
	EXTERNAL ERROR_LU,GET_TWO_INDX
!
	INTEGER J,L,ML
	INTEGER NL,NUP
!
	TWO_METHOD=TWO_PHOTON_METHOD
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
! The exact expressions for the two photon emissivity and opacity are:
!
! Eta=  (hv/4pi) (nu/nu_o) A(v) N_u [1 + Y(nu_B)]
! Chi=  (c^2/2hv^3) (hv/4pi) (nu/nu_o) A(v) .
!	    [ N_l (g_u/g_l) Y(nu_B) - N_u {1+Y(nu_B)} ]
! where
!       Y(v)=J(v)/[2hv^3/c^2]	   
!
	    ETA_CONST=CONST*FREQ/FREQ_TWO(J)
	    CHI_CONST=G_UP_TWO(J)*ETA_CONST/TWOHCSQ/FREQ**3
	    IF(TWO_METHOD .EQ. 'USE_RAD')THEN
	      IF(LST_FREQ_INDX_TWO(J) .EQ. 0)WRITE(6,*)'Using USE_RAD method for two-photon'
	      ML=GET_TWO_INDX(LST_FREQ_INDX_TWO(J),FREQ_B,NU,NCF)
	      IF(FREQ_B .LT. NU(NCF))THEN
	        READ(LUEDD,REC=IREC_ST+ML)TB,T2
	        RJ_TWO(:,J)=TB*(FREQ_B/NU(NCF))**2	!Black body extrapolation
	      ELSE
	        READ(LUEDD,REC=IREC_ST+ML-1)TA,T2
	        IF(FREQ_B .GT. T2)THEN
	          WRITE(6,*)'Error - invalid fequency -- SET_TWO_OPAC'
	          WRITE(6,*)FREQ_B,T2
	        END IF
	        READ(LUEDD,REC=IREC_ST+ML)TB,T2
	        IF(FREQ_B .LT. T2)THEN
	          WRITE(6,*)'Error - invalid fequency -- SET_TWO_OPAC'
	          WRITE(6,*)FREQ_B,T2
	        END IF
	        T1=(FREQ_B-NU(ML+1))/(NU(ML)-NU(ML+1))
	        RJ_TWO(:,J)=(1.0D0-T1)*TA+T1*TB
	      END IF
	      LST_FREQ_INDX_TWO(J)=ML
	      T1=TWOHCSQ*(FREQ_B**3)
	      RJ_TWO(:,J)=RJ_TWO(:,J)/T1 
	    ELSE IF(TWO_METHOD .EQ. 'USE_W')THEN
	      DO L=1,ND
	        T1=EXP(-HDKT*FREQ_B/T(L))
	        RJ_TWO(L,J)=W_TWO(L)*T1/(1.0D0-T1)
	      END DO
	    END IF
!
	    IF(TWO_METHOD .EQ. 'USE_RAD' .OR. TWO_METHOD .EQ. 'USE_W')THEN
	      DO L=1,ND
	        ETA(L)=ETA(L) + ETA_CONST*AY*POPS(NUP,L)*FS_RAT_UP(L,J)*(1.0D0+RJ_TWO(L,J))
	        CHI(L)=CHI(L) + CHI_CONST*AY*( 
	1                 POPS(NL,L)*FS_RAT_LOW(L,J)*RJ_TWO(L,J)/G_LOW_TWO(J)-
	1                 POPS(NUP,L)*FS_RAT_UP(L,J)/G_UP_TWO(J)*(1.0D0+RJ_TWO(L,J)) )
	      END DO
	    ELSE
!
! The expressions for CHI and ETA are approximate. Their ratio gives the
! Planck function at depth, and the expression for ETA is exact at low
! densities and for small dilution factors.
!
	      DO L=1,ND
	        ETA(L)=ETA(L) + ETA_CONST*AY*POPS(NUP,L)*FS_RAT_UP(L,J)
	        T1=EXP(-HDKT*FREQ_B/T(L))      
	        CHI(L)=CHI(L) + CHI_CONST*AY*( 
	1                 POPS(NL,L)*FS_RAT_LOW(L,J)*T1/G_LOW_TWO(J)-
	1                 POPS(NUP,L)*FS_RAT_UP(L,J)/G_UP_TWO(J) )
	      END DO
	    END IF
	  END IF
	END DO
!
	RETURN
	END
