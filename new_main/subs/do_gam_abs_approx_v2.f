!
! Simple subroutine to estimate amount of gamma-ray energy absorbed. This routine
! assumes the pure absorption approximation.
!
! It will be replaced by a more sophisticated version.
!
	SUBROUTINE DO_GAM_ABS_APPROX_V2(LOCAL_ABS_ENERGY,TOTAL_DECAY_ENERGY,KINETIC_DECAY_ENERGY,ND)
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Altered: 27-May-2016 -- CHI is now multiplied by CLUMP_FAC so that clumping is correctly accounted for.
! Altered: 29-Jan-2015 -- LOCAL_ABS_ENEGRY added to call. Fixed bug in luminosity calculation.
! Altered: 06-Jan-2015 -- KINETIC_DECAY_ENERGY included in the call.
!                         Changed to V2. Kinetic energy is assumed to be
!                         absorbed locally.
! Created: 29-May-2009
!
	INTEGER ND
	INTEGER ND_EXT,NC_EXT,NP_EXT
!
	REAL*8 LOCAL_ABS_ENERGY(ND)     	!Returned - energy absorbed locally
	REAL*8 TOTAL_DECAY_ENERGY(ND)           !Passed   - total energey EMITTED locally
	REAL*8 KINETIC_DECAY_ENERGY(ND)         !Passed   - energy local emitted as kinetic energy
	REAL*8 TA(ND)
	REAL*8 CHI(ND)
!
	REAL*8 R_EXT(3*ND-6)
!
	INTEGER I,J,ISPEC
!
	ND_EXT=3*ND-6
	NC_EXT=50
	NP_EXT=ND_EXT+NC_EXT
!
! Determine number of electron per baryon.
!
	TA(1:ND)=0.5D0                        !Number of electrons per baryon
	DO ISPEC=1,NUM_SPECIES
	  IF('HYD' .EQ. SPECIES(ISPEC))THEN
	   TA(1:ND)=0.5D0*(1.0D0+POP_SPECIES(1:ND,ISPEC)/POP_ATOM(1:ND))
	  END IF
	END DO
!
! Compute the absorbative opacity.
!
	CHI(1:ND)=0.06D0*TA(1:ND)*DENSITY(1:ND)*CLUMP_FAC(1:ND)*1.0D+10
!
! Calculate a larger R grid. We attempt to keep the outer grid spacing small.
!
	R_EXT(1)=R(1)
	R_EXT(2)=R(1)+(R(2)-R(1))/3.0D0
	R_EXT(3)=R(1)+(R(3)-R(1))/3.0D0
	R_EXT(4)=R(1)+(R(3)-R(1))/1.5D0
	I=3
	J=5
	DO I=3,ND-2
	  R_EXT(J)=R(I)
	  R_EXT(J+1)=R(I)+(R(I+1)-R(I))/3.0D0
	  R_EXT(J+2)=R(I)+(R(I+1)-R(I))/1.5D0
	  J=J+3
	END DO
	R_EXT(ND_EXT-1)=R(ND)+(R(ND-1)-R(ND))/3.0D0
	R_EXT(ND_EXT)=R(ND)
!
	DO I=1,ND_EXT-1
	  IF(R_EXT(I+1) .GE. R_EXT(I))THEN
	    WRITE(6,*)'Error -- invalid extended R grid in DO_GAM_ABS_APPROX'
	    WRITE(6,*)'Depth is',I
	    WRITE(6,*)R_EXT
	    STOP
	  END IF
	END DO
!
! Now do the transfer.
!
	CALL SUB_GAM_ABS_APPROX(R_EXT,ND_EXT,NC_EXT,NP_EXT,R,V,CHI,
	1          LOCAL_ABS_ENERGY,TOTAL_DECAY_ENERGY,KINETIC_DECAY_ENERGY,ND)
!
	RETURN
	END
! 
	SUBROUTINE SUB_GAM_ABS_APPROX(R,ND,NC,NP,SM_R,SM_V,SM_CHI,
	1              LOCAL_ABS_ENERGY,TOTAL_DECAY_ENERGY,KINETIC_DECAY_ENERGY,SM_ND)
	IMPLICIT NONE
!
	INTEGER SM_ND
	INTEGER ND,NC,NP
	REAL*8 R(ND)
!
	REAL*8 SM_R(SM_ND)
	REAL*8 SM_V(SM_ND)
	REAL*8 SM_CHI(SM_ND)
	REAL*8 LOCAL_ABS_ENERGY(SM_ND)
	REAL*8 TOTAL_DECAY_ENERGY(SM_ND)
	REAL*8 KINETIC_DECAY_ENERGY(SM_ND)
!
	INTEGER, PARAMETER :: IONE=1
!
	REAL*8 P(NP)
	REAL*8 V(ND)
!
	REAL*8 Z(ND),RJ(ND),Q(ND),F(ND)
	REAL*8 DTAU(ND),XM(ND),ETA(ND)
	REAL*8 TA(ND),TB(ND),TC(ND)
	REAL*8 SOURCE(ND),CHI(ND),DCHIDR(ND),THETA(ND)
	REAL*8 JQW(ND,NP),KQW(ND,NP)
	REAL*8 WM(ND,ND),FB(ND,ND)
!
	REAL*8 IC,T1,T2,T3,CONV_FAC,DBB,HBC_J,HBC_S,INBC
	LOGICAL THK_CONT
	CHARACTER(LEN=6) METHOD
        CHARACTER(LEN=9) INNER_BND_METH
!
	INTEGER I
	INTEGER IOS
	CHARACTER(LEN=132) STRING
        EXTERNAL JTRPWGT,HTRPWGT,KTRPWGT,NTRPWGT
!
	METHOD='LOGLOG'
	INNER_BND_METH='DIFFUSION'
	THK_CONT=.FALSE.
	DBB=0.0D0
!
	CALL MON_INTERP(  V,ND,IONE,R,ND,SM_V,        SM_ND,SM_R,SM_ND)
	CALL MON_INTERP(CHI,ND,IONE,R,ND,SM_CHI,      SM_ND,SM_R,SM_ND)
	TA(1:SM_ND)=TOTAL_DECAY_ENERGY(1:SM_ND)-KINETIC_DECAY_ENERGY(1:SM_ND)
	CALL MON_INTERP(ETA,ND,IONE,R,ND,TA,SM_ND,SM_R,SM_ND)
!
	CALL IMPAR(P,R,R(ND),NC,ND,NP)
        CALL GENANGQW(JQW,R,P,WM(1,1),WM(1,3),WM(1,5),NC,ND,NP,JTRPWGT,.FALSE.)
        CALL GENANGQW(KQW,R,P,WM(1,1),WM(1,3),WM(1,5),NC,ND,NP,KTRPWGT,.FALSE.)
!
! The emissivity is actually ETA/4PI, but we would need to multiply again to get
! the absorbed energy.
!
	SOURCE=ETA/CHI
	CALL FQCOMP_IBC_V2(TA,TB,TC,XM,DTAU,R,Z,P,Q,F,
	1      SOURCE,CHI,DCHIDR,JQW,KQW,DBB,HBC_J,HBC_S,
	1      INBC,IC,THK_CONT,INNER_BND_METH,NC,ND,NP,METHOD)
!
! RJ is the absorbed gamma-ray energy.
!
	RJ=CHI*XM
!
! Interpolate absorbed energy on to the CMFGEN grid, and add in the locally deposited
! kinetic energy.
!
	CALL MON_INTERP(TA,SM_ND,IONE,SM_R,SM_ND,RJ,ND,R,ND)
	LOCAL_ABS_ENERGY(1:SM_ND)=TA(1:SM_ND)+KINETIC_DECAY_ENERGY(1:SM_ND)
!
	DO I=1,SM_ND
	  TA(I)=TOTAL_DECAY_ENERGY(I)*SM_R(I)*SM_R(I)
	  TB(I)=LOCAL_ABS_ENERGY(I)*SM_R(I)*SM_R(I)
	  TC(I)=KINETIC_DECAY_ENERGY(I)*SM_R(I)*SM_R(I)
	END DO
	CALL LUM_FROM_ETA(TA,SM_R,SM_ND)
	CALL LUM_FROM_ETA(TB,SM_R,SM_ND)
	CALL LUM_FROM_ETA(TC,SM_R,SM_ND)
!
! The conversion factor:
!   (a) 4pi r^2 dr  --> 4pi .10^30 (as R is in units of 10^10 cm).
!   (b) The extra factor of 4 arises as ATAN(1.0D0) is pi/4.
!   (c) /Lsun to convert to solar luminosities
!
	CONV_FAC=16.0D0*ATAN(1.0D0)*1.0D+30/3.826D+33
	T1=SUM(TA)*CONV_FAC
	T2=SUM(TB)*CONV_FAC
	T3=SUM(TC)*CONV_FAC
!
	OPEN(UNIT=10,FILE='check_edep.dat',STATUS='UNKNOWN',ACTION='WRITE')
	  WRITE(10,'(A)')'!'
	  WRITE(10,'(A,ES13.5,A)')'!            Radioactive energy emitted is :',T1,' Lsun'
	  WRITE(10,'(A,ES13.5,A)')'!            Radioactive energy absorbed is:',T2,' Lsun'
	  WRITE(10,'(A,ES13.5,A)')'!Fraction of radioactive energy absorbed is:',T2/T1
	  WRITE(10,'(A,ES13.5,A)')'!      Fraction absorbed that is kinetic is:',T3/T2
	  WRITE(10,'(A)')'!'
!
	  DO I=1,SM_ND
	    WRITE(10,'(5ES14.6)')SM_V(I),LOCAL_ABS_ENERGY(I),TOTAL_DECAY_ENERGY(I),KINETIC_DECAY_ENERGY(I)
	  END DO
!
! Outout the results on the fine grid -- useful for checking.
! Get total decay energy on grid -- not just that due to gamma-rays.
!
	  CALL MON_INTERP(RJ,ND,IONE,R,ND,LOCAL_ABS_ENERGY,SM_ND,SM_R,SM_ND)
	  CALL MON_INTERP(TB,ND,IONE,R,ND,TOTAL_DECAY_ENERGY,SM_ND,SM_R,SM_ND)
	  WRITE(10,'(A)')'!'
	  WRITE(10,'(A)')'! Large grid'
	  WRITE(10,'(A)')'! V(kms)   E(non local)   E(local)   Eg(local)'
	  WRITE(10,'(A)')'!'
	  DO I=1,ND
	    WRITE(10,'(5ES14.6)')V(I),RJ(I),TB(I),ETA(I)
	  END DO
!
	CLOSE(UNIT=10)
!
	RETURN
	END
