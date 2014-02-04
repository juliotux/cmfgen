!
! Simple subroutine to estimate amount of gamma-ray energy absorbed. This routine
! assumes the pure absorption approximation.
!
! It will be replaced by a more sophisticated version.
!
	SUBROUTINE DO_GAM_ABS_APPROX(DECAY_ENERGY,ND)
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Created: 29-MAy-2009
!
	INTEGER ND
	INTEGER ND_EXT,NC_EXT,NP_EXT
!
	REAL*8 DECAY_ENERGY(ND)
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
	CHI(1:ND)=0.06D0*TA(1:ND)*DENSITY(1:ND)*1.0D+10
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
	CALL SUB_GAM_ABS_APPROX(R_EXT,ND_EXT,NC_EXT,NP_EXT,R,V,CHI,DECAY_ENERGY,ND)
!
	RETURN
	END
! 
	SUBROUTINE SUB_GAM_ABS_APPROX(R,ND,NC,NP,SM_R,SM_V,SM_CHI,DECAY_ENERGY,SM_ND)
	IMPLICIT NONE
!
	INTEGER SM_ND
	INTEGER ND,NC,NP
	REAL*8 R(ND)
!
	REAL*8 SM_R(SM_ND)
	REAL*8 SM_V(SM_ND)
	REAL*8 SM_CHI(SM_ND)
	REAL*8 DECAY_ENERGY(SM_ND)
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
	REAL*8 IC,T1,DBB,HBC_J,HBC_S,INBC
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
	CALL MON_INTERP(ETA,ND,IONE,R,ND,DECAY_ENERGY,SM_ND,SM_R,SM_ND)
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
	RJ=CHI*XM
!
	DO I=1,ND
	  TA(I)=ETA(I)*R(I)*R(I)
	  TB(I)=RJ(I)*R(I)*R(I)
	END DO
	CALL LUM_FROM_ETA(TA,R,ND)
	CALL LUM_FROM_ETA(TB,R,ND)
	T1=16.0D0*ATAN(1.0D0)*1.0D+30/3.826D+33
	OPEN(UNIT=10,FILE='check_edep.dat',STATUS='UNKNOWN',ACTION='WRITE')
	  WRITE(10,'(A,ES13.5,A)')'!Radiactive energy emitted is :',SUM(TA)/T1,' Lsun'
	  WRITE(10,'(A,ES13.5,A)')'!Radiactive energy absorbed is:',SUM(TB)/T1,' Lsun'
	  DO I=1,ND
	    WRITE(10,'(5ES14.6)')V(I),RJ(I),ETA(I)
	  END DO
!
	  WRITE(10,'(A)')'!'
	  WRITE(10,'(A)')'! Small grid'
	  WRITE(10,'(A)')'!'
	  CALL MON_INTERP(TA,SM_ND,IONE,SM_R,SM_ND,RJ,ND,R,ND)
	  DO I=1,SM_ND
	    WRITE(10,'(5ES14.6)')SM_V(I),TA(I),DECAY_ENERGY(I)
	  END DO
	CLOSE(UNIT=10)
!
	DECAY_ENERGY(1:SM_ND)=TA(1:SM_ND)
!
	RETURN
	END
