!
! Routine to adjust the Radiative equilibrium equation for adiabatic cooling.
! This routine:
!
!  (1) Increments STEQ if adiabatic cooling is to be allowed for
!  (2) Increments the variation matrix [BA] if it is being computed, and
!       if adiabatic cooling is being included.
!  (3) Evaluates the adiabatic cooling terms, splitting them into 2 terms ---
!       the dTdR term, and the velocity term. This is done for diagnostic
!       purposes.
!
 	SUBROUTINE EVAL_ADIABATIC_V3(AD_CR_V,AD_CR_DT,
	1                     POPS,AVE_ENERGY,HDKT,COMPUTE_BA,INCL_ADIABATIC,
	1                     DIAG_INDX,NUM_BNDS,NT,ND)
!
	USE MOD_CMFGEN
 	USE STEQ_DATA_MOD
 	IMPLICIT NONE
!
! Altered  21-Jun-2004 : Changed to version V3.
!                          Changed to use simple linear differencing.
!                          Changed to incorporate advection terms from rate equations
!                            in electron cooling equation.
!                          Changed to allow for changes in internal (excitation) energy.
! Altered  16-Dec-2002 : DTDR computed using forward differencing.
! Altered  12-Apr-2001 : Changed to use STEQ_DATA_MOD
!                        Changed to V2 (call changed).
! Altered  15-Mar-2001 : Minor bug fix in BA matrix at d=1.
! Created 256-Jul-1994
!
	INTEGER NT
	INTEGER ND
	INTEGER DIAG_INDX
	INTEGER NUM_BNDS
!
! Output:
!
	REAL*8 AD_CR_V(ND)
	REAL*8 AD_CR_DT(ND)
!
! Input:
!
	REAL*8 POPS(NT,ND)
	REAL*8 AVE_ENERGY(NT)
	REAL*8 HDKT
!
! Local vectors.
!
	REAL*8 A(ND)
	REAL*8 B(ND)
	REAL*8 C(ND)
	REAL*8 D(ND)
	REAL*8 GAMMA(ND)
	REAL*8 INT_EN(ND)
	REAL*8 COL_EN(ND)
	REAL*8 WORK(ND)
!
	REAL*8 ION_EN(NT)
	REAL*8 TOT_ENERGY(NT)
!
! Local variables.
!
	LOGICAL COMPUTE_BA,INCL_ADIABATIC
!
	INTEGER ERROR_LU
	REAL*8 BOLTZMANN_CONSTANT,FUN_PI
	EXTERNAL BOLTZMANN_CONSTANT,FUN_PI,ERROR_LU
!
	REAL*8 SCALE,T1,T2,PI
	INTEGER I,J,L
	INTEGER LUER
	INTEGER ISPEC
	INTEGER ID
	LOGICAL WRITE_CHK
!
! A full linearization is now obsolete, but check to make sure.
!
	LUER=ERROR_LU()
	IF(NUM_BNDS .EQ. ND)THEN
	  WRITE(LUER,*)'Error --- EVAL_ADIABATIC can''t handle a full linearization'
	  STOP
	END IF
!
! Compute the total excitation energy of each level.
!
	TOT_ENERGY(1:NT)=0.0D0
	ION_EN(1:NT)=0.0D0
	DO ISPEC=1,NUM_SPECIES
	  T1=0.0D0
	  T2=0.0D0
	  DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	    T2=T2+AVE_ENERGY(ATM(ID)%EQXzV)
	    DO I=1,ATM(ID)%NXzV
	      J=ATM(ID)%EQXzV+I-1
	      TOT_ENERGY(J)=(AVE_ENERGY(ATM(ID)%EQXzV)-AVE_ENERGY(J))+T1
	      ION_EN(J)=T2
	    END DO
	    J=ATM(ID)%EQXzV
	    T1=T1+AVE_ENERGY(J)			!Adding on ionization energy
	  END DO
	  ID=SPECIES_END_ID(ISPEC)-1
	  IF(ID .GT. 0)THEN
	    J=ATM(ID)%EQXzV
	    TOT_ENERGY(J+ATM(ID)%NXzV)=T1
	    ION_EN(J+ATM(ID)%NXzV)=T2
	  END IF
	END DO
!
! Compute the mean energy per atom. At first it is units of 10^15Hz.
!
	INT_EN(:)=0.0D0
	COL_EN(:)=0.0D0
	DO I=1,ND
	  DO J=1,NT-2
	     INT_EN(I)=INT_EN(I)+POPS(J,I)*TOT_ENERGY(J)
	     COL_EN(I)=COL_EN(I)+POPS(J,I)*ION_EN(J)
	  END DO
	END DO
	INT_EN=HDKT*INT_EN/POP_ATOM
	COL_EN=HDKT*COL_EN/POP_ATOM
!
! We now compute constants for each of the 4 terms. These make
! it simpler and cleaner for the evaluation of the linearization.
!
! For historical reasons STEQ contains Int[chi.J - eta]dv. Rather than multiply
! this term everywhere by 4pi, we divide the adiabatic cooling rate by that
! term. Not also that the R units are unimportant, since R.Chi is dimensionless.
! The 10^9 arises since T is in units of 10^4K, and V in units of 10^5 km/s.
!
	PI=FUN_PI()
	SCALE=1.0D+09*BOLTZMANN_CONSTANT()/4.0D0/PI
	DO I=1,ND
	  IF(I .EQ. ND)THEN
	    T1=R(ND-1)-R(ND)
	  ELSE
	    T1=R(I)-R(I+1)
	  END IF
	  A(I)=1.5D0*SCALE*(POP_ATOM(I)+ED(I))*V(I)/T1
	  B(I)=SCALE*(POP_ATOM(I)+ED(I))*V(I)*(3.0D0+SIGMA(I))/R(I)
	  C(I)=1.5D0*SCALE*POP_ATOM(I)*V(I)/T1
	  D(I)=SCALE*POP_ATOM(I)*V(I)/T1
	  GAMMA(I)=ED(I)/POP_ATOM(I)
	END DO
!
	IF(INCL_ADIABATIC)THEN
	  DO I=1,ND-1
 	    WORK(I)=A(I)*(T(I)-T(I+1)) + B(I)*T(I) +
	1              C(I)*T(I)*(GAMMA(I)-GAMMA(I+1)) +
	1              D(I)*(INT_EN(I)-INT_EN(I+1))
	  END DO
 	  WORK(ND)=A(ND)*(T(ND-1)-T(ND)) + B(ND)*T(ND) +
	1              C(ND)*T(ND)*(GAMMA(ND-1)-GAMMA(ND)) +
	1              D(ND)*(INT_EN(ND-1)-INT_EN(ND))
!
	  DO I=1,ND
	    STEQ_T(I)=STEQ_T(I)-WORK(I)
	  END DO
	END IF
!
	IF(INCL_ADIABATIC .AND. COMPUTE_BA)THEN
	  DO I=1,ND-1
!
! Diagonal terms.
!
	    L=DIAG_INDX
	    BA_T(NT,L,I)=BA_T(NT,L,I)-A(I)-B(I)-C(I)*(GAMMA(I)-GAMMA(I+1))
	    BA_T(NT-1,L,I)=BA_T(NT-1,L,I)-(A(I)+B(I))/(POP_ATOM(I)+ED(I))-
	1                                  C(I)*T(I)/POP_ATOM(I)
	    DO J=1,NT-2
	      BA_T(J,L,I)=BA_T(J,L,I)-HDKT*D(I)*TOT_ENERGY(J)/POP_ATOM(I)
	    END DO
!
! Upper diagonal terms.
!
	    IF(NUM_BNDS .GE. 3)THEN
	      L=DIAG_INDX+1
	      BA_T(NT,L,I)=BA_T(NT,L,I)+A(I)
	      BA_T(NT-1,L,I)=BA_T(NT-1,L,I)+C(I)*T(I)/POP_ATOM(I+1)
	      DO J=1,NT-1
	        BA_T(J,L,I)=BA_T(J,L,I)+HDKT*D(I)*TOT_ENERGY(J)/POP_ATOM(I+1)
	      END DO
	    END IF
!
	  END DO	!Loop of depth.
!
! Need to do special case of I=ND
!
	  L=DIAG_INDX
	  BA_T(NT,L,ND)=BA_T(NT,L,ND)+A(ND)-B(ND)-C(ND)*(GAMMA(ND-1)-GAMMA(ND))
	  BA_T(NT-1,L,ND)=BA_T(NT-1,L,ND)-(A(ND)+B(ND))/(POP_ATOM(ND)+ED(ND)) +
	1                                  C(ND)*T(ND)/POP_ATOM(ND)
	  DO J=1,NT-2
	    BA_T(J,L,ND)=BA_T(J,L,ND)+HDKT*D(ND)*TOT_ENERGY(J)/POP_ATOM(ND)
	  END DO
!
	  IF(NUM_BNDS .GE. 3)THEN
	    L=DIAG_INDX-1
	    BA_T(NT,L,ND)=BA_T(NT,L,ND)-A(ND)
	    BA_T(NT-1,L,ND)=BA_T(NT-1,L,ND)-C(ND)*T(ND)/POP_ATOM(ND-1)
	    DO J=1,NT-2
	      BA_T(J,L,ND)=BA_T(J,L,ND)-HDKT*D(ND)*TOT_ENERGY(J)/POP_ATOM(ND-1)
	    END DO
	  END IF
!
	END IF            !End COMPUTE_BA
!
! Now compute the adiabatic cooling rate (in ergs/cm^3/sec) for diagnostic
! purposes. The rate is output to the COOLGEN file.
!
! The factor of 4.0D-10*PI arises from the fact that A, B etc were computed
! for the radiative equilibrium equation. That equation has units a factor of 10^10/4Pi
! larger (10^10 because of opacity definition, 4Pi as we don't scale J).
!
! We split the adiabatic terms into 2 parts: The velocity term, and the
! dTdR (internal energy) term. This split was useful for diagnostic purposes,
! but has now been kept for simplicity so that GENCOOL does not need to be changed.
!
	T1=4.0D-10*PI
	DO I=1,ND-1
	  AD_CR_V(I)= B(I)*T(I)
	  AD_CR_DT(I)=A(I)*(T(I)-T(I+1))+
	1               C(I)*T(I)*(GAMMA(I)-GAMMA(I+1)) +
	1               D(I)*(COL_EN(I)-COL_EN(I+1))
	END DO
	AD_CR_V(ND)=B(ND)*T(ND)
	AD_CR_DT(ND)=A(ND)*(T(ND-1)-T(ND))+
	1               C(ND)*T(ND)*(GAMMA(ND-1)-GAMMA(ND)) +
	1               D(ND)*(COL_EN(ND-1)-COL_EN(ND))
!
	AD_CR_V=AD_CR_V*T1
	AD_CR_DT=AD_CR_DT*T1
!
	WRITE_CHK=.TRUE.
	IF(WRITE_CHK)THEN
	  OPEN(UNIT=7,FILE='ADIABAT_CHK',STATUS='UNKNOWN')
	    WRITE(7,'(A)')' '
	    WRITE(7,'(A)')'  Scaling is for STEQ_VALS(NT,:). This is cgs units scaled'
	    WRITE(7,'(A)')'  by a factor of 10^10 on 4PI. The 10^9 comes from V . T'
	    WRITE(7,'(A)')' '
	    WRITE(7,'(A,ES12.4)')' SCALE=1.0D+09*BOLTZMANN_CONSTANT()/4.0D0/PI=',SCALE
	    WRITE(7,'(A)')' T1=R(I)-R(I+1)'
	    WRITE(7,'(A)')' GAMMA=ED/POP_ATOM'
	    WRITE(7,'(A)')' A=1.5D0*SCALE*(POP_ATOM+ED)*V/T1 * (T(I)-T(I+1)'
	    WRITE(7,'(A)')' B=SCALE*(POP_ATOM+ED)*V*(3.0D0+SIGMA)/R * T(I)'
	    WRITE(7,'(A)')' C=1.5D0*SCALE*POP_ATOM*V/T1 * T(I)*(GAMMA(I)-GAMMA(I+1))'
	    WRITE(7,'(A)')' D=SCALE*POP_ATOM*V/T1* (INT_EN(I)-INT_EN(I+1))'
	    WRITE(7,'(A)')' '
	    WRITE(7,'(8X,A,6X,8(7X,A))')'R','     V',' SIGMA','     T',
	1                 'NU_ION','     A','     B','     C','     D'
	    DO I=1,ND-1
	      WRITE(7,'(ES15.7,8(1X,ES12.4))')
	1              R(I),V(I),SIGMA(I),T(I),INT_EN(I)/HDKT,
	1              A(I)*(T(I)-T(I+1)),B(I)*T(I),
	1              C(I)*T(I)*(GAMMA(I)-GAMMA(I+1)),D(I)*(INT_EN(I)-INT_EN(I+1))
	    END DO
	    WRITE(7,'(ES15.7,8(1X,ES12.4))')
	1              R(ND),V(ND),SIGMA(ND),T(ND),INT_EN(ND)/HDKT,
	1              A(ND)*(T(ND-1)-T(ND)),B(ND)*T(ND),
	1              C(ND)*T(ND)*(GAMMA(ND-1)-GAMMA(ND)),D(ND)*(INT_EN(ND-1)-INT_EN(ND))
	  CLOSE(UNIT=7)
	END IF
!
	RETURN
	END
