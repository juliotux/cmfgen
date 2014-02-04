!
! Data module for FG_J_CMF. Data placed in this module is automatically
! saved between subroutine calls..
!
	MODULE FG_J_CMF_MOD_V10
	IMPLICIT NONE
!
! The *_STORE routines are used to store the radiation field, as computed
! on the previous call to FG_J_CMF. 
! The *_PREV routines are used to store the radiation field, as computed
! for the previous frequency.
!
! The *_PREV routines are updated from the *_STORE routines when NEW_FREQ
! is .TRUE.
!
!***************************************************************************
!***************************************************************************
!
! Data variables required by both routines.
!
! Dimensionded NRAY_MAX,NP
!
	REAL*8, ALLOCATABLE :: R_RAY(:,:)
	REAL*8, ALLOCATABLE :: Z(:,:)
	REAL*8, ALLOCATABLE :: GAM(:,:)
	REAL*8, ALLOCATABLE :: DTAU(:,:)
	REAL*8, ALLOCATABLE :: AV(:,:)
	REAL*8, ALLOCATABLE :: CV(:,:)
!
! Dimensiond ND_EXT,4
!
	REAL*8, ALLOCATABLE :: V_COEF(:,:)
	REAL*8, ALLOCATABLE :: SIGMA_COEF(:,:)
	REAL*8, ALLOCATABLE :: CHI_COEF(:,:)
	REAL*8, ALLOCATABLE :: ETA_COEF(:,:)
!
	INTEGER, ALLOCATABLE :: J_PNT(:,:)
	INTEGER, ALLOCATABLE :: H_PNT(:,:)
!
! Dimensioned ND+ND_ADD=ND_EXT
!
	REAL*8, ALLOCATABLE :: R_EXT(:)
	REAL*8, ALLOCATABLE :: LOG_R_EXT(:)
	REAL*8, ALLOCATABLE :: Z_EXT(:)
	REAL*8, ALLOCATABLE :: V_EXT(:)
	REAL*8, ALLOCATABLE :: SIGMA_EXT(:)
	REAL*8, ALLOCATABLE :: ETA_EXT(:)
	REAL*8, ALLOCATABLE :: CHI_EXT(:)
	REAL*8, ALLOCATABLE :: LOG_ETA_EXT(:)
	REAL*8, ALLOCATABLE :: LOG_CHI_EXT(:)
!
! Dimensioned NRAY_MAX,NP
!
	INTEGER, ALLOCATABLE :: RAY_PNT(:,:)
!                            
! Used to compute the boundary iteration factors (HBC and NBC).
!                             
	REAL*8, ALLOCATABLE :: N_WGHTS(:)
	REAL*8, ALLOCATABLE :: H_WGHTS(:)
	REAL*8, ALLOCATABLE :: MU_VAL(:)
	INTEGER, ALLOCATABLE :: NI_RAY(:)
	INTEGER, ALLOCATABLE :: MAX_LS(:)
!
! Dimensioned NRAY_MAX
!
	REAL*8, ALLOCATABLE :: dCHIdR(:)
	REAL*8, ALLOCATABLE :: Q(:)
	REAL*8, ALLOCATABLE :: QH(:)
	REAL*8, ALLOCATABLE :: V_RAY(:)
!
	REAL*8, ALLOCATABLE :: SIGMA_RAY(:)
	REAL*8, ALLOCATABLE :: ETA_RAY(:)
	REAL*8, ALLOCATABLE :: CHI_RAY(:)
	REAL*8, ALLOCATABLE :: dCHIdR_RAY(:)
	REAL*8, ALLOCATABLE :: SOURCE_RAY(:)
!
!***************************************************************************
!***************************************************************************
!
! Variables specific to the DIFFERENCE equation approach for  solving the
! transfer equation.
!	
	REAL*8, ALLOCATABLE :: AV_PREV(:,:)
	REAL*8, ALLOCATABLE :: AV_STORE(:,:)
	REAL*8, ALLOCATABLE :: CV_PREV(:,:)
	REAL*8, ALLOCATABLE :: CV_STORE(:,:)
	REAL*8, ALLOCATABLE :: PAR_AV(:,:)
!
	REAL*8, ALLOCATABLE :: GAMH(:,:)
	REAL*8, ALLOCATABLE :: TA(:,:)
	REAL*8, ALLOCATABLE :: TB(:,:)
	REAL*8, ALLOCATABLE :: TC(:,:)
	REAL*8, ALLOCATABLE :: GB(:,:)
	REAL*8, ALLOCATABLE :: H(:,:)
!
	REAL*8, ALLOCATABLE :: XM(:)
	REAL*8, ALLOCATABLE :: U(:)
	REAL*8, ALLOCATABLE :: VB(:)
	REAL*8, ALLOCATABLE :: VC(:)
	REAL*8, ALLOCATABLE :: DIV(:)
!
	REAL*8, ALLOCATABLE :: OLDCHI(:)
	REAL*8, ALLOCATABLE :: OLDCHI_STORE(:)
!
!***************************************************************************
!***************************************************************************
!
! Variables specific to short characteristic approach of solving the
! transfer equation.
!
! Dimensioned NRAY_MAX,NP
!
	REAL*8, ALLOCATABLE :: I_P_PREV(:,:)
	REAL*8, ALLOCATABLE :: I_P_STORE(:,:)
	REAL*8, ALLOCATABLE :: I_M_PREV(:,:)
	REAL*8, ALLOCATABLE :: I_M_STORE(:,:)
	REAL*8, ALLOCATABLE :: dGAMdR(:,:)
!
	REAL*8, ALLOCATABLE :: A0(:,:)
	REAL*8, ALLOCATABLE :: A1(:,:)
	REAL*8, ALLOCATABLE :: A2(:,:)
	REAL*8, ALLOCATABLE :: A3(:,:)
	REAL*8, ALLOCATABLE :: A4(:,:)
!
	REAL*8, ALLOCATABLE :: I_P(:,:)
	REAL*8, ALLOCATABLE :: I_M(:,:)
!
! Dimensioned NRAY_MAX
!
	REAL*8, ALLOCATABLE :: EE(:)
	REAL*8, ALLOCATABLE :: E0(:)
	REAL*8, ALLOCATABLE :: E1(:)
	REAL*8, ALLOCATABLE :: E2(:)
	REAL*8, ALLOCATABLE :: E3(:)
!
	REAL*8, ALLOCATABLE :: SOURCE_PRIME(:)
	REAL*8, ALLOCATABLE :: S(:)
	REAL*8, ALLOCATABLE :: dS(:)
!
	REAL*8 PREVIOUS_FREQ
	REAL*8 VDOP_FRAC_SAV
!
	INTEGER ND_EXT
	INTEGER ND_ADD
	INTEGER NP_MAX
	INTEGER NRAY_MAX
!
	CHARACTER*10 OLD_SOLUTION_OPTIONS
	CHARACTER*10 SOLUTION_METHOD
	LOGICAL INSERT
	LOGICAL FIRST_TIME
	LOGICAL NEW_R_GRID
	LOGICAL WRITE_IP
!
	REAL*8 HN_VAL
	REAL*8 HN_VAL_PREV
!
	DATA VDOP_FRAC_SAV/-1001.0D0/ 	!Absurd value.
	DATA FIRST_TIME/.TRUE./
	DATA PREVIOUS_FREQ/0.0D0/
!
	END MODULE FG_J_CMF_MOD_V10
!
! 
!
! Routine to compute the Eddington F, G and RSQN_ON_RSQJ Eddington factors 
! for a single frequency. The transfer is done in the comoving-frame with
! first order frequency differencing. A DIFFERENCE or INTEGRAL equation 
! approach can be used. 
!
!The DIFFERENCE equation approach is at least 30% faster. However,
! the INTEGRAL equation approach is more stable. Because monotonic cubic
! interpolation is used for the Source function in the INTEGRAL equation
! approach, the intensities are guaranteed positive (or zero).
! 
! The F, G, and RSQN_ON_RSQJ Eddington factors must be supplied.
!
! NB:
!	F = K / J
!	G=  N / H
!	RSQN_ON_RSQJ(I) = RSQ_N(I)/( RSQ_J(I)+ RQS_J(I+1))
!
! where
!	RSQ_N=0.25*N(I)*(R(I)+R(I+1))**2
!
! NB: Only one of G and RSQN_ON_RSQJ is defined at a given depth. This
!     avoids having to test what mode I am using for the Eddington factors.
!
!     IF N_TYPE='G_ONLY' G is defined at all depths, and
!       RSQN_ON_RSQJ=0 at all depths.
!     IF N_TYPE='N_ON_J' RSQN_ON_RSQJ is defined at all 
!       depths, and G=0 at all depths.
!     IF N_TYPE='MIXED' one of G or RSQN_ON_RSQJ is 
!       non-zero, and is the value to be used in MOM_J_CMF
!
! Routine also returns I+, so that observers flux can be computed. Note because
! the outer boundary can be thick, I+ is actually I+ - I- and hence is a flux
! like variable (This I+ = 2v at outer boundary).
!
! The radiation field at the previous (bluer) frequency is stored locally.
! Coupling to this frequency is controlled by the variable INIT. If INIT is
! true, the atmosphere is treated with the assumption that V=0.
!
! INIT must be TRUE on the very first call.
!
!
	SUBROUTINE FG_J_CMF_V10(ETA,CHI,ESEC,V,SIGMA,R,P,
	1                  JNU,HNU,KNU,NNU,RSQN_ON_RSQJ,
	1                  JQW,HQW,KQW,NQW,
	1                  IN_HBC,HBC,NBC,
	1                  IPLUS_P,FREQ,dLOG_NU,DIF,DBB,IC,
	1                  VDOP_VEC,VDOP_FRAC,
	1                  METHOD,SOLUTION_OPTIONS,
	1                  THK,INIT,NEW_FREQ,N_TYPE,NC,NP,ND)
	USE FG_J_CMF_MOD_V10
	IMPLICIT NONE
!
! Altered 16-Feb-2006  : Setting of RMAX altered so that it is not dependent on the
!                          number of grid points. Diferent setting for SN and winds.
! Altered 26-Feb-2004  : Changed to check for creating a new R_GRID etc.
!                          Done to handle automatic grid changes.
! Altered 13-Sep-2002  : Check before doing PSQ SQRT installed (for INTEL compiler).
! Altered 03-Feb-2002  : Bug fix. Can now treat THK=.FALSE. and NP=NC+ND.
!                          Routine was not handlling rays with NI=1 & 2 correctly.
!                          Altered RMAX for extension: Bow variable and not always 10.
!                          Needs further work.
! Altered 28-Feb-2002  : Extensively altered in early 2000:
!                        Minor bug fixes & cleaning.
!                        Routine now allows variable VDOP_FRAC. Useful for
!                          testing purposes. This facility should NOT be used with
!                          CMFGEN.
! Created 02-Nov-2000  : Based on FG_J_CMF_V9
!
	INTEGER NC,NP,ND
	REAL*8 ETA(ND),CHI(ND),ESEC(ND)
	REAL*8 V(ND),SIGMA(ND),R(ND),P(NP)
!
! NB: J,H,K,N refer to the first 4 moments of the radiation field.
!     QW denotes quadrature weight.
!
	REAL*8 JQW(ND,NP),HQW(ND,NP),KQW(ND,NP),NQW(ND,NP)
	REAL*8 JNU(ND),HNU(ND),KNU(ND),NNU(ND),RSQN_ON_RSQJ(ND)
	REAL*8 IN_HBC,HBC,NBC
	REAL*8 IPLUS_P(NP)
!
! VDOP_VEC(I) is the minimum DOPPLER width for all species at depth I. It will include
! both a turbulent, and and thermal contribution for the ionization species with the 
! highest mass. VDOP_FRAC is used to set the minimum velocity step size along a ray.
!
	REAL*8 VDOP_VEC(ND)
	REAL*8 VDOP_FRAC
!
	REAL*8 DBB,IC,FREQ,dLOG_NU
	CHARACTER*(*) SOLUTION_OPTIONS
	CHARACTER*6 METHOD
	CHARACTER*6 N_TYPE
	LOGICAL DIF		!Use diffusion approximation
!
! Use "Thick" boundary condition. at outer boundary. Only noted when INIT 
! is true. All subsequent frequencies will use the same boundary condition
! independent of the passed value (Until INIT is set to TRUE again).
!
	LOGICAL THK
!
! First frequency -- no frequency coupling.

	LOGICAL INIT
!
! Upon leaving this routine the radiation field along each ray is stored. This
! will provide the blue wing information necessary for the next frequency.
! This routine may, however, be used in an iterative loop. In this case the
! "blue wing" information should remain unaltered between calls.
! NEW_FREQ indicates that a new_frequency is being passed, and hence the "blue
! wing" information should be updated.
!
	LOGICAL NEW_FREQ	
!
	INTEGER ND_ADD_MAX
	PARAMETER (ND_ADD_MAX=24)
!
	INTEGER, SAVE :: FREQ_CNT
	INTEGER LU_IP
	INTEGER ACCESS_F
        INTEGER IOS,REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC
!
! The following arrays do not need to be stored, and hence can be created 
! each time.
!
	REAL*8 CV_BOUND(NP)		!Outer boundary V
	REAL*8 I_M_IN_BND(NP)		!Inner boundary
	REAL*8 IBOUND(NP)		!Incident intensity on outer boundary.
!
	INTEGER N_ERR_MAX,FG_ERR_CNT
	PARAMETER (N_ERR_MAX=1000)
	REAL*8 FG_ERR_ON_FREQ
	INTEGER FG_ERR_TYPE
	COMMON /FG_J_CMF_ERR/FG_ERR_ON_FREQ(N_ERR_MAX),
	1                    FG_ERR_TYPE(N_ERR_MAX),FG_ERR_CNT
	LOGICAL NEG_AV_VALUE,NEG_NNU_VALUE,BAD_NNU_VALUE
!
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
! Local variables.
!
	REAL*8, PARAMETER :: ONE=1
	INTEGER, PARAMETER :: NINS=4
!
	LOGICAL, PARAMETER :: LFALSE=.FALSE.
	LOGICAL, PARAMETER :: LTRUE=.TRUE.
!                
	INTEGER NI_SMALL
	INTEGER I,J,K,LS
	INTEGER IDMIN,IDMAX
	INTEGER NI
!
	REAL*8 DBC
	REAL*8 I_CORE
	REAL*8 T1,T2
	REAL*8 DELTA_Z
	REAL*8 ALPHA
	REAL*8 ESEC_POW
	REAL*8 BETA
	REAL*8 VINF
	REAL*8 RMAX,DEL_R_FAC
	REAL*8 MU,dZ,PSQ
	REAL*8 DEL_R
!
! Change the following statement to TRUE if running on a VECTOR machine.
!
	LOGICAL, PARAMETER :: VECTOR_MACHINE=.FALSE.
!
!
!
! Check to see whether we have a new R grid, or the solution options have
! changd. This can only happen when INIT is TRUE.
!
	NEW_R_GRID=.FALSE.
	IF(INIT .AND. .NOT. FIRST_TIME)THEN
	  DO I=1,ND
	    IF(R(I) .NE. R_EXT(ND_ADD+I))THEN
	      NEW_R_GRID=.TRUE.
	      WRITE(171,*)'Updating RGRID in FG_J_CMF_V10'
	      EXIT
	    END IF
	  END DO
	  IF(VDOP_FRAC .NE. VDOP_FRAC_SAV
	1       .OR.  OLD_SOLUTION_OPTIONS .NE. SOLUTION_OPTIONS)NEW_R_GRID=.TRUE.
	END IF
	
! Deallocate all allocated rays if we are using a diferent solution technique.
! This option will only be used when testing, since in CMFGEN we will always use
! the same atmospheric structure.
!
	IF( ALLOCATED(R_EXT) .AND. NEW_R_GRID)THEN
	  DEALLOCATE ( R_EXT )
	  DEALLOCATE ( LOG_R_EXT )
	  DEALLOCATE ( V_EXT )
	  DEALLOCATE ( Z_EXT )
	  DEALLOCATE ( SIGMA_EXT )
	  DEALLOCATE ( ETA_EXT )
	  DEALLOCATE ( CHI_EXT )
	  DEALLOCATE ( LOG_ETA_EXT )
	  DEALLOCATE ( LOG_CHI_EXT )
!
	  DEALLOCATE ( Z )
	  DEALLOCATE ( RAY_PNT )
	  DEALLOCATE ( R_RAY )
	  DEALLOCATE ( NI_RAY )
	  DEALLOCATE ( J_PNT )
	  DEALLOCATE ( H_PNT )
	  DEALLOCATE ( N_WGHTS )
	  DEALLOCATE ( H_WGHTS )
	  DEALLOCATE ( MU_VAL )
	  DEALLOCATE ( V_COEF )
	  DEALLOCATE ( SIGMA_COEF )
	  DEALLOCATE ( ETA_COEF )
	  DEALLOCATE ( CHI_COEF )
	  DEALLOCATE ( V_RAY )
	  DEALLOCATE ( SIGMA_RAY )
	  DEALLOCATE ( ETA_RAY )
	  DEALLOCATE ( CHI_RAY )
	  DEALLOCATE ( dCHIdR_RAY )
	  DEALLOCATE ( SOURCE_RAY )
	  DEALLOCATE ( dCHIdR )
	  DEALLOCATE ( Q )
	  DEALLOCATE ( MAX_LS )
	  DEALLOCATE ( AV )
	  DEALLOCATE ( CV )
	  DEALLOCATE ( DTAU )
	  DEALLOCATE ( GAM )
!
! These arrays are only used when using the DIFFERENCE apporach.
!
	  IF( ALLOCATED(AV_PREV) )THEN
	    DEALLOCATE ( AV_PREV )
	    DEALLOCATE ( AV_STORE )
	    DEALLOCATE ( PAR_AV )
	    DEALLOCATE ( CV_PREV )
	    DEALLOCATE ( CV_STORE )
	    DEALLOCATE ( GAMH )
!
	    DEALLOCATE ( TA )
	    DEALLOCATE ( TB )
	    DEALLOCATE ( TC )
	    DEALLOCATE ( GB )
	    DEALLOCATE ( H )
	    DEALLOCATE ( OLDCHI )
	    DEALLOCATE ( OLDCHI_STORE )
!
	    DEALLOCATE ( XM )
	    DEALLOCATE ( U )
	    DEALLOCATE ( VB )
	    DEALLOCATE ( VC )
	    DEALLOCATE ( QH )
	    DEALLOCATE ( DIV )
	  END IF
!
! These arrays are only used when using the INTEGERAL apporach.
!
	  IF( ALLOCATED(I_P) )THEN
	    DEALLOCATE ( I_P )
	    DEALLOCATE ( I_M )
!
	    DEALLOCATE ( I_P_PREV )
	    DEALLOCATE ( I_P_STORE )
	    DEALLOCATE ( I_M_PREV )
	    DEALLOCATE ( I_M_STORE )
	    DEALLOCATE ( dGAMdR )
!
	    DEALLOCATE ( A0 )
	    DEALLOCATE ( A1 )
	    DEALLOCATE ( A2 )
	    DEALLOCATE ( A3 )
	    DEALLOCATE ( A4 )
!
	    DEALLOCATE (EE)
	    DEALLOCATE (E0)
	    DEALLOCATE (E1)
	    DEALLOCATE (E2)
	    DEALLOCATE (E3)
!
	    DEALLOCATE (SOURCE_PRIME)
	    DEALLOCATE (S)
	    DEALLOCATE (dS)
	  END IF
	END IF
	VDOP_FRAC_SAV=VDOP_FRAC
!
! 
!
! Set up the revised grid to improve computational accuracy. Unles we are
! carrying out tests, these need only be constructed once.
!
	IF(FIRST_TIME .OR. .NOT. ALLOCATED(R_EXT) )THEN
!
	  OLD_SOLUTION_OPTIONS=SOLUTION_OPTIONS
	  IF(SOLUTION_OPTIONS .EQ. 'DIFF/INS')THEN
	    SOLUTION_METHOD='DIFFERENCE'
	    INSERT=.TRUE.
	  ELSE IF(SOLUTION_OPTIONS(1:4) .EQ. 'DIFF')THEN
	    SOLUTION_METHOD='DIFFERENCE'
	    INSERT=.FALSE.
	  ELSE IF(SOLUTION_OPTIONS .EQ. 'INT/INS')THEN
	    SOLUTION_METHOD='INTEGRAL'
	    INSERT=.TRUE.
	  ELSE IF(SOLUTION_OPTIONS(1:3) .EQ. 'INT')THEN
	    SOLUTION_METHOD='INTEGRAL'
	    INSERT=.FALSE.
	  ELSE
	    J=ERROR_LU()
	    WRITE(J,*)'Error in FG_J_CMF_V10: Invalid solution option'
	    WRITE(J,*)SOLUTION_OPTIONS
	    STOP
	  END IF
!
          ND_ADD=0
          IF(THK)ND_ADD=ND_ADD_MAX
          ND_EXT=ND+ND_ADD
!
	  ALLOCATE ( R_EXT(ND_EXT) )
	  ALLOCATE ( LOG_R_EXT(ND_EXT) )
	  ALLOCATE ( V_EXT(ND_EXT) )
	  ALLOCATE ( Z_EXT(ND_EXT) )
	  ALLOCATE ( SIGMA_EXT(ND_EXT) )
	  ALLOCATE ( ETA_EXT(ND_EXT) )
	  ALLOCATE ( CHI_EXT(ND_EXT) )
	  ALLOCATE ( LOG_ETA_EXT(ND_EXT) )
	  ALLOCATE ( LOG_CHI_EXT(ND_EXT) )
!
! Compute the extended R grid, excluding inserted points.
!
	  DO I=1,ND
	    R_EXT(ND_ADD+I)=R(I)
	  END DO
	  IF(THK)THEN
	    IF(V(ND) .LT. 10.0D0 .AND. R(1)/R(ND) .GE. 9.99D0)THEN
	      RMAX=10.0D0*R(1)		!Stellar wind
	    ELSE IF(V(ND) .GT. 10.0D0 .OR. V(1) .GT. 2.0D+04)THEN
	      RMAX=1.5D0*R(1)		!SN model
	    ELSE
	      RMAX=MIN(10.0D0,SQRT(R(1)/R(ND)))*R(1)
	    END IF
	    ALPHA=R(1)+(R(1)-R(2))                               
	    DEL_R_FAC=EXP( LOG(RMAX/ALPHA)/(ND_ADD-3) )
	    R_EXT(1)=RMAX
	    R_EXT(4)=RMAX/DEL_R_FAC
	    R_EXT(2)=R_EXT(1)-0.1*(R_EXT(1)-R_EXT(4))
	    R_EXT(3)=R_EXT(1)-0.4*(R_EXT(1)-R_EXT(4))
	    DO I=5,ND_ADD-1
	      R_EXT(I)=R_EXT(I-1)/DEL_R_FAC
	    END DO
	    R_EXT(ND_ADD)=ALPHA
	    J=ERROR_LU()
	  END IF
!
! Compute VEXT and R_EXT. We assume a BETA velocity law at large R.
!
	  V_EXT(ND_ADD+1:ND_EXT)=V(1:ND)
	  SIGMA_EXT(ND_ADD+1:ND_EXT)=SIGMA(1:ND)
	  IF(THK)THEN
	    BETA=(SIGMA(1)+1.0D0)*(R(1)/R(ND)-1.0D0)
            VINF=V(1)/(1-R(ND)/R(1))**BETA
	    DO I=1,ND_ADD
	      V_EXT(I)=VINF*(1.0D0-R_EXT(ND_EXT)/R_EXT(I))**BETA
	      SIGMA_EXT(I)=BETA/(R_EXT(I)/R_EXT(ND_EXT)-1.0D0)-1.0D0
	    END DO
	    WRITE(J,'(A)')' '
	    WRITE(J,*)'Using thick boundary condition in FG_J_CMF_V10'
	    WRITE(J,'(2(A,ES16.8,3X))')' R(1)=',R(1),'RMAX=',RMAX
	    WRITE(J,'(2(A,ES16.8,3X))')' V(1)=',V(1),'VMAX=',V_EXT(1)
	  END IF
!
! Work out the maximim number of points per ray. This will allow us to allocate
! the require memory.
!
	  NRAY_MAX=0.0
	  DO LS=1,NP
	    NI_SMALL=ND_EXT-(LS-NC-1)
	    IF(LS .LE. NC+1)NI_SMALL=ND_EXT
!
	    DO I=1,NI_SMALL
	      IF(R_EXT(I) .EQ. P(LS))THEN
	        Z_EXT(I)=0.0D0
	      ELSE
	        Z_EXT(I)=SQRT(R_EXT(I)*R_EXT(I)-P(LS)*P(LS))
	      END IF
	   END DO
!
	    K=1
	    T2=VDOP_FRAC*MINVAL(VDOP_VEC(1:NI_SMALL-ND_ADD)) 
	    DO I=1,NI_SMALL-1
	      T1=(Z_EXT(I)*V_EXT(I)/R_EXT(I)-Z_EXT(I+1)*V_EXT(I+1)/R_EXT(I+1))/T2
	      IF(.NOT. THK .AND. NI_SMALL .EQ. 2)T1=MAX(2.01D0,T1)
	      IF(T1 .GT. 1)K=K+INT(T1)
	      K=K+1
	    END DO
            NRAY_MAX=MAX(NRAY_MAX,K)
	  END DO
	  J=ERROR_LU()
	  WRITE(J,'(A,I6)')'Maximum number of pints along ray in FG_J_CMF_V10 is',NRAY_MAX
!
	  ALLOCATE ( Z(NRAY_MAX,NP) )
	  ALLOCATE ( RAY_PNT(NRAY_MAX,NP) )
	  ALLOCATE ( R_RAY(NRAY_MAX,NP) )
	  ALLOCATE ( NI_RAY(NP) )
!
! 
!
! Define the grid along each ray. NB: RAY_PNT is used to indicate the
! interpolation interval on the original (EXT) grid.
!
	  DO LS=1,NP
	    NI_SMALL=ND_EXT-(LS-NC-1)
	    IF(LS .LE. NC+1)NI_SMALL=ND_EXT
!
	    DO I=1,NI_SMALL
	      IF(R_EXT(I) .EQ. P(LS))THEN
	        Z_EXT(I)=0.0D0
	      ELSE
	        Z_EXT(I)=SQRT(R_EXT(I)*R_EXT(I)-P(LS)*P(LS))
	      END IF
	    END DO
!
	    K=1
	    Z(1,LS)=Z_EXT(1)
	    RAY_PNT(1,LS)=1
	    T2=VDOP_FRAC*MINVAL(VDOP_VEC(1:NI_SMALL-ND_ADD)) 
	    DO I=1,NI_SMALL-1
	      T1=(Z_EXT(I)*V_EXT(I)/R_EXT(I)-Z_EXT(I+1)*V_EXT(I+1)/R_EXT(I+1))/T2
	      IF(.NOT. THK .AND. NI_SMALL .EQ. 2)T1=MAX(2.01D0,T1)
	      IF(T1 .GT. 1)THEN
	        DELTA_Z=(Z_EXT(I+1)-Z_EXT(I))/(INT(T1)+1)
	        DO J=1,INT(T1)
	          K=K+1
	          Z(K,LS)=Z(K-1,LS)+DELTA_Z
	          RAY_PNT(K,LS)=I
	        END DO
	      END IF
	      K=K+1
	      Z(K,LS)=Z_EXT(I+1)
	      RAY_PNT(K,LS)=I
	    END DO
	    NI_RAY(LS)=K
!
	    PSQ=P(LS)*P(LS)
	    DO I=1,NI_RAY(LS)
	      R_RAY(I,LS)=SQRT(Z(I,LS)*Z(I,LS)+PSQ)
	    END DO
	    R_RAY(NI_RAY(LS),LS)=R_EXT(NI_SMALL)
!
	  END DO
!
! 
!
! J_PNT and H_PNT are used to indicate the postion of J(I) amd H(I) along the ray
! so that J and H can be computed on our regular radius grid.
!
	  ALLOCATE ( J_PNT(ND,NP) ); J_PNT(:,:)=0
	  ALLOCATE ( H_PNT(ND,NP) ); H_PNT(:,:)=0
!
	  DO LS=1,NP
	    NI_SMALL=ND-(LS-NC-1);   IF(LS .LE. NC+1)NI_SMALL=ND
	    K=1
	    DO I=1,NI_SMALL
	      DO WHILE(J_PNT(I,LS) .EQ. 0)
	        IF(R(I) .LE. R_RAY(K,LS) .AND. R(I) .GE. R_RAY(K+1,LS))THEN
	          IF( (R_RAY(K,LS)-R(I)) .LT. (R(I)-R_RAY(K+1,LS)) )THEN
	            J_PNT(I,LS)=K
	          ELSE
	            J_PNT(I,LS)=K+1
	          END IF
	        ELSE
	          K=K+1
	        END IF
	      END DO
	    END DO
	  END DO
!
	  IF( SOLUTION_METHOD .EQ. 'INTEGRAL')THEN
!
! CV is defined on the nodes. Need to interpolate to the midpoint of R.
!
	    DO LS=1,NP
	      NI_SMALL=ND-(LS-NC-1);   IF(LS .LE. NC+1)NI_SMALL=ND
	      K=1
	      DO I=1,NI_SMALL-1
	        DO WHILE(H_PNT(I,LS) .EQ. 0)
	          T1=0.5D0*(R(I)+R(I+1))
	          IF(T1 .LE. R_RAY(K,LS) .AND. T1 .GE. R_RAY(K+1,LS))THEN
	            H_PNT(I,LS)=K
	          ELSE
	            K=K+1
	          END IF
	        END DO
	      END DO
	    END DO
	  ELSE
!
! CV is defined at the mid-points, but these don't, in general, correspond to
! the midpoints of the R grid.
!
	    DO LS=1,NP
	      K=1
	      DO I=1,ND-1
	        DO WHILE(H_PNT(I,LS) .EQ. 0)
	          T1=0.5D0*(R(I)+R(I+1))
	          IF(T1 .LE. 0.5D0*(R_RAY(K,LS)+R_RAY(K+1,LS)) .AND. T1 .GE. 0.5D0*(R_RAY(K+1,LS)+R_RAY(K+2,LS)) )THEN
	            H_PNT(I,LS)=K
	          ELSE
	            K=K+1
	          END IF
	        END DO
	      END DO
	    END DO
	  END IF
!
! 
!
	  ALLOCATE ( N_WGHTS(NP) )
	  ALLOCATE ( H_WGHTS(NP) )
	  ALLOCATE ( MU_VAL(NP) )
!
	  ALLOCATE ( V_COEF(ND_EXT,4) )
	  ALLOCATE ( SIGMA_COEF(ND_EXT,4) )
	  ALLOCATE ( ETA_COEF(ND_EXT,4) )
	  ALLOCATE ( CHI_COEF(ND_EXT,4) )
!
!***************************************************************************
!***************************************************************************
!
	  ALLOCATE ( V_RAY(NRAY_MAX) )
	  ALLOCATE ( SIGMA_RAY(NRAY_MAX) )
	  ALLOCATE ( ETA_RAY(NRAY_MAX) )
	  ALLOCATE ( CHI_RAY(NRAY_MAX) )
	  ALLOCATE ( dCHIdR_RAY(NRAY_MAX) )
	  ALLOCATE ( SOURCE_RAY(NRAY_MAX) )
	  ALLOCATE ( dCHIdR(NRAY_MAX) )
	  ALLOCATE ( Q(NRAY_MAX))
	  ALLOCATE ( MAX_LS(NRAY_MAX) )
!
	  ALLOCATE ( AV(NRAY_MAX,NP) )
	  ALLOCATE ( CV(NRAY_MAX,NP) )
	  ALLOCATE ( DTAU(NRAY_MAX,NP) )
	  ALLOCATE ( GAM(NRAY_MAX,NP) )
!
!***************************************************************************
!***************************************************************************
!
! Allocation of variables required by both DIFFERENCE solution approach.
!
	  IF(SOLUTION_METHOD .EQ. 'DIFFERENCE')THEN
	    ALLOCATE ( AV_PREV(NRAY_MAX,NP) )
	    ALLOCATE ( AV_STORE(NRAY_MAX,NP) )
	    ALLOCATE ( PAR_AV(NRAY_MAX,NP) )
	    ALLOCATE ( CV_PREV(NRAY_MAX,NP) )
	    ALLOCATE ( CV_STORE(NRAY_MAX,NP) )
	    ALLOCATE ( GAMH(NRAY_MAX,NP) )
!
	    ALLOCATE ( TA(NRAY_MAX,NP) )
	    ALLOCATE ( TB(NRAY_MAX,NP) )
	    ALLOCATE ( TC(NRAY_MAX,NP) )
	    ALLOCATE ( GB(NRAY_MAX,NP) )
	    ALLOCATE ( H(NRAY_MAX,NP) )
	    ALLOCATE ( OLDCHI(NP) )
	    ALLOCATE ( OLDCHI_STORE(NP) )
!
	    ALLOCATE ( XM(NRAY_MAX) )
	    ALLOCATE ( U(NRAY_MAX) )
	    ALLOCATE ( VB(NRAY_MAX) )
	    ALLOCATE ( VC(NRAY_MAX) )
	    ALLOCATE ( QH(NRAY_MAX) )
	    ALLOCATE ( DIV(NRAY_MAX) )
	  END IF
!
!***************************************************************************
!***************************************************************************
!
! Allocation of variables required by INTEGRAL solution.
!
	  IF(SOLUTION_METHOD .EQ. 'INTEGRAL')THEN
!
	    ALLOCATE (I_P(NRAY_MAX,NP))
	    ALLOCATE (I_M(NRAY_MAX,NP))
!
	    ALLOCATE ( I_P_PREV(NRAY_MAX,NP) )
	    ALLOCATE ( I_P_STORE(NRAY_MAX,NP) )
	    ALLOCATE ( I_M_PREV(NRAY_MAX,NP) )
	    ALLOCATE ( I_M_STORE(NRAY_MAX,NP) )
	    ALLOCATE ( dGAMdR(NRAY_MAX,NP) )
!
	    ALLOCATE ( A0(NRAY_MAX,NP) )
	    ALLOCATE ( A1(NRAY_MAX,NP) )
	    ALLOCATE ( A2(NRAY_MAX,NP) )
	    ALLOCATE ( A3(NRAY_MAX,NP) )
	    ALLOCATE ( A4(NRAY_MAX,NP) )
!
	    ALLOCATE (EE(NRAY_MAX))
	    ALLOCATE (E0(NRAY_MAX))
	    ALLOCATE (E1(NRAY_MAX))
	    ALLOCATE (E2(NRAY_MAX))
	    ALLOCATE (E3(NRAY_MAX))
!
	    ALLOCATE (SOURCE_PRIME(NRAY_MAX))
	    ALLOCATE (S(NRAY_MAX))
	    ALLOCATE (dS(NRAY_MAX))
!
	  END IF
!
	  FIRST_TIME=.FALSE.
	ELSE
	  IF(OLD_SOLUTION_OPTIONS .NE. SOLUTION_OPTIONS)THEN
	    J=ERROR_LU()
	    WRITE(J,*)'Error in FG_J_CMF_V10'
	    WRITE(J,*)'Can''t switch SOLUTION_OPTIONS while runing code'
	    WRITE(J,*)'New setting:',SOLUTION_OPTIONS
	    WRITE(J,*)'Old setting:',OLD_SOLUTION_OPTIONS
	    STOP
	  END IF
	END IF
!
!
!
	NEG_AV_VALUE=.FALSE.
	NEG_NNU_VALUE=.FALSE.
	BAD_NNU_VALUE=.FALSE.
!         
! Perform initializations. 
!
	IF(INIT)THEN
!
! Insert extra points into radius grid. Not all points will be used along a ray.
! We only insert additional points in the interval between Z(NI)=0 and Z(NI-1),
! and between Z(NI-1) and  Z(NI-2).
!
	  FG_ERR_ON_FREQ(:)=0.0D0
	  FG_ERR_TYPE(:)=0
	  FG _ERR_CNT=0
	  IF(SOLUTION_METHOD .EQ. 'DIFFERENCE')THEN
	    AV_PREV(:,:)=0.0D0
	    CV_PREV(:,:)=0.0D0
          ELSE IF(SOLUTION_METHOD .EQ. 'INTEGRAL')THEN
	    I_P_PREV(:,:)=0.0D0
	    I_M_PREV(:,:)=0.0D0
	  END IF
!
	  DO LS=1,NP
	    IF(R(1) .EQ. P(LS))THEN
	      MU_VAL(LS)=0.0D0
	    ELSE
	      MU_VAL(LS)=SQRT(R(1)*R(1)-P(LS)*P(LS))/R(1)
	    END IF
	    IPLUS_P(LS)=0.0D0
	  END DO
	  CALL HTRPWGT(MU_VAL,H_WGHTS,NP)
	  CALL NTRPWGT(MU_VAL,N_WGHTS,NP)
!
	  LOG_R_EXT(1:ND_EXT)=LOG(R_EXT(1:ND_EXT))
	  CALL MON_INT_FUNS_V2(V_COEF,V_EXT,LOG_R_EXT,ND_EXT)
	  CALL MON_INT_FUNS_V2(SIGMA_COEF,SIGMA_EXT,LOG_R_EXT,ND_EXT)
!
	  DO LS=1,NP
!
	    DO I=1,NI_RAY(LS)
	      K=RAY_PNT(I,LS)
	      T1=LOG(R_RAY(I,LS)/R_EXT(K))
	      V_RAY(I)=((V_COEF(K,1)*T1+V_COEF(K,2))*T1+V_COEF(K,3))*T1+V_COEF(K,4)
	      SIGMA_RAY(I)=((SIGMA_COEF(K,1)*T1+SIGMA_COEF(K,2))*T1+SIGMA_COEF(K,3))*T1+SIGMA_COEF(K,4)
	    END DO
!
!
! Compute GAMMA. This section is straight from the subroutine GAMMA, except
! That _EXT has been added to V, SIGMA, and R.
!
! We assume (1)	SIGMAd+1/2 = 0.5*( SIGMAd+1+SIGMAd )
!  	    (2)	Vd+1/2=0.5*( Vd + Vd+1 )
! Note that V is in km/s and SIGMA=(dlnV/dlnR-1.0)
!
	    NI=NI_RAY(LS)
	    DO I=1,NI_RAY(LS)
	      MU=Z(I,LS)/R_RAY(I,LS)
	      T1=3.33564D-06*V_RAY(I)/R_RAY(I,LS)
	      GAM(I,LS)=T1*( 1.0D0+SIGMA_RAY(I)*(MU**2) )
	    END DO
!
	    IF(SOLUTION_METHOD .EQ. 'INTEGRAL' .AND. NI_RAY(LS) .NE. 1)THEN
	      DO I=1,NI_RAY(LS)
	        MU=Z(I,LS)/R_RAY(I,LS)
	        T1=3.33564D-06*V_RAY(I)/R_RAY(I,LS)
	        dGAMdR(I,LS)=GAM(I,LS)*SIGMA_RAY(I)/R_RAY(I,LS)
	        J=MAX(I-1,1); K=MIN(NI,I+1)
	        dGAMdR(I,LS)=dGAMdR(I,LS)+T1*(MU**2)*
	1           (SIGMA_RAY(J)-SIGMA_RAY(K))/(R_RAY(J,LS)-R_RAY(K,LS))
	      END DO
	    ELSE IF(SOLUTION_METHOD .EQ. 'DIFFERENCE')THEN
	      DO I=1,NI_RAY(LS)-1
	        MU=(Z(I,LS)+Z(I+1,LS))/(R_RAY(I,LS)+R_RAY(I+1,LS))
	        GAMH(I,LS)=(V_RAY(I+1)+V_RAY(I))*3.33564D-06*
	1                   (  1.0D0+0.5D0*(MU**2)*
	1                   (SIGMA_RAY(I)+SIGMA_RAY(I+1))  )/
	1                   (R_RAY(I,LS)+R_RAY(I+1,LS))
	     END DO
	     GAMH(NI,LS)=0.0
	    END IF
!
	  END DO		!LS Loop
!
! Determine maximum ray index having I points. Used in the
! vector section of the DIFFERENCE approach.
!
	  DO I=1,NRAY_MAX
	    DO LS=1,NP
	      IF(NI_RAY(LS) .GE. I)MAX_LS(I)=LS
	    END DO
	  END DO
	  HN_VAL=0.0D0
!
	ELSE IF(NEW_FREQ)THEN
          IF(SOLUTION_METHOD .EQ. 'DIFFERENCE')THEN
	    OLDCHI(1:NP)=OLDCHI_STORE(1:NP)
	    AV_PREV(:,:)=AV_STORE(:,:)
	    CV_PREV(:,:)=CV_STORE(:,:)
	  ELSE IF(SOLUTION_METHOD .EQ. 'INTEGRAL')THEN
	    I_P_PREV(:,:)=I_P_STORE(:,:)
	    I_M_PREV(:,:)=I_M_STORE(:,:)
	  END IF
	  IF(FREQ .GE. PREVIOUS_FREQ)THEN
	    WRITE(ERROR_LU(),*)'Error in FG_J_CMF_V10'
	    WRITE(ERROR_LU(),*)'Frequencies must be monotonically decreasng'
	    STOP
	  END IF
	  HN_VAL_PREV=HN_VAL
	ELSE IF(FREQ .NE. PREVIOUS_FREQ)THEN
	   WRITE(ERROR_LU(),*)'Error in FG_J_CMF_V10'
	   WRITE(ERROR_LU(),*)
	1      'Frequencies must not change if NEW_FREQ=.FALSE.'
	   STOP
	END IF
	PREVIOUS_FREQ=FREQ
!
! 
!
	CALL TUNE(1,'FG_BIG_SEC')
!
! Compute CHI_EXT, and ETA_EXT. CHI_EXT could be saved, as it doesn't change
! during the iteration procedure. ETA does however change (since it depends
! of J) and thus ETA_EXT must be re-computed on each entry.
!
! The first checks whether we may have negative line opacities due
! to stimulated emission. In such a case we simply assume an 1/r^2
! extrapolation.
!
! We also interpolate in ESEC, since ESEC (in the absence of negative 
! absorption) provides a lower bound to the opacity. NB: When CHI is much
! larger then ESEC its variation with r dominates, and it is possible to
! extrapolate CHI below ESEC.
!
	CALL TUNE(1,'FG_CHI_BEG')
	IF(ND_ADD .NE. 0)THEN
	  IDMIN=1; IDMAX=4
	  IF(CHI(IDMIN) .LE. ESEC(IDMIN) .OR. CHI(IDMAX) .LE. ESEC(IDMAX))THEN
	    ESEC_POW=LOG(ESEC(IDMAX)/ESEC(IDMIN))/LOG(R(IDMIN)/R(IDMAX))
	    IF(ESEC_POW .LT. 2.)ESEC_POW=2
	    DO I=1,ND_ADD
	      CHI_EXT(I)=CHI(IDMIN)*(R(IDMIN)/R_EXT(I))**ESEC_POW
	      dCHIdR(I)= -ESEC_POW*CHI_EXT(I)/R_EXT(I)
	    END DO
	  ELSE
	    ALPHA=LOG( (CHI(IDMAX)-ESEC(IDMAX)) / (CHI(IDMIN)-ESEC(IDMIN)) )
	1          /LOG(R(IDMIN)/R(IDMAX))
	    IF(ALPHA .LT. 2.)ALPHA=2.0
	    ESEC_POW=LOG(ESEC(IDMAX)/ESEC(IDMIN))/LOG(R(IDMIN)/R(IDMAX))
	    IF(ESEC_POW .LT. 2.)ESEC_POW=2
   	    DO I=1,ND_ADD
	      T1=(CHI(IDMIN)-ESEC(IDMIN))*(R(IDMIN)/R_EXT(I))**ALPHA
	      T2=ESEC(IDMIN)*(R(IDMIN)/R_EXT(I))**ESEC_POW
	      CHI_EXT(I)=T1+T2
	      dCHIdR(I)=(-ALPHA*T1-ESEC_POW*T2)/R_EXT(I)
	    END DO
	  END IF
	  DO I=ND_ADD+1,ND_EXT
	    CHI_EXT(I)=CHI(I-ND_ADD)
	  END DO
!
! We limit alpha to 3.5 to avoid excess envelope emission. If alpha were
! 3 we would get a logarithmic flux divergence as we increase the volume.
!
	  ALPHA=LOG(ETA(4)/ETA(1))/LOG(R(1)/R(4))
	  IF(ALPHA .LT. 3.5)ALPHA=3.5
	  DO I=1,ND_ADD
	    ETA_EXT(I)=ETA(1)*(R(1)/R_EXT(I))**ALPHA
	    IF(ETA_EXT(I) .LE. 1.0D-280)ETA_EXT(I)=1.0D-280
	  END DO
	  DO I=ND_ADD+1,ND_EXT
	    ETA_EXT(I)=ETA(I-ND_ADD)
	  END DO
	ELSE
	  CHI_EXT(1:ND)=CHI(1:ND)		!NB: In this can ND=ND_EXT
	  ETA_EXT(1:ND)=ETA(1:ND)
	END IF
!
	LOG_CHI_EXT(1:ND_EXT)=LOG(CHI_EXT(1:ND_EXT))
	CALL MON_INT_FUNS_V2(CHI_COEF,LOG_CHI_EXT,LOG_R_EXT,ND_EXT)
	LOG_ETA_EXT(1:ND_EXT)=LOG(ETA_EXT(1:ND_EXT))
	CALL MON_INT_FUNS_V2(ETA_COEF,LOG_ETA_EXT,LOG_R_EXT,ND_EXT)
	CALL TUNE(2,'FG_CHI_BEG')
!
! 
!***************************************************************************
!***************************************************************************
!
	IF(SOLUTION_METHOD .EQ. 'DIFFERENCE')THEN
!
! Zero AV and CV matrices.
!
	  AV(:,:)=0.0D0
	  CV(:,:)=0.0D0
!
!	  CALL TUNE(1,'LS_LOOP')
!
	  DO LS=1,NP
	    NI=NI_RAY(LS)
!
! NB: dCHIdR = dLOG(CHI)/dLOG(R) * CHI/R
!
	    DO I=1,NI_RAY(LS)
	      K=RAY_PNT(I,LS)
	      IF(R_RAY(I,LS) .EQ. R_EXT(K))THEN
	        CHI_RAY(I)=CHI_EXT(K)
	        ETA_RAY(I)=ETA_EXT(K)
	        dCHIdR_RAY(I)=CHI_COEF(K,3)*CHI_RAY(I)/R_RAY(I,LS)
	      ELSE
	        T1=LOG(R_RAY(I,LS)/R_EXT(K))
	        T2=((CHI_COEF(K,1)*T1+CHI_COEF(K,2))*T1+CHI_COEF(K,3))*T1+CHI_COEF(K,4)
	        CHI_RAY(I)=EXP(T2)
	        T2=((ETA_COEF(K,1)*T1+ETA_COEF(K,2))*T1+ETA_COEF(K,3))*T1+ETA_COEF(K,4)
	        ETA_RAY(I)=EXP(T2)
	        T2=(3.0D0*CHI_COEF(K,1)*T1+2.0D0*CHI_COEF(K,2))*T1+CHI_COEF(K,3)
	        dCHIdR_RAY(I)=T2*CHI_RAY(I)/R_RAY(I,LS)
	       END IF
	     END DO
!
	    SOURCE_RAY(1:NI)=ETA_RAY(1:NI)/CHI_RAY(1:NI)
!
! Even in the THK case we assume that the intensity incident on the
! outer boundary is zero. In THK case we effectively get IBOUND not
! equal zero at the model atmosphere boundary as a consequence of the
! extension.
!
	    IBOUND(LS)=0.0D0
!
! 
! By setting PF(1)=0 when evaluating SOURCE we ensure a pure continuum 
! calculation for the first frequency.
!
	    IF(INIT)THEN
	      DO I=1,NI
	        Q(I)=0.0D0
	        QH(I)=0.0D0
	      END DO                
	      OLDCHI(LS)=CHI_RAY(NI)
	    ELSE IF(NEW_FREQ)THEN
	      DO I=1,NI-1
	        QH(I)=GAMH(I,LS)*2.0D0/((CHI_RAY(I)+CHI_RAY(I+1))*dLOG_NU)
	        Q(I)=GAM(I,LS)/(CHI_RAY(I)*dLOG_NU)
	      END DO
	      QH(NI)=0.0D0
	      Q(NI)=GAM(NI,LS)/(CHI_RAY(NI)*dLOG_NU)
	    END IF
!
	    IF(DIF .AND. LS .LE. NC)THEN
	      T1=0.0D0
	      IF(.NOT. INIT)T1=GAM(NI,LS)/(CHI_RAY(NI)*dLOG_NU)	     !Q(NI)
	      DBC=DBB*SQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)/CHI_RAY(NI)
	1            *(1.0D0+T1*(1.0D0-CHI_RAY(NI)/OLDCHI(LS)))
	    END IF
	    OLDCHI_STORE(LS)=CHI_RAY(NI)
!
	    IF(NEW_FREQ)THEN
!
! Compute the optical depth increments. This code is from TAU, and NORDTAU. We
! check that the Euler-Mauclarin correction is not too large. This is mainly
! done to prevent negative optical depths. The check is not necessary when
! we re using monotonic interpolation. 
!
	      IF(METHOD .EQ. 'ZERO')THEN
	        DO I=1,NI-1
	          dZ=Z(I,LS)-Z(I+1,LS)
	          DTAU(I,LS)=0.5D0*(CHI_RAY(I)+CHI_RAY(I+1))*dZ
	        END DO  
	      ELSE IF(INSERT .OR. METHOD(4:6) .EQ. 'MON')THEN
	        DO I=1,NI-1
	          dZ=Z(I,LS)-Z(I+1,LS)
	          DTAU(I,LS)=0.5D0*dZ*( CHI_RAY(I)+CHI_RAY(I+1) +
	1            dZ*( dCHIdR_RAY(I+1)*Z(I+1,LS)/R_RAY(I+1,LS) -
	1            dCHIdR_RAY(I)*Z(I,LS)/R_RAY(I,LS) )/6.0D0 )
	        END DO
	      ELSE
	        DO I=1,NI-1
	          dZ=Z(I,LS)-Z(I+1,LS)
	          DTAU(I,LS)=0.5D0*dZ*( CHI_RAY(I)+CHI_RAY(I+1) +
	1            dZ*( dCHIdR_RAY(I+1)*Z(I+1,LS)/R_RAY(I+1,LS) -
	1            dCHIdR_RAY(I)*Z(I,LS)/R_RAY(I,LS) )/6.0D0 )
	          IF( CHI_RAY(I) .LT. CHI_RAY(I+1) )THEN
	            DTAU(I,LS)=MAX(CHI_RAY(I)*dZ,DTAU(I,LS))
	            DTAU(I,LS)=MIN(CHI_RAY(I+1)*dZ,DTAU(I,LS))
	          ELSE
	            DTAU(I,LS)=MIN(CHI_RAY(I)*dZ,DTAU(I,LS))
	            DTAU(I,LS)=MAX(CHI_RAY(I+1)*dZ,DTAU(I,LS))
	          END IF
	       END DO
	      END IF
!
	      CALL TUVGHD_RH(TA(1,LS),TB(1,LS),TC(1,LS),
	1            U,VB,VC,GB(1,LS),H(1,LS),XM,
	1            Q,QH,DTAU(1,LS),SOURCE_RAY,DIF,DBC,IC,LS,NC,NI)
	      XM(1)=-IBOUND(LS)
!
! Update PAR_AV matrix. U, VC, VC, AV_PREV, and CV_PREV do not change on
! successive calls with the same frequency. Thus PAR_AV can be used on
! successive calls with the same frequency.
!
	      PAR_AV(1,LS)=U(1)*AV_PREV(1,LS)
	      DO I=2,NI-1
	         PAR_AV(I,LS)=U(I)*AV_PREV(I,LS)-
	1             (VB(I)*CV_PREV(I-1,LS)+VC(I)*CV_PREV(I,LS)) 
	      END DO
	      PAR_AV(NI,LS)=U(NI)*AV_PREV(NI,LS)-VB(NI)*CV_PREV(NI-1,LS)
!
	      TC(NI,LS)=0				!As used.
	      DIV(1)=1.0/(TC(1,LS)+TB(1,LS))
	      TC(1,LS)=TC(1,LS)*DIV(1)
	      TB(1,LS)=TB(1,LS)*DIV(1)
	      DO I=2,NI
	        DIV(I)=1.0D0/(TA(I,LS)*TB(I-1,LS)+TB(I,LS)+TC(I,LS))
	        TB(I,LS)=(TA(I,LS)*TB(I-1,LS)+TB(I,LS))*DIV(I)
	        TC(I,LS)=TC(I,LS)*DIV(I)
	      END DO
	      TB(1:NI,LS)=DIV(1:NI)
	    ELSE
!
	      XM(1)=0.0D0
	      DO I=2,NI-1
	        XM(I)=-0.5D0*SOURCE_RAY(I)*(DTAU(I-1,LS)+DTAU(I,LS))
	      END DO
!
	      IF(LS .GT. NC)THEN
	        XM(NI)=0.5*DTAU(NI-1,LS)*SOURCE_RAY(NI)
	      ELSE IF(DIF)THEN
	        XM(NI)=DBC
	      ELSE
	        XM(NI)=IC
	     END IF
!
	    END IF
!
! Update AV matrix.
!
	    DO I=1,NI
	      AV(I,LS)=XM(I)+PAR_AV(I,LS)
	    END DO
!
	  END DO				!LS
!
! 
! 
! Solve for the radiation field along ray for this frequency.
!
	  IF(VECTOR_MACHINE)THEN
!
! This section is for a VECTOR machine. Loop over inner index is outer
! loop to remove a dependency. This is inefficient on scaler machines
! as array is not accessed sequentially.
!
!
! Forward substitution.
!
	    AV(1,1:NP)=AV(1,1:NP)*TB(1,1:NP)
	    DO I=2,NRAY_MAX
	      DO LS=1,MAX_LS(I)
	        IF(I .LE. NI_RAY(LS))
	1           AV(I,LS)=(AV(I,LS)+TA(I,LS)*AV(I-1,LS))*TB(I,LS)
	      END DO
	    END DO
C
C Backward substitution.
C
	    DO LS=1,NP
	      AV(NI_RAY(LS),LS)=-AV(NI_RAY(LS),LS)
	    END DO
	    DO I=NRAY_MAX,1,-1
	      DO LS=1,MAX_LS(I)
	        IF(I .LT. NI_RAY(LS))
	1          AV(I,LS)=TC(I,LS)*AV(I+1,LS)-AV(I,LS)
	      END DO
	    END DO
	  ELSE
C
C This section of code is for a scaler machine where the dependence
C on a previous computation does not matter. More efficient than previous
C code as array is accessed in correct manner.
C
C Forward substitution.
C
	    DO LS=1,NP
	      AV(1,LS)=AV(1,LS)*TB(1,LS)
	      DO I=2,NI_RAY(LS)
	        AV(I,LS)=(AV(I,LS)+TA(I,LS)*AV(I-1,LS))*TB(I,LS)
	      END DO
	    END DO
C
C Backward substitution.
C
	    DO LS=1,NP
	      AV(NI_RAY(LS),LS)=-AV(NI_RAY(LS),LS)
	      DO I=NI_RAY(LS)-1,1,-1
	        AV(I,LS)=TC(I,LS)*AV(I+1,LS)-AV(I,LS)
	      END DO
	    END DO
	  END IF
!	  CALL TUNE(2,'LS_LOOP')
!
!
!
! Verify validity of AV values. If these go negative, we set them +ve 
! but a factor of 10 smaller. We also ensure that the AV are not 
! extremely close to zero by comparing with neighboring AV values.
! The CV (fluxes) are computed using the revised values.
!
	  DO LS=1,NP
	   NI=NI_RAY(LS)
!
	    K=2
	    DO I=1,NI
	      T1=1.0D-08*ABS(AV(K,LS))
	      IF(AV(I,LS) .LE. T1)THEN
	        AV(I,LS)=MAX(0.01D0*ABS(AV(I,LS)),T1)
	        NEG_AV_VALUE=.TRUE.
	      ELSE
	        K=K+1
	      END IF
	      IF(I .EQ. 1)K=1
	    END DO
!
! Update C vector (i.e. flux variable).
!                          
	    DO I=1,NI_RAY(LS)-1
	      CV(I,LS)=GB(I,LS)*(AV(I,LS)-AV(I+1,LS))+H(I,LS)*CV_PREV(I,LS)
	    END DO
!
! Note that V=AV(1)-IBOUND.
!
	    IF(ND_ADD .EQ. 0)THEN
	      CV_BOUND(LS)=AV(1,LS)-IBOUND(LS)
	    ELSE IF(LS .EQ. NP)THEN
	      CV_BOUND(LS)=0.0D0
	    ELSE
	      K=J_PNT(1,LS)
	      CV_BOUND(LS)=0.5D0*(CV(K+1,LS)+CV(K-1,LS))
	    END IF
	    I_M_IN_BND(LS)=AV(NI,LS)-(AV(NI,LS)-AV(NI-1,LS))/DTAU(NI-1,LS)
	  END DO
!
	  AV_STORE(:,:)=AV(:,:)
	  CV_STORE(:,:)=CV(:,:)
!
	END IF			!Solution method
!
!************************************************************************
!************************************************************************
!
	CALL TUNE(1,'FG_SOL_INT')
	IF(SOLUTION_METHOD .EQ. 'INTEGRAL')THEN
!
! Enter loop to perform integration along each ray.
!
	  DO LS=1,NP
!
	    IF(NI_RAY(LS) .EQ. 1)THEN
	      I_P(1,LS)=0.0D0
	      I_M(1,LS)=0.0D0
	      GOTO 1000
	    END IF
!
! NB: dCHIdR = dLOG(CHI)/dLOG(R) * CHI/R
!
!	    CALL TUNE(1,'FG_CHI_Q')
	    NI=NI_RAY(LS)
	    DO I=1,NI_RAY(LS)
	      K=RAY_PNT(I,LS)
	      IF(R_RAY(I,LS) .EQ. R_EXT(K))THEN
	        CHI_RAY(I)=CHI_EXT(K)
	        ETA_RAY(I)=ETA_EXT(K)
	        dCHIdR_RAY(I)=CHI_COEF(K,3)*CHI_RAY(I)/R_RAY(I,LS)
	      ELSE
	        T1=LOG(R_RAY(I,LS)/R_EXT(K))
	        T2=((CHI_COEF(K,1)*T1+CHI_COEF(K,2))*T1+CHI_COEF(K,3))*T1+CHI_COEF(K,4)
	        CHI_RAY(I)=EXP(T2)
	        T2=((ETA_COEF(K,1)*T1+ETA_COEF(K,2))*T1+ETA_COEF(K,3))*T1+ETA_COEF(K,4)
	        ETA_RAY(I)=EXP(T2)
	        T2=(3.0D0*CHI_COEF(K,1)*T1+2.0D0*CHI_COEF(K,2))*T1+CHI_COEF(K,3)
	        dCHIdR_RAY(I)=T2*CHI_RAY(I)/R_RAY(I,LS)
	      END IF
	    END DO
!
	    IF(DIF .AND. LS .LE. NC)THEN
	      I_CORE=( ETA_RAY(NI)+
	1         DBB*SQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND) )/CHI_RAY(NI)
	    ELSE 
	      I_CORE=IC
	    END IF
!
! By setting PF(1)=0 when evaluating SOURCE we ensure a pure continuum 
! calculation for the first frequency.
!
	    IF(INIT)THEN
	      Q(1:NI)=0.0D0
	      SOURCE_RAY(1:NI)=ETA_RAY(1:NI)/CHI_RAY(1:NI)
	    ELSE
	      Q(1:NI)=GAM(1:NI,LS)/dLOG_NU
	      CHI_RAY(1:NI)=CHI_RAY(1:NI)+Q(1:NI)
	      dCHIdR_RAY(1:NI)=dCHIdR_RAY(1:NI)+dGAMdR(1:NI,LS)/dLOG_NU
	      Q(1:NI)=Q(1:NI)/CHI_RAY(1:NI)
	      SOURCE_RAY(1:NI)=ETA_RAY(1:NI)/CHI_RAY(1:NI)
	    END IF
!	    CALL TUNE(2,'FG_CHI_Q')
!
!
!
!	    CALL TUNE(1,'FG_NEW_F')
	    IF(NEW_FREQ)THEN
!                       
! Compute the optical depth increments. This code is from TAU, and NORDTAU. We
! check that the Euler-Mauclarin correction is not too large. This is mainly
! done to prevent negative optical depths. The check is still necessary when
! we re using monotonic interpolation, since the monotonic interpolation only
! applies to CHI. CHI_RAY contains an additional term due to the frequency
! derivative. As non relativistic, DTAU is the same for both directions.
!
	      IF(METHOD .EQ. 'ZERO')THEN
	        DO I=1,NI-1
	          dZ=Z(I,LS)-Z(I+1,LS)
	          DTAU(I,LS)=0.5D0*(CHI_RAY(I)+CHI_RAY(I+1))*dZ
	        END DO
	      ELSE
	        DO I=1,NI-1
	          dZ=Z(I,LS)-Z(I+1,LS)
	          DTAU(I,LS)=0.5D0*dZ*( CHI_RAY(I)+CHI_RAY(I+1) +
	1            dZ*( dCHIdR_RAY(I+1)*Z(I+1,LS)/R_RAY(I+1,LS) -
	1            dCHIdR_RAY(I)*Z(I,LS)/R_RAY(I,LS) )/6.0D0 )
	          IF( CHI_RAY(I) .LT. CHI_RAY(I+1) )THEN
	            DTAU(I,LS)=MAX(CHI_RAY(I)*dZ,DTAU(I,LS))
	            DTAU(I,LS)=MIN(CHI_RAY(I+1)*dZ,DTAU(I,LS))
	          ELSE
	            DTAU(I,LS)=MIN(CHI_RAY(I)*dZ,DTAU(I,LS))
	            DTAU(I,LS)=MAX(CHI_RAY(I+1)*dZ,DTAU(I,LS))
	          END IF
	        END DO
	      END IF
!
! Compute the functions used to evaluate the weighted integral over the
! polynomial fit to the source function.
!
! NB:  En(I)= EXP(-DTAU) {Integral[0 to DTAU] t^n EXP(t) dt }/ DTAU^n
!
	      DO I=1,NI-1
	        T1=DTAU(I,LS)
	        EE(I)=0.0D0
	        IF(T1 .LT. 700.0D0)EE(I)=EXP(-T1)
	        IF(T1 .GT. 0.5)THEN
	          E0(I)=1.0D0-EE(I)
	          E1(I)=1.0D0-E0(I)/T1
	          E2(I)=1.0D0-2.0D0*E1(I)/T1
	          E3(I)=1.0D0-3.0D0*E2(I)/T1
	        ELSE IF(T1 .GT. 0.1)THEN
	          E3(I)=0.25D0*T1*( 1.0D0-0.20*T1*
	1               (1.0D0-T1/6.0D0*(1.0D0-T1/7.0D0*
	1               (1.0D0-T1/8.0D0*(1.0D0-T1/9.0D0*
	1               (1.0D0-T1/10.0D0*(1.0D0-T1/11.0D0*
	1               (1.0D0-T1/12.0D0*(1.0D0-T1/13.0D0)))))))) )
	          E2(I)=T1*( 1.0D0-E3(I) )/3.0D0
	          E1(I)=T1*( 1.0D0-E2(I) )/2.0D0
	          E0(I)=T1*( 1.0D0-E1(I) )
	        ELSE
	          E3(I)=0.25D0*T1*( 1.0D0-0.20*T1*
	1               (1.0D0-T1/6.0D0*(1.0D0-T1/7.0D0*
	1               (1.0D0-T1/8.0D0*(1.0D0-T1/9.0D0) ))))
	          E2(I)=T1*( 1.0D0-E3(I) )/3.0D0
	          E1(I)=T1*( 1.0D0-E2(I) )/2.0d0
	          E0(I)=T1*( 1.0D0-E1(I) )
	        END IF
	      END DO
!
              DO I=1,NI-1
	        A0(I,LS)=EE(I)
	        A1(I,LS)=E0(I)-3.0D0*E2(I)+2.0D0*E3(I)
	        A2(I,LS)=3.0D0*E2(I)-2.0D0*E3(I)
	        A3(I,LS)=DTAU(I,LS)*(E1(I)-2.0D0*E2(I)+E3(I))
	        A4(I,LS)=DTAU(I,LS)*(E3(I)-E2(I))
	      END DO
	    END IF
!	    CALL TUNE(2,'FG_NEW_F')
!
! ******************* INWARD DIRECTED RAYS *********************************
!
! Compute the Source function for inward directed rays, and find the
! monotonic interpolating polynomial.
!
!	    CALL TUNE(1,'FG_INT_INT')
!
	    SOURCE_PRIME(1:NI)=SOURCE_RAY(1:NI)+Q(1:NI)*I_M_PREV(1:NI,LS)
	    DO I=1,NI-1
	      S(I)=(SOURCE_PRIME(I+1)-SOURCE_PRIME(I))/DTAU(I,LS)
	    END DO
!
! Now compute the derivatives node I. 
!
	    dS(1)=S(1) +(S(1)-S(2))*DTAU(1,LS)/(DTAU(1,LS)+DTAU(2,LS))
	    DO I=2,NI-1
	      dS(I)=(S(I-1)*DTAU(I,LS)+S(I)*DTAU(I-1,LS))/
	1                       (DTAU(I-1,LS)+DTAU(I,LS))
	    END DO
	    dS(NI)=S(NI-1)+(S(NI-1)-S(NI-2))*DTAU(NI-1,LS)/
	1                       (DTAU(NI-2,LS)+DTAU(NI-1,LS))
!             
! Adjust first derivatives so that function is monotonic  in each interval.
!
	    dS(1)=( SIGN(ONE,S(1))+SIGN(ONE,dS(1)) )*
	1                      MIN(ABS(S(1)),0.5*ABS(dS(1)))
	    DO I=2,NI-1
	      dS(I)=( SIGN(ONE,S(I-1))+SIGN(ONE,S(I)) )*
	1               MIN(ABS(S(I-1)),ABS(S(I)),0.5*ABS(dS(I)))
	    END DO
	    dS(NI)=( SIGN(ONE,S(NI-1))+SIGN(ONE,dS(NI)) )*
	1               MIN(ABS(S(NI-1)),0.5*ABS(dS(NI)))
!
            I_M(1,LS)=0.0D0
	    DO I=1,NI-1
	      I_M(I+1,LS)=I_M(I,LS)*A0(I,LS)+ (
	1             SOURCE_PRIME(I)*A1(I,LS)
	1        +    SOURCE_PRIME(I+1)*A2(I,LS)
	1        +    dS(I)*A3(I,LS)
	1        +    dS(I+1)*A4(I,LS) )
	    END DO
!
! ******************* OUTWARD DIRECTED RAYS *********************************
!
! Compute the Source function for outward directed rays, and find the
! monotonic interpolating polynomial.
!
	    SOURCE_PRIME(1:NI)=SOURCE_RAY(1:NI)+Q(1:NI)*I_P_PREV(1:NI,LS)
	    DO I=1,NI-1
	      S(I)=(SOURCE_PRIME(I+1)-SOURCE_PRIME(I))/DTAU(I,LS)
	    END DO
!
! Now compute the derivatives at node I. 
!
	    dS(1)=S(1) +(S(1)-S(2))*DTAU(1,LS)/(DTAU(1,LS)+DTAU(2,LS))
	    DO I=2,NI-1
	      dS(I)=(S(I-1)*DTAU(I,LS)+S(I)*DTAU(I-1,LS))/
	1                  (DTAU(I-1,LS)+DTAU(I,LS))
	    END DO
	    dS(NI)=S(NI-1)+(S(NI-1)-S(NI-2))*DTAU(NI-1,LS)/
	1                  (DTAU(NI-2,LS)+DTAU(NI-1,LS))
!             
! Adjust the first derivatives so that function is monotonic in each interval.
!
	    dS(1)=( SIGN(ONE,S(1))+SIGN(ONE,dS(1)) )*
	1                      MIN(ABS(S(1)),0.5*ABS(dS(1)))
	    DO I=2,NI-1
	      dS(I)=( SIGN(ONE,S(I-1))+SIGN(ONE,S(I)) )*
	1          MIN(ABS(S(I-1)),ABS(S(I)),0.5*ABS(dS(I)))
	    END DO
	    dS(NI)=( SIGN(ONE,S(NI-1))+SIGN(ONE,dS(NI)) )*
	1            MIN(ABS(S(NI-1)),0.5*ABS(dS(NI)))
!
	    IF(LS .LE. NC)THEN
	      I_P(NI,LS)=I_CORE
	    ELSE
	      I_P(NI,LS)=I_M(NI,LS)
	    END IF
	    DO I=NI-1,1,-1
	      I_P(I,LS)=I_P(I+1,LS)*A0(I,LS)+ (
	1             SOURCE_PRIME(I+1)*A1(I,LS)
	1        +    SOURCE_PRIME(I)*A2(I,LS)
	1        -    dS(I+1)*A3(I,LS)
	1        -    dS(I)*A4(I,LS) )
	    END DO
C
C Note that V=AV(1)-IBOUND.
C
1000	    CONTINUE
	    K=J_PNT(1,LS)
	    CV_BOUND(LS)=0.5D0*(I_P(K,LS)-I_M(K,LS))
	    I_M_IN_BND(LS)=I_M(NI,LS)
!
!	    CALL TUNE(2,'FG_INT_INT')
	  END DO			!LS
!
! Compute the mean intensity like variable U at each grid point. Used to
! compute J.
!
	  DO LS=1,NP
	    DO I=1,NI_RAY(LS)
	      AV(I,LS)=0.5D0*( I_P(I,LS)+I_M(I,LS) )
	    END DO
	  END DO
!
! We define CV on the gregular grid. We perform the interpolation onto
! the original R grid as we compute H.
!
	  DO LS=1,NP
	    DO I=1,NI_RAY(LS)
	      CV(I,LS)= 0.5D0*( I_P(I,LS)-I_M(I,LS) )
	    END DO
	  END DO
!
	  I_P_STORE(:,:)=I_P(:,:)
	  I_M_STORE(:,:)=I_M(:,:)
!
!	  IF(NEW_FREQ)THEN
!	     WRITE(117,*)FREQ
!	     DO LS=1,NP
!	       WRITE(117,'(5ES15.6)')P(LS),I_M_PREV(1,LS),I_P_PREV(1,LS), I_M_PREV(2,LS),I_P_PREV(2,LS)
!	     END DO
!	  END IF
!
	END IF				!Solution method
	CALL TUNE(2,'FG_SOL_INT')
!
!
!***************************************************************************
!***************************************************************************
!
!
! Zero boundary conditions.
!
	HBC=0.0D0			!H/J at model outer boundary.
	NBC=0.0D0			!N/J at model outer boundary.
	IN_HBC=0.0D0
!
! Initialize intensity matrices.
!
	JNU(:)=0.0D0			!1:ND
	HNU(:)=0.0D0
	KNU(:)=0.0D0
	NNU(:)=0.0D0
	RSQN_ON_RSQJ(:)=0.0D0
!
	CALL TUNE(1,'J_LOOP')
	DO LS=1,NP
	  NI_SMALL=ND-(LS-NC-1);   IF(LS .LE. NC+1)NI_SMALL=ND
	  DO I=1,NI_SMALL
	    K=J_PNT(I,LS)
	    JNU(I)=JNU(I)+JQW(I,LS)*AV(K,LS)
	    KNU(I)=KNU(I)+KQW(I,LS)*AV(K,LS)
	  END DO
	END DO
!
! Perform the interpolation of CV onto the mid points of the R grid, and
! compute H & N on this grid. With the INTEGRAL method, we will always
! need to interpolate, since, V is always defined on the grid points.
!
! For the DIFFERNECE method, we will only need to interpolate when we have
! included additional points along each ray.
!
	IF(SOLUTION_METHOD .EQ. 'INTEGRAL')THEN
	  DO LS=1,NP
	    NI_SMALL=ND-(LS-NC-1);   IF(LS .LE. NC+1)NI_SMALL=ND
	    DO I=1,NI_SMALL-1
	      K=H_PNT(I,LS)
	      T2=0.5D0*(R(I)+R(I+1))
	      T1=(T2-R_RAY(K,LS))/(R_RAY(K+1,LS)-R_RAY(K,LS))
	      HNU(I)=HNU(I)+HQW(I,LS)*((1.0D0-T1)*CV(K,LS)+T1*CV(K+1,LS))
	      NNU(I)=NNU(I)+NQW(I,LS)*((1.0D0-T1)*CV(K,LS)+T1*CV(K+1,LS))
	    END DO
	  END DO			!End do LS
	ELSE
	  DO LS=1,NP
	    NI_SMALL=ND-(LS-NC-1);   IF(LS .LE. NC+1)NI_SMALL=ND
	    DO I=1,NI_SMALL-1
	      K=H_PNT(I,LS)
	      T2=0.5D0*(R(I)+R(I+1))
	      T1=( 2.0*T2-(R_RAY(K,LS)+R_RAY(K+1,LS)) )/(R_RAY(K+2,LS)-R_RAY(K,LS))
	      HNU(I)=HNU(I)+HQW(I,LS)*((1.0D0-T1)*CV(K,LS)+T1*CV(K+1,LS))
	      NNU(I)=NNU(I)+NQW(I,LS)*((1.0D0-T1)*CV(K,LS)+T1*CV(K+1,LS))
	    END DO
	  END DO			!End do LS
	END IF
!
	DO LS=1,NP
	  NI=NI_RAY(LS)
	  HBC=HBC+CV_BOUND(LS)*H_WGHTS(LS)
	  NBC=NBC+CV_BOUND(LS)*N_WGHTS(LS)
	  IPLUS_P(LS)=2.0D0*CV_BOUND(LS)
	  IF(LS .LE. NC+1)THEN
	    IN_HBC=IN_HBC + JQW(ND,LS)*I_M_IN_BND(LS)*Z(NI,LS)/R_RAY(NI,LS)
	  END IF
	END DO
!
	CALL TUNE(2,'J_LOOP')
!
!
!
! Compute the boundary Eddington factors.
!
	HN_VAL=HBC+SIGMA(1)*NBC
	HBC=HBC/JNU(1)
	NBC=NBC/JNU(1)
	IN_HBC=IN_HBC/(2.0D0*JNU(ND)-IC)
!
!	NBC=(HN_VAL_PREV-HN_VAL)/JNU(1)a
!	IN_HBC=HN_VAL_PREV
!
!	HBC=JNU(2)/JNU(1)
!
! Compute the Eddington F and G factors, which are defined by
! F=K/J and G=N/H. The F and G factors are returned in KNU and NNU
! respectively.
!
	DO I=1,ND
	  KNU(I)=KNU(I)/JNU(I)
	END DO
!
! Now compute the "G" Eddington factors. Because H may be zero we have to
! be careful. Depending on N_TYPE we can specify N in terms of J, N in
! terms of H, or N in terms of J and H.
!
	IF(N_TYPE .EQ. 'N_ON_J')THEN
	  DO I=1,ND-1
	    ALPHA=0.25D0*(R(I)+R(I+1))*(R(I)+R(I+1))
	    RSQN_ON_RSQJ(I)=ALPHA*NNU(I)/
	1                  (R(I)*R(I)*JNU(I)+R(I+1)*R(I+1)*JNU(I+1))
	    NNU(I)=0.0D0
	  END DO
	ELSE IF(N_TYPE .EQ. 'MIXED')THEN
	  RSQN_ON_RSQJ(1:ND-1)=0.0D0
	  NNU(1)=NNU(1)/HNU(1)
	  DO I=2,ND-1
	    IF(HNU(I) .NE. 0)THEN
	      T1=NNU(I)/HNU(I)
	    ELSE
	      T1=100.0D0
	    END IF
	    IF(T1 .GT. 1.1 .OR. T1 .LT. 0.05)THEN
	      ALPHA=0.25D0*(R(I)+R(I+1))*(R(I)+R(I+1))
	      RSQN_ON_RSQJ(I)=ALPHA*NNU(I)/
	1                  (R(I)*R(I)*JNU(I)+R(I+1)*R(I+1)*JNU(I+1))
	      NNU(I)=0.0D0
	    ELSE
	      NNU(I)=T1
	    END IF
	  END DO
	ELSE IF(N_TYPE .EQ. 'G_ONLY')THEN
!
! Compute G Eddington factor storing in N. We also check the validity of 
! the Eddington factor in case strange values are occurring because H 
! is near zero (switching sign?).
!
	  DO I=1,ND-1
	    IF(HNU(I) .NE. 0)THEN
	      NNU(I)=NNU(I)/HNU(I)
	      IF(NNU(I) .GT. 1.0D0)THEN
                BAD_NNU_VALUE=.TRUE.
	        NNU(I)=1.0D0
	      ELSE IF( NNU(I) .LT. 0.01) THEN
	        NEG_NNU_VALUE=.TRUE.
	        NNU(I)=0.01D0
	      END IF
	    ELSE
	      NNU(I)=0.0D0		!Later replaced by average.
	      J=ERROR_LU()
	      WRITE(J,'(1X,A,1PE16.8)')'HNU zero for frequency:',FREQ
	      WRITE(J,'(1X,A,I4)')'Error occurred at depth:',I
	    END IF
	  END DO
!
	ELSE
	  J=ERROR_LU()
	  WRITE(J,*)'Unrecognized N_TYPE in FG_J_CMF_V4'
	  WRITE(J,*)'N_TYPE=',N_TYPE
	  STOP
	END IF
!
! Store frequencies at which errors occurred, and give an indication of the
! error.
!
	J=1
	DO WHILE (J .LE. N_ERR_MAX .AND. 
	1            (NEG_AV_VALUE .OR. BAD_NNU_VALUE .OR. NEG_NNU_VALUE) )
	  IF(FG_ERR_ON_FREQ(J) .EQ. FREQ)THEN
	    FG_ERR_TYPE(J)=0
	    IF(NEG_AV_VALUE)FG_ERR_TYPE(J)=FG_ERR_TYPE(J)+1
	    IF(NEG_NNU_VALUE)FG_ERR_TYPE(J)=FG_ERR_TYPE(J)+5
	    IF(BAD_NNU_VALUE)FG_ERR_TYPE(J)=FG_ERR_TYPE(J)+10
	    NEG_AV_VALUE=.FALSE.
	    BAD_NNU_VALUE=.FALSE.
	    NEG_NNU_VALUE=.FALSE.
	  ELSE IF(J .EQ. FG_ERR_CNT+1)THEN
	    FG_ERR_CNT=J
	    FG_ERR_TYPE(J)=0
	    FG_ERR_ON_FREQ(J)=FREQ
	    IF(NEG_AV_VALUE)FG_ERR_TYPE(J)=FG_ERR_TYPE(J)+1
	    IF(NEG_NNU_VALUE)FG_ERR_TYPE(J)=FG_ERR_TYPE(J)+5
	    IF(BAD_NNU_VALUE)FG_ERR_TYPE(J)=FG_ERR_TYPE(J)+10
	    NEG_AV_VALUE=.FALSE.
	    BAD_NNU_VALUE=.FALSE.
	    NEG_NNU_VALUE=.FALSE.
	  END IF
	  J=J+1
	END DO
!
	DO I=1,ND
	  IF(JNU(I) .LT. 0)THEN
	    WRITE(111,*)'FREQ=',FREQ
	    WRITE(111,'(3X,A,5X,4(5X,A,5X))')'I','JNU','ETA','CHI','ESEC'
	    DO J=1,ND
	      WRITE(111,'(X,I5,4ES16.6)')J,JNU(J),ETA(J),CHI(J),ESEC(J)
	    END DO
	    J=ERROR_LU()
	    WRITE(J,*)'Error on FG_J_CMF_V10 --- negative mean intensities.'
	    WRITE(J,*)'Check out file FOR111 for aditional information.'
	    WRITE(J,*)'Halting code execution.'
	    STOP
	  END IF
	END DO
!
        ACCESS_F=5
	LU_IP=56
	IF(INIT)THEN
	  INQUIRE(FILE='IP_FG_DATA',EXIST=WRITE_IP)
          IF(WRITE_IP)THEN
            CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
            I=WORD_SIZE*(NP+1)/UNIT_SIZE
            CALL WRITE_DIRECT_INFO_V3(NP,I,'20-Aug-2000','IP_FG_DATA',LU_IP)
            OPEN(UNIT=LU_IP,FILE='IP_FG_DATA',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='UNKNOWN',RECL=I,IOSTAT=IOS)
	    FREQ_CNT=0
            WRITE(LU_IP,REC=3)ACCESS_F,FREQ_CNT,NP
            WRITE(LU_IP,REC=ACCESS_F)(P(I),I=1,NP)
	    WRITE_IP=.TRUE.
	  END IF
        END IF
	IF(WRITE_IP)THEN
	  T1=0.0D0
          IF(FREQ_CNT .NE. 0)READ(LU_IP,REC=ACCESS_F+FREQ_CNT)(IBOUND(LS),LS=1,NP),T1
	  IF(T1 .NE. FREQ)FREQ_CNT=FREQ_CNT+1
          WRITE(LU_IP,REC=3)ACCESS_F,FREQ_CNT,NP
          WRITE(LU_IP,REC=ACCESS_F+FREQ_CNT)(CV_BOUND(LS),LS=1,NP),FREQ
        END IF
	CALL TUNE(2,'FG_BIG_SEC')
!
	RETURN
	END
