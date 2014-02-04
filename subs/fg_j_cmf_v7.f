C
C Data module for FG_J_CMF. Data placed in this module is automatically
C saved between subroutine calls..
C
	MODULE FG_J_CMF_MOD_V7
C
C The *_STORE routines are used to store the radiation field, as computed
C on the previous call to FG_J_CMF. 
C The *_PREV routines are used to store the radiation field, as computed
C for the previous frequency.
C
C The *_PREV routines are updated from the *_STORE routines when NEW_FREQ
C is .TRUE.
C
	REAL*8, ALLOCATABLE :: OLDCHI(:)
	REAL*8, ALLOCATABLE :: OLDCHI_STORE(:)
	REAL*8, ALLOCATABLE :: AV_PREV(:,:)
	REAL*8, ALLOCATABLE :: AV_STORE(:,:)
	REAL*8, ALLOCATABLE :: CV_PREV(:,:)
	REAL*8, ALLOCATABLE :: CV_STORE(:,:)
	REAL*8, ALLOCATABLE :: R_RAY(:,:)
	REAL*8, ALLOCATABLE :: Z(:,:)
	REAL*8, ALLOCATABLE :: GAM(:,:)
	REAL*8, ALLOCATABLE :: GAMH(:,:)
C
	REAL*8, ALLOCATABLE :: R_INS(:,:)
	REAL*8, ALLOCATABLE :: V_INS(:,:)
	REAL*8, ALLOCATABLE :: SIGMA_INS(:,:)
	REAL*8, ALLOCATABLE :: CHI_INS(:,:)
	REAL*8, ALLOCATABLE :: dCHIdR_INS(:,:)
	REAL*8, ALLOCATABLE :: SOURCE_INS(:,:)

	REAL*8, ALLOCATABLE :: R_EXT(:)
	REAL*8, ALLOCATABLE :: V_EXT(:)
	REAL*8, ALLOCATABLE :: SIGMA_EXT(:)
C                            
C Used to compute the boundary iteration factors (HBC and NBC).
C                             
	REAL*8, ALLOCATABLE :: N_WGHTS(:)
	REAL*8, ALLOCATABLE :: H_WGHTS(:)
	REAL*8, ALLOCATABLE :: MU_VAL(:)
	INTEGER, ALLOCATABLE :: NI_RAY(:)
	INTEGER, ALLOCATABLE :: MAX_LS(:)
C
	LOGICAL FIRST_TIME
	DATA FIRST_TIME/.TRUE./
C
	INTEGER ND_EXT
	INTEGER ND_ADD
	INTEGER ND_MAX
	INTEGER NP_MAX
	INTEGER NRAY_MAX
C
	END MODULE FG_J_CMF_MOD_V7
C
C 
C
C Routine to compute the Eddington F, G and RSQN_ON_RSQJ Eddington factors 
C for single frequency transfer in the comoving-frame. Based on FG_J_CMF. 
C First order differencing in frequency is used.
C
C The F, G, and RSQN_ON_RSQJ Eddington factors must be supplied.
C
C NB:
C	F = K / J
C	G=  N / H
C	RSQN_ON_RSQJ(I) = RSQ_N(I)/( RSQ_J(I)+ RQS_J(I+1))
C
C where
C	RSQ_N=0.25*N(I)*(R(I)+R(I+1))**2
C
C NB: Only one of G and RSQN_ON_RSQJ is defined at a given depth. This
C     avoids having to test what mode I am using for the Eddington factors.
C
C     IF N_TYPE='G_ONLY' G is defined at all depths, and
c       RSQN_ON_RSQJ=0 at all depths.
C     IF N_TYPE='N_ON_J' RSQN_ON_RSQJ is defined at all 
C       depths, and G=0 at all depths.
C     IF N_TYPE='MIXED' one of G or RSQN_ON_RSQJ is 
C       non-zero, and is the value to be used in MOM_J_CMF
C
C Routine also returns I+, so that observers flux can be computed. Note because
C the outer boundary can be thick, I+ is actually I+ - I- and hence is a flux
C like variable (This I+ = 2v at outer boundary).
C
C The radiation field at the previous (bluer) frequency is stored locally.
C Coupling to this frequency is controlled by the variable INIT. If INIT is
C true, the atmosphere is treated with the assumption that V=0.
C
C INIT must be TRUE on the very first call.
C
	SUBROUTINE FG_J_CMF_V7(ETA,CHI,ESEC,V,SIGMA,R,P,
	1                  JNU,HNU,KNU,NNU,RSQN_ON_RSQJ,
	1                  JQW,HQW,KQW,NQW,
	1                  IN_HBC,HBC,NBC,
	1                  IPLUS_P,FL,dLOG_NU,DIF,DBB,IC,METHOD,
	1                  THK,INIT,NEW_FREQ,N_TYPE,NC,NP,ND)
	USE FG_J_CMF_MOD_V7
	IMPLICIT NONE
C
C Altered: 20-Oct-1998   Bug fix: dCHIdR was not being computed for extended
C                          region when CHI < ESEC at outer boundary.
C Altered: 13-Oct-1998   Bug fix: Incorrect computation of H in region of
C                          inserted points.
C Altered: 29-Sep-1997   Scaler loop code installed to improve speed on a
C                          an Alphastation.
C Altered: 21-Aug-1997   ESEC now longer assumed to vary as 1/r^2, but as a
C                           power law (steep or steeper than 1/r^2).
C Altered: 19_Dec-1996   Euler-MauClaurin correction to DTAU is limited to
C                           avoide negative DTAU's.
C Altered: 10-Sep-1996   3.5 used as limit on ETA rather than 2.0 to avoid
C                           excess outer envelope emission.
C Altered: 16-Jun-1996   dCHIdR zeroed if METHOD='ZERO' (needed for 
C                           interpolation section.
C Altered: 03-Feb-1996   Treatment of negative AV values altered.
C Altered: 25-Mar-1995   THK option reinstalled --- previously boundary was
C                          always assumed to be thick. TAU_BOUND deleted
C                          call. Version changed to _V5.
C                        HBC and NBC are now scalers.
C
C Altered: 21-Mar-1995   Factor of 2 missing from IPLUS definition. 
C                           Bug introduced when extended outer boundary.
C Altered:  Mid Feb-1995  !Improved outer boundary condition installed.
C     to    Mid Mar-1995     Boundary is now extended a factor of 10 in radius,
C                            and I- is set to zero at this boundary. Solved
C                            problem in outer b.c. caused by the very different
C                            line and continuum tau's and source functions, and
C                            their radial dependence.
C                         To handle low tau's, we use a new THOMAS routine
C                            (after Rybicki and Hummer). Several simple routines
C                            replaced by direct insertion.
C                         Optimization of code speed. Quantities such a Z
C                            are now only computed once (when INIT is set).
C                         Error reporting was inserted.via common block. 
C                            This allows reporting of those frequencies on which
C                            error occur (eg -ve intensities, zero flux).
C	                    (This comment inserted 21-Mar-95)
C
C Altered:   09-Mar-1995  RSQ_N_ONJ installed in effort to solve problem caused
C                           by zero H, and hence undefined G values.
C
C Altered   22-Dec-1994         !IPLU_P inserted
C Finalized 11-Nov-1994		!HTRPWGT included, and cleaned.
C Created   12-Oct-1994
C
	INTEGER NC,NP,ND
	REAL*8 ETA(ND),CHI(ND),ESEC(ND)
	REAL*8 V(ND),SIGMA(ND),R(ND),P(NP)
C
C NB: J,H,K,N refer to the first 4 moments of the radiation field.
C     QW denotes quadrature weight.
C
	REAL*8 JQW(ND,NP),HQW(ND,NP),KQW(ND,NP),NQW(ND,NP)
	REAL*8 JNU(ND),HNU(ND),KNU(ND),NNU(ND),RSQN_ON_RSQJ(ND)
	REAL*8 IN_HBC,HBC,NBC
	REAL*8 IPLUS_P(NP)
C
	REAL*8 DBB,IC,FL,dLOG_NU
	CHARACTER*6 METHOD
	CHARACTER*6 N_TYPE
	LOGICAL DIF		!Use diffusion approximation
C
C Use "Thick" boundary cond. at outer boundary. Only noted when INIT 
C is true. All subsequent frequencies will use the same boundary condition
C independent of the passed value (Until INIT is set to TRUE again).
C
	LOGICAL THK
C
C First frequency -- no frequency coupling.

	LOGICAL INIT
C
C Upon leaving this routine the radiation field along each ray is stored. This
C will provide the blue wing information necessary for the next frequency.
C This routine may, however, be used in an iterative loop. In this case the
C "blue wing" information should remain unaltered between calls.
C NEW_FREQ indicates that a new_frequency is being passed, and hence the "blue
C wing" information should be updated.
C
	LOGICAL NEW_FREQ	
C
	INTEGER ND_ADD_MAX
	PARAMETER (ND_ADD_MAX=24)
C
C The following arrays do not need to be stored, and hence can be crteated 
C each time.
C
	REAL*8 TA(ND+ND_ADD_MAX+6,NP)
	REAL*8 TB(ND+ND_ADD_MAX+6,NP)
	REAL*8 TC(ND+ND_ADD_MAX+6,NP)
	REAL*8 AV(ND+ND_ADD_MAX+6,NP)
	REAL*8 CV(ND+ND_ADD_MAX+6,NP)
	REAL*8 GB(ND+ND_ADD_MAX+6,NP)
	REAL*8 H(ND+ND_ADD_MAX+6,NP)
C
	REAL*8 DTAU(ND+ND_ADD_MAX+6)
	REAL*8 XM(ND+ND_ADD_MAX+6)
	REAL*8 SOURCE(ND+ND_ADD_MAX+6)
	REAL*8 U(ND+ND_ADD_MAX+6)
	REAL*8 VB(ND+ND_ADD_MAX+6)
	REAL*8 VC(ND+ND_ADD_MAX+6)
	REAL*8 Q(ND+ND_ADD_MAX+6)
	REAL*8 QH(ND+ND_ADD_MAX+6)
	REAL*8 dCHIdR(ND+ND_ADD_MAX+6)
C
	REAL*8 ETA_EXT(ND+ND_ADD_MAX+6)
	REAL*8 CHI_EXT(ND+ND_ADD_MAX+6)
C
	REAL*8 V_RAY(ND+ND_ADD_MAX+6)
	REAL*8 SIGMA_RAY(ND+ND_ADD_MAX+6)
	REAL*8 SOURCE_RAY(ND+ND_ADD_MAX+6)
	REAL*8 CHI_RAY(ND+ND_ADD_MAX+6)
	REAL*8 dCHIdR_RAY(ND+ND_ADD_MAX+6)
C
	REAL*8 DIV(ND+ND_ADD_MAX+6)
	REAL*8 DTAU_BND(NP)
C
	INTEGER N_ERR_MAX,FG_ERR_CNT
	PARAMETER (N_ERR_MAX=1000)
	REAL*8 FG_ERR_ON_FREQ
	INTEGER FG_ERR_TYPE
	COMMON /FG_J_CMF_ERR/FG_ERR_ON_FREQ(N_ERR_MAX),
	1                    FG_ERR_TYPE(N_ERR_MAX),FG_ERR_CNT
	LOGICAL NEG_AV_VALUE,NEG_NNU_VALUE,BAD_NNU_VALUE
C
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
C
C Local variables.
C
	INTEGER, PARAMETER :: NINS=4
	LOGICAL, PARAMETER :: LFALSE=.FALSE.
	LOGICAL, PARAMETER :: LTRUE=.TRUE.
C                 
	INTEGER NI_SMALL
	INTEGER I,J,K,LS
	INTEGER NI
	REAL*8 DBC
	REAL*8 IBOUND			!Incident intensity on outer boundary.
	REAL*8 T1,T2
	REAL*8 ALPHA
	REAL*8 ESEC_POW
	REAL*8 BETA
	REAL*8 VINF
	REAL*8 RMAX,DEL_R_FAC
	REAL*8 CV_BOUND
	REAL*8 MU,dZ,PSQ
C
	REAL*8 DEL_R
C
C Change the following statement to TRUE if running on a VECTOR machine.
C
	LOGICAL, PARAMETER :: VECTOR_MACHINE=.FALSE.
C
	IF(FIRST_TIME)THEN
	  ND_ADD=0
	  IF(THK)ND_ADD=ND_ADD_MAX
	  ND_EXT=ND+ND_ADD
C
	  ND_MAX=ND+ND_ADD+6
	  NP_MAX=NP
	  ALLOCATE ( OLDCHI(NP) )
	  ALLOCATE ( OLDCHI_STORE(NP) )
	  ALLOCATE ( AV_PREV(ND_MAX,NP) )
	  ALLOCATE ( AV_STORE(ND_MAX,NP) )
	  ALLOCATE ( CV_PREV(ND_MAX,NP) )
	  ALLOCATE ( CV_STORE(ND_MAX,NP) )
	  ALLOCATE ( R_RAY(ND_MAX,NP) )
	  ALLOCATE ( Z(ND_MAX,NP) )
	  ALLOCATE ( GAM(ND_MAX,NP) )
	  ALLOCATE ( GAMH(ND_MAX,NP) )
C
	  ALLOCATE ( R_EXT(ND_EXT) )
	  ALLOCATE ( V_EXT(ND_EXT) )
	  ALLOCATE ( SIGMA_EXT(ND_EXT) )
C
C Must be dimensioned ND-1 for MON_INT_INS_V1 as their are only ND-1 intervals.
C
	  ALLOCATE ( R_INS(ND-1,NINS) )
	  ALLOCATE ( V_INS(ND-1,NINS) )
	  ALLOCATE ( SIGMA_INS(ND-1,NINS) )
	  ALLOCATE ( CHI_INS(ND-1,NINS) )
	  ALLOCATE ( dCHIdR_INS(ND-1,NINS) )
	  ALLOCATE ( SOURCE_INS(ND-1,NINS) )
C
	  ALLOCATE ( N_WGHTS(NP) )
	  ALLOCATE ( H_WGHTS(NP) )
	  ALLOCATE ( MU_VAL(NP) )
	  ALLOCATE ( NI_RAY(NP) )
	  ALLOCATE ( MAX_LS(ND_MAX) )
C
	  FIRST_TIME=.FALSE.
	END IF
C
	IF(NP .NE. NP_MAX)THEN
	  WRITE(ERROR_LU(),*)
	1    'Inconsistent NP and NP_MAX dimension in FG_J_CMF'
	  STOP
	END IF
	IF(ND+ND_ADD+6 .NE. ND_MAX)THEN
	  WRITE(ERROR_LU(),*)
	1    'Inconsistent ND and ND_MAX dimension in FG_J_CMF'
	  STOP
	END IF
C
	NEG_AV_VALUE=.FALSE.
	NEG_NNU_VALUE=.FALSE.
	BAD_NNU_VALUE=.FALSE.
C         
C Perform initializations. 
C
	IF(INIT)THEN
C
C Insert extra points into radius grid. Not all points will be used along a ray.
C We will only insert 2 additional points in the interval between Z=0 and 
C Z(1st) and 1 additional point in the next interval.
C
	  R_INS(:,:)=0.0D0
	  DO I=1,ND-1
	    DEL_R=R(I)-R(I+1)
	    R_INS(I,1)= R(I+1)+DEL_R/1.5D0
	    R_INS(I,2)= R(I+1)+DEL_R/3.0D0
	    R_INS(I,3)= R(I+1)+0.16D0*DEL_R		!4/25
	    R_INS(I,4)= R(I+1)+0.04D0*DEL_R		!1/25
	  END DO
C
C Perform monotonic interpolations in V and SIGMA. the dCHIR... arrays are
C not used as the derivatives are not computed (last passed variable).
C
	  CALL MON_INT_INS_V1(V_INS,R_INS,NINS,V,R,ND,
	1                   LFALSE,LFALSE,dCHIdR,dCHIdR_INS,LFALSE)
	  CALL MON_INT_INS_V1(SIGMA_INS,R_INS,NINS,SIGMA,R,ND,
	1                   LFALSE,LFALSE,dCHIdR,dCHIdR_INS,LFALSE)
C
	  FG_ERR_ON_FREQ(:)=0.0D0
	  FG_ERR_TYPE(:)=0
	  FG _ERR_CNT=0
	  AV_PREV(:,:)=0.0D0
	  CV_PREV(:,:)=0.0D0
C
	  DO LS=1,NP
	    MU_VAL(LS)=SQRT(R(1)*R(1)-P(LS)*P(LS))/R(1)
	    IPLUS_P(LS)=0.0D0
	  END DO
	  CALL HTRPWGT(MU_VAL,H_WGHTS,NP)
	  CALL NTRPWGT(MU_VAL,N_WGHTS,NP)
C
C Compute the extended R grid, excluding inserted points.
C
	  DO I=1,ND
	    R_EXT(ND_ADD+I)=R(I)
	  END DO
	  IF(THK)THEN
	    RMAX=10.0D0*R(1)
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
	  END IF
C
C Compute VEXT and R_EXT. We assume a BETA velocity law at large R.
C
	  V_EXT(ND_ADD+1:ND_EXT)=V(1:ND)
	  SIGMA_EXT(ND_ADD+1:ND_EXT)=SIGMA(1:ND)
	  IF(THK)THEN
	    BETA=(SIGMA(1)+1.0D0)*(R(1)/R(ND)-1.0D0)
            VINF=V(1)/(1-R(ND)/R(1))**BETA
	    DO I=1,ND_ADD
	      V_EXT(I)=VINF*(1-R_EXT(ND_EXT)/R_EXT(I))**BETA
	      SIGMA_EXT(I)=BETA/(R_EXT(I)/R_EXT(ND_EXT)-1.0D0)-1.0D0
	    END DO
	  END IF
C
	  DO LS=1,NP
	    PSQ=P(LS)*P(LS)
C
C Choose the points needed for each ray.
C
	     R_RAY(1:ND_EXT,LS)=R_EXT(1:ND_EXT)
C
C Only use V and SIGMA to comute GAMMA. Thus don't need to save there values
C
	     V_RAY(1:ND_EXT)=V_EXT(1:ND_EXT)
	     SIGMA_RAY(1:ND_EXT)=SIGMA_EXT(1:ND_EXT)
C
C We now insert some additional points around the end of ray as required.
C
	    IF(LS .LE. NC+1)THEN
	      NI_RAY(LS)=ND_EXT
	    ELSE
	      NI_RAY(LS)=ND_EXT-(LS-NC-1)
	      IF(NI_RAY(LS) .GT. ND_ADD+3)THEN
C
C In this case we insert 4 points between the last 2 points, and 2 points between 
C the second last pair.
C
	        J=NI_RAY(LS)
                R_RAY(J+6,LS)=R_RAY(J,LS)
                R_RAY(J+1,LS)=R_RAY(J-1,LS)
	        R_RAY(J-1,LS)=R_INS(J-2-ND_ADD,1)
	        R_RAY(J,LS)=R_INS(J-2-ND_ADD,2)
	        R_RAY(J+2,LS)=R_INS(J-1-ND_ADD,1)
	        R_RAY(J+3,LS)=R_INS(J-1-ND_ADD,2)
 	        R_RAY(J+4,LS)=R_INS(J-1-ND_ADD,3)
	        R_RAY(J+5,LS)=R_INS(J-1-ND_ADD,4)
C
C Do the same thing for V and SIGMA
C
                V_RAY(J+6)=V_RAY(J)
                V_RAY(J+1)=V_RAY(J-1)
	        V_RAY(J-1)=V_INS(J-2-ND_ADD,1) 
	        V_RAY(J)=V_INS(J-2-ND_ADD,2) 
	        V_RAY(J+2)=V_INS(J-1-ND_ADD,1) 
	        V_RAY(J+3)=V_INS(J-1-ND_ADD,2) 
	        V_RAY(J+4)=V_INS(J-1-ND_ADD,3) 
	        V_RAY(J+5)=V_INS(J-1-ND_ADD,4)
C
                SIGMA_RAY(J+6)=SIGMA_RAY(J)
                SIGMA_RAY(J+1)=SIGMA_RAY(J-1)
	        SIGMA_RAY(J-1)=SIGMA_INS(J-2-ND_ADD,1) 
	        SIGMA_RAY(J)=SIGMA_INS(J-2-ND_ADD,2) 
	        SIGMA_RAY(J+2)=SIGMA_INS(J-1-ND_ADD,1) 
	        SIGMA_RAY(J+3)=SIGMA_INS(J-1-ND_ADD,2) 
	        SIGMA_RAY(J+4)=SIGMA_INS(J-1-ND_ADD,3) 
	        SIGMA_RAY(J+5)=SIGMA_INS(J-1-ND_ADD,4)
C
	        NI_RAY(LS)=NI_RAY(LS)+6
	      ELSE IF(NI_RAY(LS) .EQ. ND_ADD+3)THEN
C
C In this case we insert 4 points between the last 2 points only. We don't do
C the insertion in the next zone, as the R(ND-1) and R(ND) are generally close
C together anyway.
C
	        J=NI_RAY(LS)
                R_RAY(J+4,LS)=R_RAY(J,LS)
	        R_RAY(J,LS)=R_INS(J-1-ND_ADD,1) 
	        R_RAY(J+1,LS)=R_INS(J-1-ND_ADD,2) 
	        R_RAY(J+2,LS)=R_INS(J-1-ND_ADD,3) 
	        R_RAY(J+3,LS)=R_INS(J-1-ND_ADD,4)
C
                V_RAY(J+4)=V_RAY(J)
	        V_RAY(J)=V_INS(J-1-ND_ADD,1) 
	        V_RAY(J+1)=V_INS(J-1-ND_ADD,2) 
	        V_RAY(J+2)=V_INS(J-1-ND_ADD,3) 
	        V_RAY(J+3)=V_INS(J-1-ND_ADD,4)
C
                SIGMA_RAY(J+4)=SIGMA_RAY(J)
	        SIGMA_RAY(J)=SIGMA_INS(J-1-ND_ADD,1) 
	        SIGMA_RAY(J+1)=SIGMA_INS(J-1-ND_ADD,2) 
	        SIGMA_RAY(J+2)=SIGMA_INS(J-1-ND_ADD,3) 
	        SIGMA_RAY(J+3)=SIGMA_INS(J-1-ND_ADD,4)
C
	        NI_RAY(LS)=NI_RAY(LS)+4
	      END IF
	    END IF
C
	    NI=NI_RAY(LS)
	    DO I=1,NI_RAY(LS)
	      Z(I,LS)=SQRT(R_RAY(I,LS)*R_RAY(I,LS)-P(LS)*P(LS))
	    END DO
C
C Compute GAMMA. This section is straight from the subroutine GAMMA, except
C That _EXT has been added to V, SIGMA, and R.
C
C We assume (1)	SIGMAd+1/2 = 0.5*( SIGMAd+1+SIGMAd )
C  	    (2)	Vd+1/2=0.5*( Vd + Vd+1 )
C Note that V is in km/s and SIGMA=(dlnV/dlnR-1.0)
C
	    DO I=1,NI_RAY(LS)-1
	      MU=Z(I,LS)/R_RAY(I,LS)
	      GAM(I,LS)=3.33564D-06*V_RAY(I)/R_RAY(I,LS)*
	1                   (  1.0D0+SIGMA_RAY(I)*(MU**2)  )
	      MU=(Z(I,LS)+Z(I+1,LS))/(R_RAY(I,LS)+R_RAY(I+1,LS))
	      GAMH(I,LS)=(V_RAY(I+1)+V_RAY(I))*3.33564D-06*
	1                   (  1.0D0+0.5D0*(MU**2)*
	1                   (SIGMA_RAY(I)+SIGMA_RAY(I+1))  )/
	1                   (R_RAY(I,LS)+R_RAY(I+1,LS))
	    END DO
	    NI=NI_RAY(LS)
	    GAMH(NI,LS)=0.0
	    MU=Z(NI,LS)/R_RAY(NI,LS)
	    GAM(NI,LS)=3.33564D-06*V_RAY(NI)*
	1              ( 1.0D0+SIGMA_RAY(NI)*(MU**2) )/R_RAY(NI,LS)
C
	  END DO		!LS Loop
C
C Determine maximum ray index having I points.
C
	  NRAY_MAX=MAXVAL(NI_RAY(1:NP))
	  DO I=1,NRAY_MAX
	    DO LS=1,NP
	      IF(NI_RAY(LS) .GE. I)MAX_LS(I)=LS
	    END DO
	  END DO
C
	ELSE IF(NEW_FREQ)THEN
	  OLDCHI(1:NP)=OLDCHI_STORE(1:NP)
	  AV_PREV(:,:)=AV_STORE(:,:)
	  CV_PREV(:,:)=CV_STORE(:,:)
	END IF
C
C Insure that an initialization has been performed.
C
	IF(ND_ADD .LT. 0)THEN
	  WRITE(ERROR_LU(),*)'Error in FG_J_CMF_V5 -- no',
	1                       '  initialization call'
	  STOP
	END IF
C
C 
C
C Initialize intensity matrices.
C
	JNU(:)=0.0D0			!1:ND
	HNU(:)=0.0D0
	KNU(:)=0.0D0
	NNU(:)=0.0D0
	RSQN_ON_RSQJ(:)=0.0D0
C
C Compute CHI_EXT, and ETA_EXT. CHI_EXT could be saved, as it doesn't change
C during the iteration procedure. ETA does however change (since it depends
C of J) and thus ETA_EXT must be re-computed on each entry.
C
C The first checks whether we may have negative line opacities due
C to stimulated emission. In such a case we simply assume an 1/r^2
C extrapolation.
C
C We also interpolate in ESEC, since ESEC (in the absence of negative 
C absorption) provides a lower bound to the opacity. NB: When CHI is much
C larger then ESEC its variation with r dominates, and it is possible to
C extrapolate CHI below ESEC.
C
	IF(ND_ADD .NE. 0)THEN
	  IF(CHI(1) .LE. ESEC(1) .OR. CHI(4) .LE. ESEC(4))THEN
	    ESEC_POW=LOG(ESEC(4)/ESEC(1))/LOG(R(1)/R(4))
	    IF(ESEC_POW .LT. 2.)ESEC_POW=2
	    DO I=1,ND_ADD
	      CHI_EXT(I)=CHI(1)*(R(1)/R_EXT(I))**ESEC_POW
	      dCHIdR(I)= -ESEC_POW*CHI_EXT(I)/R_EXT(I)
	    END DO
	  ELSE
	    ALPHA=LOG( (CHI(4)-ESEC(4)) / (CHI(1)-ESEC(1)) )
	1          /LOG(R(1)/R(4))
	    IF(ALPHA .LT. 2.)ALPHA=2.0
	    ESEC_POW=LOG(ESEC(4)/ESEC(1))/LOG(R(1)/R(4))
	    IF(ESEC_POW .LT. 2.)ESEC_POW=2
	    DO I=1,ND_ADD
	      T1=(CHI(1)-ESEC(1))*(R(1)/R_EXT(I))**ALPHA
	      T2=ESEC(1)*(R(1)/R_EXT(I))**ESEC_POW
	      CHI_EXT(I)=T1+T2
	      dCHIdR(I)=(-ALPHA*T1-ESEC_POW*T2)/R_EXT(I)
	    END DO
	  END IF
	  DO I=ND_ADD+1,ND_EXT
	    CHI_EXT(I)=CHI(I-ND_ADD)
	  END DO
C
C We limit alpha to 3.5 to avoid excess enevelope emission. If alpha were
C 3 we would get a logarithmic flux divergence as we increase the volume.
C
	  ALPHA=LOG(ETA(4)/ETA(1))/LOG(R(1)/R(4))
	  IF(ALPHA .LT. 3.5)ALPHA=3.5
	  DO I=1,ND_ADD
	    ETA_EXT(I)=ETA(1)*(R(1)/R_EXT(I))**ALPHA
	  END DO
	  DO I=ND_ADD+1,ND_EXT
	    ETA_EXT(I)=ETA(I-ND_ADD)
	  END DO
	ELSE
	  CHI_EXT(1:ND)=CHI(1:ND)		!NB: In this can ND=ND_EXT
	  ETA_EXT(1:ND)=ETA(1:ND)
	END IF
C
	SOURCE(1:ND_EXT)=ETA_EXT(1:ND_EXT)/CHI_EXT(1:ND_EXT)
C
C Perform monotonic interpolations in CHI and ETA. The dCHIR... arrays are
C not used as the derivatives are not computed (last passed variable).
C Note: We pass dCHIdR(ND_ADD+1) to MON_INT_INS as this will store the
C derivatives in the correct location for the EXTENED ray.
C
	CALL MON_INT_INS_V1(SOURCE_INS,R_INS,NINS,ETA,R,ND,
	1                   LTRUE,LTRUE,dCHIdR(ND_ADD+1),dCHIdR_INS,LFALSE)
	CALL MON_INT_INS_V1(CHI_INS,R_INS,NINS,CHI,R,ND,
	1                   LTRUE,LTRUE,dCHIdR(ND_ADD+1),dCHIdR_INS,LTRUE)
	SOURCE_INS(:,:)=SOURCE_INS(:,:)/CHI_INS(:,:)
C
	IF(METHOD .EQ. 'ZERO')THEN
	  dCHIdR(:)=0.0D0
	  dCHIdR_INS(:,:)=0.0D0
	END IF
C
C Zero boundary conditions.
C
	HBC=0.0D0			!H/J at model outer boundary.
	NBC=0.0D0			!N/J at model outer boundary.
	IN_HBC=0.0D0
C
C Zero AV and CV matrices.
C
	AV(:,:)=0.0D0
	CV(:,:)=0.0D0
C
C Enter loop to perform integration along each ray.
C
	DO LS=1,NP
	  NI=NI_RAY(LS)
	  NI_SMALL=ND-(LS-NC-1)
	  IF(LS .LE. NC+1)THEN
	     NI=ND_EXT
             NI_SMALL=ND
	     CHI_RAY(1:NI)=CHI_EXT(1:NI)
	     SOURCE_RAY(1:NI)=SOURCE(1:NI)
	     dCHIdR_RAY(1:NI)=dCHIdR(1:NI)
	  ELSE
C
C Include additional points to accomodate large change in UV near z=0.
C	
	    IF(NI_RAY(LS) .GT. ND_ADD+9)THEN
C
C In this case we insert 4 points between the last 2 points, and 2 points between 
C the second last pair.
C
	      J=NI-6
	      CHI_RAY(1:J)=CHI_EXT(1:J)
              CHI_RAY(J+6)=CHI_RAY(J)
              CHI_RAY(J+1)=CHI_RAY(J-1)
	      CHI_RAY(J-1)=CHI_INS(J-2-ND_ADD,1) 
	      CHI_RAY(J)=CHI_INS(J-2-ND_ADD,2) 
	      CHI_RAY(J+2)=CHI_INS(J-1-ND_ADD,1) 
	      CHI_RAY(J+3)=CHI_INS(J-1-ND_ADD,2) 
	      CHI_RAY(J+4)=CHI_INS(J-1-ND_ADD,3) 
	      CHI_RAY(J+5)=CHI_INS(J-1-ND_ADD,4)
C
	      SOURCE_RAY(1:J)=SOURCE(1:J)		!Set  NI, and NI-1
              SOURCE_RAY(J+6)=SOURCE_RAY(J)
              SOURCE_RAY(J+1)=SOURCE_RAY(J-1)
	      SOURCE_RAY(J-1)=SOURCE_INS(J-2-ND_ADD,1) 
	      SOURCE_RAY(J)=SOURCE_INS(J-2-ND_ADD,2) 
	      SOURCE_RAY(J+2)=SOURCE_INS(J-1-ND_ADD,1) 
	      SOURCE_RAY(J+3)=SOURCE_INS(J-1-ND_ADD,2) 
	      SOURCE_RAY(J+4)=SOURCE_INS(J-1-ND_ADD,3) 
	      SOURCE_RAY(J+5)=SOURCE_INS(J-1-ND_ADD,4)
C
	      dCHIdR_RAY(1:J)=dCHIdR(1:J)		!Set  NI, and NI-1
              dCHIdR_RAY(J+6)=dCHIdR_RAY(J)
              dCHIdR_RAY(J+1)=dCHIdR_RAY(J-1)
	      dCHIdR_RAY(J-1)=dCHIdR_INS(J-2-ND_ADD,1) 
	      dCHIdR_RAY(J)=dCHIdR_INS(J-2-ND_ADD,2) 
	      dCHIdR_RAY(J+2)=dCHIdR_INS(J-1-ND_ADD,1) 
	      dCHIdR_RAY(J+3)=dCHIdR_INS(J-1-ND_ADD,2) 
	      dCHIdR_RAY(J+4)=dCHIdR_INS(J-1-ND_ADD,3) 
	      dCHIdR_RAY(J+5)=dCHIdR_INS(J-1-ND_ADD,4)
C
	    ELSE IF(NI_RAY(LS) .EQ. ND_ADD+7)THEN
C
C In this case we insert 4 points between the last 2 points only (see earlier).
C
	      J=NI-4
	      CHI_RAY(1:J)=CHI_EXT(1:J)
              CHI_RAY(J+4)=CHI_RAY(J)
	      CHI_RAY(J)=CHI_INS(J-1-ND_ADD,1) 
	      CHI_RAY(J+1)=CHI_INS(J-1-ND_ADD,2) 
	      CHI_RAY(J+2)=CHI_INS(J-1-ND_ADD,3) 
	      CHI_RAY(J+3)=CHI_INS(J-1-ND_ADD,4)
C
	      SOURCE_RAY(1:J)=SOURCE(1:J)
              SOURCE_RAY(J+4)=SOURCE_RAY(J)
	      SOURCE_RAY(J)=SOURCE_INS(J-1-ND_ADD,1) 
	      SOURCE_RAY(J+1)=SOURCE_INS(J-1-ND_ADD,2) 
	      SOURCE_RAY(J+2)=SOURCE_INS(J-1-ND_ADD,3) 
	      SOURCE_RAY(J+3)=SOURCE_INS(J-1-ND_ADD,4)
C
	      dCHIdR_RAY(1:J)=dCHIdR(1:J)
              dCHIdR_RAY(J+4)=dCHIdR_RAY(J)
	      dCHIdR_RAY(J)=dCHIdR_INS(J-1-ND_ADD,1) 
	      dCHIdR_RAY(J+1)=dCHIdR_INS(J-1-ND_ADD,2) 
	      dCHIdR_RAY(J+2)=dCHIdR_INS(J-1-ND_ADD,3) 
	      dCHIdR_RAY(J+3)=dCHIdR_INS(J-1-ND_ADD,4)
C
	    END IF
	  END IF
C
C Even in the THK case we assume that the intensity incident on the
C outer boundary is zero. In THK case we effectively get IBOUND not
C equal zero at the model atmosphere boundary as a consequence of the
C extension.
C
	  IBOUND=0.0D0
C
C 
C
C By setting PF(1)=0 when evaluating SOURCE we ensure a pure continuum 
C calculation for the first frequency.
C
	  IF(INIT)THEN
	    DO I=1,NI
	      Q(I)=0.0D0
	      QH(I)=0.0D0
	    END DO                
	    OLDCHI(LS)=CHI_RAY(NI)
	  ELSE
	    DO I=1,NI-1
	      QH(I)=GAMH(I,LS)*2.0D0/((CHI_RAY(I)+CHI_RAY(I+1))*dLOG_NU)
	      Q(I)=GAM(I,LS)/(CHI_RAY(I)*dLOG_NU)
	    END DO
	    QH(NI)=0.0D0
	    Q(NI)=GAM(NI,LS)/(CHI_RAY(NI)*dLOG_NU)
	  END IF
	  IF(DIF .AND. LS .LE. NC)THEN
	      DBC=DBB*SQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)/CHI_RAY(NI)
	1   *(1.0D0+Q(NI)*(1.0D0-CHI_RAY(NI)/OLDCHI(LS)))
	  END IF
	  OLDCHI_STORE(LS)=CHI_RAY(NI)
C
C Compute the optical depth increments. This code is from TAU, and NORDTAU. We
C check that the Euler-Mauclarin correction is not too large. This is mainly
C done to prevent negative optical depths.
C
	  IF(METHOD .EQ. 'ZERO')THEN
	    DO I=1,NI-1
	      DTAU(I)=0.5D0*(CHI_RAY(I)+CHI_RAY(I+1))*(Z(I,LS)-Z(I+1,LS))
	    END DO
	  ELSE
	     DO I=1,NI-1
	        dZ=Z(I,LS)-Z(I+1,LS)
	        DTAU(I)=0.5D0*dZ*(CHI_RAY(I)+CHI_RAY(I+1))
	        T1=dZ*dZ*( dCHIdR_RAY(I+1)*Z(I+1,LS)/R_RAY(I+1,LS)-
	1              dCHIdR_RAY(I)*Z(I,LS)/R_RAY(I,LS) )/12.0D0
                T1=SIGN(  MIN( ABS(0.8D0*DTAU(I)),ABS(T1) ),T1  )
	        DTAU(I)=DTAU(I)+T1
	     END DO
	  END IF
	  DTAU_BND(LS)=DTAU(NI-1)
C
	  CALL TUVGHD_RH(TA(1,LS),TB(1,LS),TC(1,LS),
	1           U,VB,VC,GB(1,LS),H(1,LS),XM,
	1           Q,QH,DTAU,SOURCE_RAY,DIF,DBC,IC,LS,NC,NI)
	  XM(1)=-IBOUND				!Needs fixing
C
C Update AV matrix.
C
	  AV(1,LS)=XM(1)+U(1)*AV_PREV(1,LS)
	  DO I=2,NI-1
	      AV(I,LS)=XM(I)+( U(I)*AV_PREV(I,LS)-
	1             (VB(I)*CV_PREV(I-1,LS)+VC(I)*CV_PREV(I,LS)) )
	  END DO
	  AV(NI,LS)=XM(NI)+( U(NI)*AV_PREV(NI,LS)-VB(NI)*CV_PREV(NI-1,LS) )
C
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
C
	END DO				!LS
C
C 
C 
C Solve for the radiation field along ray for this frequency.
C
	IF(VECTOR_MACHINE)THEN
C
C This section is for a VECTOR machine. Loop over inner index is outer
C loop to remove a dependency. This in inefficient on scaler machines
C as array is not accessed sequentially.
C
C
C Forward substitution.
C
	  AV(1,1:NP)=AV(1,1:NP)*TB(1,1:NP)
	  DO I=2,NRAY_MAX
	    DO LS=1,MAX_LS(I)
	      IF(I .LE. NI_RAY(LS))
	1         AV(I,LS)=(AV(I,LS)+TA(I,LS)*AV(I-1,LS))*TB(I,LS)
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
	1        AV(I,LS)=TC(I,LS)*AV(I+1,LS)-AV(I,LS)
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
C
	DO LS=1,NP
C
C Verify validity of AV values. If these go negative, we set them +ve 
C but a factor of 10 smaller. We also ensure that the AV are not 
C extremely close to zero by comparing with neigboring AV values.
C The CV (fluxes) are computed using the revised values.
C
	    K=2
	    DO I=1,NI_RAY(LS)
	      T1=1.0D-06*ABS(AV(K,LS))
	      IF(AV(I,LS) .LE. T1)THEN
	        AV(I,LS)=MAX(0.1D0*ABS(AV(I,LS)),T1)
	        NEG_AV_VALUE=.TRUE.
	      ELSE
	        K=K+1
	      END IF
	      IF(I .EQ. 1)K=1
	    END DO
C
C Update C vector (i.e. flux variable).
C
	    NI=NI_RAY(LS)
	    NI_SMALL=ND-(LS-NC-1)
	    IF(LS .LE. NC+1)NI_SMALL=ND
C                          
	    DO I=1,NI_RAY(LS)-1
	      CV(I,LS)=GB(I,LS)*(AV(I,LS)-AV(I+1,LS))+H(I,LS)*CV_PREV(I,LS)
	    END DO
C
C We have to be careful here as we have to wory about those points that
C have been inserted.
C
	  IF(LS .LE. NC+1)THEN
C
C No insertions were done for the core rays.
C
	    JNU(1:ND)=JNU(1:ND)+JQW(1:ND,LS)*AV(1+ND_ADD:ND+ND_ADD,LS)
	    KNU(1:ND)=KNU(1:ND)+KQW(1:ND,LS)*AV(1+ND_ADD:ND+ND_ADD,LS)
	    J=ND-1
	    HNU(1:J)=HNU(1:J)+HQW(1:J,LS)*CV(1+ND_ADD:J+ND_ADD,LS)
	    NNU(1:J)=NNU(1:J)+NQW(1:J,LS)*CV(1+ND_ADD:J+ND_ADD,LS)
C
	  ELSE
	    DO I=1,NI_SMALL-2
	      JNU(I)=JNU(I)+JQW(I,LS)*AV(I+ND_ADD,LS)
	      KNU(I)=KNU(I)+KQW(I,LS)*AV(I+ND_ADD,LS)
	    END DO
	    IF(NI .GE. ND_ADD+7)THEN
	      JNU(NI_SMALL-1)=JNU(NI_SMALL-1)+JQW(NI_SMALL-1,LS)*AV(NI-5,LS)
	      KNU(NI_SMALL-1)=KNU(NI_SMALL-1)+KQW(NI_SMALL-1,LS)*AV(NI-5,LS)
	      JNU(NI_SMALL)=JNU(NI_SMALL)+JQW(NI_SMALL,LS)*AV(NI,LS)
	      KNU(NI_SMALL)=KNU(NI_SMALL)+KQW(NI_SMALL,LS)*AV(NI,LS)
	    ELSE IF(NI .EQ. ND_ADD+2)THEN
	      JNU(NI_SMALL-1)=JNU(NI_SMALL-1)+JQW(NI_SMALL-1,LS)*AV(NI-1,LS)
	      KNU(NI_SMALL-1)=KNU(NI_SMALL-1)+KQW(NI_SMALL-1,LS)*AV(NI-1,LS)
	      JNU(NI_SMALL)=JNU(NI_SMALL)+JQW(NI_SMALL,LS)*AV(NI,LS)
	      KNU(NI_SMALL)=KNU(NI_SMALL)+KQW(NI_SMALL,LS)*AV(NI,LS)
	    ELSE
	      IF(NI_SMALL .NE. 1)THEN
	        J=ERROR_LU()
	        WRITE(J,*)'Error in FG_J_CMF_V6 --- invalid code logic'
	        STOP
	      END IF
	      JNU(NI_SMALL)=JNU(NI_SMALL)+JQW(NI_SMALL,LS)*AV(NI,LS)
	      KNU(NI_SMALL)=KNU(NI_SMALL)+KQW(NI_SMALL,LS)*AV(NI,LS)
	    END IF
C
	    DO I=1,NI_SMALL-3
	      HNU(I)=HNU(I)+HQW(I,LS)*CV(I+ND_ADD,LS)
	      NNU(I)=NNU(I)+NQW(I,LS)*CV(I+ND_ADD,LS)
	    END DO
	    IF(NI .GE. ND_ADD+9)THEN
	      J=NI_SMALL-2
	      HNU(J)=HNU(J)+HQW(J,LS)*CV(NI-7,LS)
	      NNU(J)=NNU(J)+NQW(J,LS)*CV(NI-7,LS)
	      J=NI_SMALL-1
	      HNU(J)=HNU(J)+HQW(J,LS)*CV(NI-4,LS)
	      NNU(J)=NNU(J)+NQW(J,LS)*CV(NI-4,LS)
	    ELSE IF(NI .EQ. ND_ADD+7)THEN
	      J=NI_SMALL-2
	      HNU(J)=HNU(J)+HQW(J,LS)*CV(NI-6,LS)
	      NNU(J)=NNU(J)+NQW(J,LS)*CV(NI-6,LS)
	      J=NI_SMALL-1
	      HNU(J)=HNU(J)+HQW(J,LS)*CV(NI-4,LS)
	      NNU(J)=NNU(J)+NQW(J,LS)*CV(NI-4,LS)
	    ELSE IF(NI .EQ. ND_ADD+2)THEN
	      J=NI_SMALL-1
	      HNU(J)=HNU(J)+HQW(J,LS)*CV(NI-1,LS)
	      NNU(J)=NNU(J)+NQW(J,LS)*CV(NI-1,LS)
	    END IF
	  END IF
	
C
C The quadrature weights used to compute HBC and NBC are now reliable.
C Note that V=AV(1)-IBOUND.
C
	  IF(ND_ADD .EQ. 0)THEN
	     CV_BOUND=AV(1,LS)-IBOUND
	  ELSE
	    IF(NI_SMALL .EQ. 1)THEN
	      CV_BOUND=0.0D0
	    ELSE
	      CV_BOUND=0.5D0*(CV(ND_ADD+1,LS)+CV(ND_ADD,LS))
	    END IF
	  END IF
	  HBC=HBC+CV_BOUND*H_WGHTS(LS)
	  NBC=NBC+CV_BOUND*N_WGHTS(LS)
C          
	  IF(LS .LE. NC+1)THEN
	    IN_HBC=IN_HBC + 
	1        JQW(ND,LS)*( AV(NI,LS)-(AV(NI,LS)-AV(NI-1,LS))
	1          /DTAU_BND(LS) )*Z(NI,LS)/R_RAY(NI,LS)
 	  END IF
C
	  DO I=1,NI
	   AV_STORE(I,LS)=AV(I,LS)
	   CV_STORE(I,LS)=CV(I,LS)
	  END DO
	  IPLUS_P(LS)=2.0D0*CV_BOUND          !To compute observed flux
C 
	END DO			!End do LS
C
C
C
C Compute the boundary Eddington factors.
C
	HBC=HBC/JNU(1)
	NBC=NBC/JNU(1)
	IN_HBC=IN_HBC/(2.0D0*JNU(ND)-IC)
C
C Compute the Eddington F and G factors, which are defined by
C F=K/J and G=N/H. The F and G factors are returned in KNU and NNU
C respectively.
C
	DO I=1,ND
	  KNU(I)=KNU(I)/JNU(I)
	END DO
C
C Now compute the "G" Eddington factors. Because H may be zero we have to
C be careful. Depending on N_TYPE we can specify N in terms of J, N in
C terms of H, or N in terms of J and H.
C
	IF(N_TYPE .EQ. 'N_ON_J')THEN
	  DO I=1,ND-1
	    ALPHA=0.25D0*(R(I)+R(I+1))*(R(I)+R(I+1))
	    RSQN_ON_RSQJ(I)=ALPHA*NNU(I)/
	1                  (R(I)*R(I)*JNU(I)+R(I+1)*R(I+1)*JNU(I+1))
	    NNU(I)=0.0D0
	  END DO
	ELSE IF(N_TYPE .EQ. 'MIXED')THEN
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
C
C Compute G Eddington factor storing in N.
C
	  DO I=1,ND-1
	    IF(HNU(I) .NE. 0)THEN
	      NNU(I)=NNU(I)/HNU(I)
	    ELSE
	      NNU(I)=0.0D0		!Later replaced by average.
	      J=ERROR_LU()
	      WRITE(J,'(1X,A,1PE16.8)')'HNU zero for frequency:',FL
	      WRITE(J,'(1X,A,I4)')'Error occurred at depth:',I
	    END IF
	  END DO
C
C Check the validity of the Eddington factor in case strange values are
C occurring because H is near zero (switching sign?).
C
	  DO I=2,ND-1
	    IF(NNU(I) .GT. 1.10)THEN
              BAD_NNU_VALUE=.TRUE.
	    ELSE IF( NNU(I) .LT. 0.05) THEN
	      NEG_NNU_VALUE=.TRUE.
	    END IF
	  END DO
	  DO I=2,ND-1
	    IF(NNU(I) .GT. 1.1 .OR. NNU(I) .LT. 0.05)THEN
	      J=1
	      DO WHILE( ( (NNU(I+J) .GT. 1.0) .OR. (NNU(I+J) .LT. 0.05) )
	1       .AND. (J .LT. ND-1) )
	        J=J+1
	      END DO
	      DO K=I,I+J-1
	        NNU(K)=NNU(I-1)+(K-I+1.0D0)/(J+1.0D0)*(NNU(J)-NNU(I-1))
	      END DO
	    END IF
	  END DO
	ELSE
	  J=ERROR_LU()
	  WRITE(J,*)'Unrecognized N_TYPE in FG_J_CMF_V4'
	  WRITE(J,*)'N_TYPE=',N_TYPE
	  STOP
	END IF
C
C Store frequencies at which errors occurred, and give an indication of the
C error.
C
	J=1
	DO WHILE (J .LE. N_ERR_MAX .AND. 
	1            (NEG_AV_VALUE .OR. BAD_NNU_VALUE .OR. NEG_NNU_VALUE) )
	  IF(FG_ERR_ON_FREQ(J) .EQ. FL)THEN
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
	    FG_ERR_ON_FREQ(J)=FL
	    IF(NEG_AV_VALUE)FG_ERR_TYPE(J)=FG_ERR_TYPE(J)+1
	    IF(NEG_NNU_VALUE)FG_ERR_TYPE(J)=FG_ERR_TYPE(J)+5
	    IF(BAD_NNU_VALUE)FG_ERR_TYPE(J)=FG_ERR_TYPE(J)+10
	    NEG_AV_VALUE=.FALSE.
	    BAD_NNU_VALUE=.FALSE.
	    NEG_NNU_VALUE=.FALSE.
	  END IF
	  J=J+1
	END DO
C
	DO I=1,ND
	  IF(JNU(I) .LT. 0)THEN
	    DO J=1,ND
	      WRITE(111,*)FL
	      WRITE(111,*)JNU(I)
	    END DO
	    STOP
	  END IF
	END DO
C
	RETURN
	END
 
