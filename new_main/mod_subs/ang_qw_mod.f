!
! Data module defining Quadrature Weights, and related variables,  
! for CMFGEN. 
!
	MODULE ANG_QW_MOD
!
! Created 02-May-2004
!
! Angular qudrature weights for J (AQW), H, and K. These are
! frequency indepoendent.
!
	REAL*8, ALLOCATABLE :: P(:)              !NP
	REAL*8, ALLOCATABLE :: AQW(:,:)          !ND,NP
	REAL*8, ALLOCATABLE :: HQW(:,:)          !ND,NP
	REAL*8, ALLOCATABLE :: KQW(:,:)          !ND,NP
	REAL*8, ALLOCATABLE :: NQW(:,:)          !ND,NP
!
! Defined at the mid points of the mesh.
!
	REAL*8, ALLOCATABLE :: HMIDQW(:,:)       !ND,NP
	REAL*8, ALLOCATABLE :: NMIDQW(:,:)       !ND,NP
!
! If required, these arrays shouLd have size NDEXT*NPEXT
!
	REAL*8, ALLOCATABLE :: PEXT(:)            !NPEXT
	REAL*8, ALLOCATABLE :: AQWEXT(:,:)
	REAL*8, ALLOCATABLE :: HQWEXT(:,:)
	REAL*8, ALLOCATABLE :: KQWEXT(:,:)
	REAL*8, ALLOCATABLE :: NQWEXT(:,:)
	REAL*8, ALLOCATABLE :: HMIDQWEXT(:,:)
	REAL*8, ALLOCATABLE :: NMIDQWEXT(:,:)
!
! Parameters, vectors, and arrays for computing the observed flux.
! Was 2000 (changed 20-Apr-2009).
!
        INTEGER, PARAMETER :: NST_CMF=6000
        REAL*8  NU_STORE(NST_CMF)
!
        INTEGER NP_OBS_MAX
        INTEGER NP_OBS
        REAL*8 V_AT_RMAX                !Used if we extend the atmosphere.
        REAL*8 RMAX_OBS
        REAL*8 HFLUX_AT_OB		!In comoving frame
        REAL*8 HFLUX_AT_IB
!
! We allocate memory for the following vectors as we use them for the regular
! flux computation, and when extra depth points are inserted (ACCURATE=.TRUE.)
!
        REAL*8, ALLOCATABLE :: IPLUS_STORE(:,:)
        REAL*8, ALLOCATABLE :: P_OBS(:)
        REAL*8, ALLOCATABLE :: IPLUS(:)
        REAL*8, ALLOCATABLE :: MU_AT_RMAX(:)
        REAL*8, ALLOCATABLE :: HQW_AT_RMAX(:)
!
	END MODULE ANG_QW_MOD
