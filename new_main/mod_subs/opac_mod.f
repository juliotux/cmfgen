	MODULE OPAC_MOD
	IMPLICIT NONE
!
	REAL*8, ALLOCATABLE :: ETA(:)
	REAL*8, ALLOCATABLE :: ETA_CONT(:)
	REAL*8, ALLOCATABLE :: ETA_C_EVAL(:)
	REAL*8, ALLOCATABLE :: ETA_NOSCAT(:)
	REAL*8, ALLOCATABLE :: ETA_NOSCAT_EVAL(:)
	REAL*8, ALLOCATABLE :: ETA_CLUMP(:)
	REAL*8, ALLOCATABLE :: ETA_MECH(:)
!
	REAL*8, ALLOCATABLE :: CHI(:)
	REAL*8, ALLOCATABLE :: CHI_CONT(:)
	REAL*8, ALLOCATABLE :: CHI_C_EVAL(:)
	REAL*8, ALLOCATABLE :: CHI_NOSCAT(:)
	REAL*8, ALLOCATABLE :: CHI_NOSCAT_EVAL(:)
	REAL*8, ALLOCATABLE :: CHI_CLUMP(:)
	REAL*8, ALLOCATABLE :: TCHI(:)
!
	REAL*8, ALLOCATABLE :: CHI_SCAT(:)
	REAL*8, ALLOCATABLE :: CHI_RAY(:)
	REAL*8, ALLOCATABLE :: ESEC(:)
	REAL*8, ALLOCATABLE :: ESEC_CLUMP(:)
	REAL*8, ALLOCATABLE :: CHI_SCAT_CLUMP(:)
!
! Use to relate opacities at the current frequency to that at last frequency.
! They store the opacities at the previous frequency. Note that CHI_SCAT is
! only constant when we have pure electron scattering --- it varies if we allow
! for Rayleigh scattering.
!
	REAL*8, ALLOCATABLE ::  CHI_PREV(:)
	REAL*8, ALLOCATABLE ::  CHI_NOSCAT_PREV(:)
	REAL*8, ALLOCATABLE ::  CHI_SCAT_PREV(:)
	REAL*8, ALLOCATABLE ::  ETA_PREV(:)
!
! To allow the variation of non-coherent electron scattering to be treated
! in a partially coherent approximation.
!
	REAL*8, ALLOCATABLE :: ES_COH_VEC(:)
	REAL*8, ALLOCATABLE :: THETA(:)
	REAL*8, ALLOCATABLE :: ZETA(:)
	REAL*8, ALLOCATABLE :: SOURCE(:)
!
	REAL*8, ALLOCATABLE :: EMHNUKT(:)
	REAL*8, ALLOCATABLE :: EMHNUKT_CONT(:)
!
	REAL*8, ALLOCATABLE :: XRAY_LUM_0P1(:)
	REAL*8, ALLOCATABLE :: XRAY_LUM_1KEV(:)
	REAL*8, ALLOCATABLE :: XRAY_LUM_TOT(:)
!
	REAL*8, ALLOCATABLE :: VCHI(:,:)              !Variation of CHI array.
        REAL*8, ALLOCATABLE :: VETA(:,:)              !Variation of ETA array.
        REAL*8, ALLOCATABLE :: VCHI_SAV(:,:)
        REAL*8, ALLOCATABLE :: VETA_SAV(:,:)
        REAL*8, ALLOCATABLE :: VCHI_ALL(:,:)
	REAL*8, ALLOCATABLE :: VCHI_ALL_SAV(:,:)
        REAL*8, ALLOCATABLE :: VETA_ALL(:,:)
	REAL*8, ALLOCATABLE :: VETA_ALL_SAV(:,:)
!
	END MODULE OPAC_MOD
