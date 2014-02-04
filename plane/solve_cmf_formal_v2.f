!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      SUBROUTINE SOLVE_CMF_FORMAL_V2(CHI,ETA,IP,FREQ,NU_DNU,
     *                     INNER_BND_METH,B_NUE,dBDTAU,ND,NP,NC)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Determine the transfer matrix along the given ip ray.  Formulation
! is the Short Characteristic solution from Olson & Kunsasz (87) JQSRT
! and Hauschildt (92) JQSRT
!
! Written   11-94 DLM
! Altered 2/21/96 DLM Updated to F90 standard
! Altered 5/9/97  DLM installed MODULE for rpz and allocation
!                     of many variables
! Altered 5/16/97 DLM installed three possible inner boundary conditions
!                     'HOLLOW'    - hollow core with I+=I-
!                     'GRAY'      - DIFFUSION approximation for gray atmosphere: Mihalas (1980) 237 p 574
!                     'DIFFUSION' - DIFFUSION approximation with I=B(nu,T)+mu*dBdtau(nu,T): Mihalas (1978) p 51
! Altered 6/12/97 DLM Name changed from solve_rel_formal.f
! Altered 31/12/04 DJH : Cleaned
! Altered 22-Jan-2010: Changed call to V3 and altered boundary condition options:
!                     'ZERO_FLUX' - hollow core with I+=I-
!                     'HOLLOW'    - hollow core with I+=I- but with allowance for velocity shifts.
!                     'GRAY'      - DIFFUSION approximation for gray atmosphere: Mihalas (1980) 237 p 574
!                     'DIFFUSION' - DIFFUSION approximation with I=B(nu,T)+mu*dBdtau(nu,T): Mihalas (1978) p 51
! Altered 09-May-2010: Changed to call OPTDEPTH_V4.
!                      CUR_LOC is noe set in calling routine. Since independent of ip, allows this routine to be
!                        in parallel loop.
! Altered 07-Apr-2103: Source function is now interpolated using linear interpolation for last two intervals for
!                        outgoing rays. Problems were arising because of rapidly varying source functon and very 
!                        unequal step sizes.
!
!--------------------------------------------------------------------
!
      USE MOD_SPACE_GRID_V2
      IMPLICIT NONE
!
! Grid size variables
!
      INTEGER :: NP,ND,NC
      INTEGER :: I,IP
!
! Opacity and emissivity variables
!
      real*8, dimension(nd) :: chi,eta
!
! Frequency variable
!
      real*8 :: freq
      real*8 :: nu_dnu
!
! Boundary conditions
!
      real*8 :: B_nue,dBdtau
!
      character(len=*) :: INNER_BND_METH
!
! Local variables
!
      LOGICAL, PARAMETER :: L_TRUE=.TRUE.
      LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
      real*8 ibound
      real*8 ee,e0,e1,e2,alpha,beta,gamma,t1
      real*8, dimension(nd) :: chi_tau
      real*8, dimension(nd) :: source_prime
      real*8, dimension(nd) :: tau_loc
      real*8, dimension(nd) :: dtau_loc
      real*8 new_freq
!
      integer ist,iend,imid
      integer k
!
      integer :: iz
      integer :: nzz
      integer :: luer,error_lu
      external error_lu
!
!--------------------------------------------------------------------
!
! Determine transfer variables.
!
!   Returns:
!     source_prime=(eta+alpha*Iprev)/(alpha+chi_prime)
!     chi_tau=alpha+chi_prime
!
!     where:
!       chi_prime=chi+3*b_*
!       alpha=freq(k)*b_*/(freq(k-1)-freq(k))
!       b_*=gamma*((1-mu_*^2)*beta/r+gamma^2*mu_**(mu_*+beta)dbetadr)
!
!         where: _* is _p or _m for ray in plus or minus direction
!
!	WRITE(6,*)ND,NU_DNU
!	WRITE(6,*)IP,RAY(IP)%NZ
!	WRITE(6,*)ETA(1),CHI(1)
!	WRITE(6,*)RAY(IP)%B_M(1)
!	WRITE(6,*)RAY(IP)%I_M_PREV(1)
!	WRITE(6,*)CHI_TAU(1)
!	WRITE(6,*)SOURCE_PRIME(1)
!
      CALL REL_VARIABLES(ND,CHI,ETA,NU_DNU,RAY(IP)%B_M,RAY(IP)%I_M_PREV,CHI_TAU,SOURCE_PRIME)
!
! Determine "optical depth" from outside to center (s_m,mu_m)
!
      CALL OPTDEPTH_V4(DTAU_LOC,CHI_TAU,RAY(IP)%NZ,IP,L_FALSE,'ZERO')
!
! Calculate transfer in inward direction, I- (mu=-1)
!
      RAY(IP)%I_m(1)=0.0d0
!
!---------------------------------------------------------------
!
! Formal integral from OUTSIDE to INSIDE
!
!---------------------------------------------------------------
!
      nzz=ray(ip)%nz
      do iz=2,ray(ip)%nz-1
!
        t1=dtau_loc(iz-1)
        if(t1 .gt. 0.01D0)then
          ee=exp(-dtau_loc(iz-1))
          e0=1.0d0-ee
          e1=dtau_loc(iz-1)-e0
          e2=dtau_loc(iz-1)*dtau_loc(iz-1)-2.0d0*e1
        else
          e2=t1*t1*t1*(1.0D0-0.25D0*t1*(1.0d0-0.2D0*t1*(1.0D0-t1/6.0D0*
     1              (1.0D0-t1/7.0D0))))/3.0D0
          e1=0.5d0*(t1*t1-e2)
          e0=t1-e1
          ee=1.0d0-e0
        end if
!
        alpha=e0+(e2-(dtau_loc(iz)+2.0d0*dtau_loc(iz-1))*e1)/(dtau_loc(iz-1)*(dtau_loc(iz)+dtau_loc(iz-1)))
        beta=((dtau_loc(iz)+dtau_loc(iz-1))*e1-e2)/(dtau_loc(iz-1)*dtau_loc(iz))
        gamma=(e2-dtau_loc(iz-1)*e1)/(dtau_loc(iz)*(dtau_loc(iz)+dtau_loc(iz-1)))
     	t1=alpha*source_prime(iz-1)+beta*source_prime(iz)+gamma*source_prime(iz+1)
        if(t1 .lt. 0.0d0)then
	  t1=source_prime(iz)+dtau_loc(iz)*(source_prime(iz)-source_prime(iz-1))/dtau_loc(iz-1)
     	  t1=alpha*source_prime(iz-1)+beta*source_prime(iz)+gamma*t1
	end if
	ray(ip)%I_m(iz)=ray(ip)%I_m(iz-1)*ee+t1
!
      end do
!
! If only 2 points along the ray or at inner boundary, must use linear
! source function
!
      t1=dtau_loc(nzz-1)
      if(t1 .gt. 0.01d0)then
        ee=exp(-dtau_loc(nzz-1))
        e0=1.0d0-ee
        e1=dtau_loc(nzz-1)-e0
        e2=dtau_loc(nzz-1)*dtau_loc(nzz-1)-2.0d0*e1
      else
        e2=t1*t1*t1*(1.0D0-0.25D0*t1*(1.0D0-0.2D0*t1*(1.0D0-t1/6.0D0*
     1              (1.0D0-t1/7.0D0))))/3.0D0
        e1=0.5d0*(t1*t1-e2)
        e0=t1-e1
        ee=1.0d0-e0
      end if
!
      alpha=e0-e1/dtau_loc(nzz-1)
      beta=e1/dtau_loc(nzz-1)
      ray(ip)%I_m(nzz)=ray(ip)%I_m(nzz-1)*ee+
     *     alpha*source_prime(nzz-1)+beta*source_prime(nzz)
!
! If we are using the holow boundary condition, we need to store the
! intensity at the inner boundary. This flux will be used to compute
! the incident intensity at another frequency. We can't directly use
! the current flux, since the radiation will be redshifted by the 
! expansion.
!
      if(ip .le. nc .and. INNER_BND_METH .eq. 'HOLLOW')then
	ray(ip)%i_in_bnd_store(cur_loc)=ray(ip)%I_m(nzz)
!
! Now get flux incident at inner boudary from other side of envelope.
!
	ist=cur_loc+1
        iend=cur_loc+n_store
        if(freq_store(n_store-1) .EQ. 0.0D0)THEN
          ist=0; iend=cur_loc
        end if
        k=mod(ist,n_store)
        new_freq=freq*ray(ip)%freq_conv_fac
	if(new_freq .ge. freq_store(k) .and. freq_store(n_store-1) .eq. 0.0D0)THEN
          ibound=ray(ip)%i_in_bnd_store(0)
	else if(new_freq .GT. freq_store(k))then
          write(6,*)'Error in solve_cmf_formal_v2: invalid frequency range'
	  write(6,*)freq_store(k),new_freq,freq
	  write(6,*)k,ist,iend
          stop
        else
          do while( iend-ist .gt. 1)
            imid=(ist+iend)/2
            k=mod(imid,n_store)
            if(new_freq .gt. freq_store(k))then
              iend=imid
            else
              ist=imid
            end if
          end do
          iend=mod(iend,n_store); ist=mod(ist,n_store)
          t1=(new_freq-freq_store(iend))/(freq_store(ist)-freq_store(iend))
          ibound=t1*ray(ip)%i_in_bnd_store(ist)+(1.0d0-t1)*ray(ip)%i_in_bnd_store(iend)
        end if
      end if
!
!---------------------------------------------------------------
!
! Calculate transfer in outward direction, I+ (mu=1)
!
! At inner boundary use I+ = I- for core and non-core rays
!
      if(ip .gt. nc+1 .or. INNER_BND_METH .eq. 'ZERO_FLUX')then
        ray(ip)%I_p(nzz)=ray(ip)%I_m(nzz)
      else if(INNER_BND_METH .eq. 'HOLLOW')then
        ray(ip)%I_p(nzz)=ibound
      elseif((INNER_BND_METH .eq. 'GRAY') .or. (INNER_BND_METH .eq. 'DIFFUSION'))then
        ray(ip)%I_p(nzz)=B_nue+ray(ip)%mu_p(nzz)*dBdtau
      else
	luer=error_lu()
        write(luer,*)'Unknown inner boundary condition in SOLVE_CMF_FORMAL'
        write(luer,*)'Passed boundary condition is ',trim(INNER_BND_METH)
        stop
      end if
!
! Determine transfer variables.
!
      call rel_variables(nd,chi,eta,nu_dnu,ray(ip)%b_p,
     *       ray(ip)%I_p_prev,chi_tau,source_prime)
!
! Determine "optical depth" from outside to center (s_m,mu_m)
!
      CALL OPTDEPTH_V4(DTAU_LOC,CHI_TAU,RAY(IP)%NZ,IP,L_TRUE,'ZERO')
!
!---------------------------------------------------------------
!
! Formal integral from INSIDE to OUTSIDE
!
!---------------------------------------------------------------
!
      do iz=nzz-1,3,-1
!
        t1=dtau_loc(iz)
        if(t1 .gt. 0.01d0)then
          ee=exp(-dtau_loc(iz))
          e0=1.0d0-ee
          e1=dtau_loc(iz)-e0
          e2=dtau_loc(iz)*dtau_loc(iz)-2.0d0*e1
        else
           e2=t1*t1*t1*(1.0D0-0.25D0*t1*(1.0D0-0.2D0*t1*(1.0D0-t1/6.0D0*
     1              (1.0D0-t1/7.0D0))))/3.0D0
           e1=0.5d0*(t1*t1-e2)
           e0=t1-e1
           ee=1.0d0-e0
         end if
!
        alpha=(e2-dtau_loc(iz)*e1)/
     *       (dtau_loc(iz-1)*(dtau_loc(iz)+dtau_loc(iz-1)))
        beta=((dtau_loc(iz)+dtau_loc(iz-1))*e1-e2)/(dtau_loc(iz-1)*dtau_loc(iz))
        gamma=e0+(e2-(dtau_loc(iz-1)+2.0d0*dtau_loc(iz))*e1)/(dtau_loc(iz)*(dtau_loc(iz)+dtau_loc(iz-1)))
        t1=alpha*source_prime(iz-1)+beta*source_prime(iz)+gamma*source_prime(iz+1)
        if(t1 .lt. 0.0d0)then
	  t1=source_prime(iz)+dtau_loc(iz-1)*(source_prime(iz)-source_prime(iz+1))/dtau_loc(iz)
     	  t1=alpha*t1+beta*source_prime(iz)+gamma*source_prime(iz+1)
	end if
        ray(ip)%I_p(iz)=ray(ip)%I_p(iz+1)*ee+t1
!
      end do
!
! If only 2 points along the ray or at outer boundary, must use linear
! source function
!
      do iz=min(nzz-1,2),1,-1
        t1=dtau_loc(iz)
        if(t1 .gt. 0.01d0)then
          ee=exp(-dtau_loc(iz))
          e0=1.0d0-ee
          e1=dtau_loc(iz)-e0
          e2=dtau_loc(iz)*dtau_loc(iz)-2.0d0*e1
        else
          e2=t1*t1*t1*(1.0D0-0.25D0*t1*(1.0d0-0.2D0*t1*(1.0D0-t1/6.0D0*
     1                (1.0D0-t1/7.0D0))))/3.0D0
          e1=0.5d0*(t1*t1-e2)
          e0=t1-e1
          ee=1.0d0-e0
        end if
!
        beta=e1/dtau_loc(iz)
        gamma=e0-e1/dtau_loc(iz)
        ray(ip)%I_p(iz)=ray(ip)%I_p(iz+1)*ee+beta*source_prime(iz)+gamma*source_prime(iz+1)
       end do
!     if(ip .le. nc)write(166,*)freq,ray(ip)%I_p(nzz),ray(ip)%I_m(nzz)
!
      return
      end
