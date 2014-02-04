C
C Subroutine to correct for HI and HII absorption in the UV
C Is called by PLT_SPEC as an option, and in turn calls HIABS and HIIABS
C
	subroutine uvabs(origfreq,origflux,Norig,Norig_max,
     1                   t_in_k,v_turb,log_ntot,log_h2_ntot,
     1                   v_r,min_res_kms,wave_max,wave_min,
     1                   hi_abs,h2_abs)
C
	implicit none

	include 'constants.inc'

	integer*4 IOS             ! used for Input/Output errors.

	integer*4 Norig           ! number of points in initial array
	integer*4 Norig_max       ! Maximum number of points in initial array
!
	logical hi_abs             ! correct for HI absorption?
	logical h2_abs             ! correct for HII absorption?
C
	integer*4 i,j
	integer*4 Nmod            ! number of points in model array
	integer*4 Ntmp            ! number of points in temporary array
C
	integer*4 Norig_new
C
	real*8 origwave(Norig_max)
	real*8 origflux(Norig_max)
	real*8 origfreq(Norig_max)
!
	real*8, allocatable ::  modwave(:)
	real*8, allocatable ::  modflux(:)
	real*8, allocatable ::  tempwave(:)
	real*8, allocatable ::  tempflux(:)
	real*8, allocatable ::  Hwave(:)
	real*8, allocatable ::  Hflux(:)
	real*8, allocatable ::  Hfreq(:)
!
	real*8 t1
C
C       HIABS variables
C
	real*8 v_turb        ! turbulent velocity in km/s
	real*8 v_r           ! radial velocity of star w.r.t. ISM
	real*8 log_ntot      ! log base 10 of H column density
	real*8 log_h2_ntot   ! log base 10 of HII column density
	real*8 t_in_k        ! temperature in K
	real*8 wave_max      ! maximum wavelength
	real*8 wave_min      ! minimum wavelength
	real*8 hut_res       ! resolution of H.U.T. in angstroms
	real*8 model_res     ! model resolution in angstroms
	real*8 kernal_res    ! kernal to convolve model spectrm with res.
	real*8 min_res_kms   ! minimum resolution to work with
	logical fft          ! 1 = convolve via fft, 0 = straight convolution
	integer*4 alter_min  ! min and max indices of original array to be
        integer*4 alter_max  ! changed.
	integer*4 larger_index,smaller_index
!
	call nu_to_lambda(origfreq,origwave,Norig,forward)

	alter_min = larger_index(origwave,Norig,wave_min)
	alter_max = smaller_index(origwave,Norig,wave_max)
	Nmod = alter_max - alter_min + 1

	allocate (modwave(Nmod))
	allocate (modflux(Nmod))
	modwave(1:Nmod) = origwave(alter_min:alter_max) ! preserve original 
	modflux(1:Nmod) = origflux(alter_min:alter_max) ! arrays for later

C       Here we look at the model data to determine the minimum spacing (in
C	lambda) to get our model_res.  This is presumably the spacing found
C	near the lines.  However, the user can dictate the minimum resolution
C       to use to keep the array from getting too large.

	call get_model_res(modwave,nmod,model_res)
	t1=wave_min*min_res_kms/c_kms
	if (model_res .lt. t1)model_res = t1
	t1=0.25D0*wave_min*v_turb/c_kms
	if (model_res .gt. t1)model_res = t1

	Ntmp=(modwave(Nmod)-modwave(1))/model_res+1
	Ntmp=max(Ntmp,NMOD)
	allocate (tempwave(Ntmp))
	allocate (tempflux(Ntmp))
	write(6,*)'Begnning linearization'
        call linearize(tempwave,tempflux,Ntmp,
	1                modwave,modflux,Nmod,model_res)
	write(6,*)'End linearization'
	if (v_r .ne. zero) then       ! doppler shift model spectrum
         tempwave(1:Ntmp) = (one + v_r/c_kms)*tempwave(1:Ntmp)
	end if

C       set up original H absorption arrays and calculate HI and HII 
C       absorption fluxes.  Here we constrain the H array to have the 
C       same spacing in lambda as the model data.

	allocate (Hwave(Ntmp))
	allocate (Hfreq(Ntmp))
	allocate (Hflux(Ntmp))
	call Habsinit(Hwave,Hflux,Hfreq,tempwave,Ntmp)
	if (hi_abs) then
	   call HIabs(Hwave,Hflux,Ntmp,v_turb,log_ntot,t_in_k)
	endif
	if (h2_abs) then
	   call H2abs(Hwave,Hflux,Ntmp,v_turb,log_h2_ntot,t_in_k)
	endif
	call is_line_abs(Hwave,Hflux,Ntmp,v_turb,log_h2_ntot,t_in_k)
!
	write(6,*)'Done IS abs'

C       multiply fluxes

	tempflux(1:Ntmp) = tempflux(1:Ntmp)*Hflux(1:Ntmp)

! store in original model array

	if (v_r .ne. zero) then       ! doppler shift back 
         tempwave(1:Ntmp) = tempwave(1:Ntmp)/(one + v_r/c_kms)
	end if

C      Remap modified array back onto original wavelength array
!
	Norig_new=Norig-(alter_max-alter_min+1)+Ntmp
	if(Norig_new .gt. Norig_max)then
	  write(6,*)'Error in UVABS --- input array too small'
	  write(6,*)'Norig_new=',Norig_new
	  write(6,*)'Norig_max=',Norig_max
	  goto 999
	end if
	IF(ALTER_MAX .NE. Norig)THEN
	  origflux(alter_min+Ntmp:Norig_new)=origflux(alter_max+1:Norig)
	  origfreq(alter_min+Ntmp:Norig_new)=origfreq(alter_max+1:Norig)
	END IF
	tempwave(1:Ntmp)=1.0D-02*C_KMS/tempwave(1:Ntmp)
	origflux(alter_min:alter_min+Ntmp-1)=tempflux(1:Ntmp)
	origfreq(alter_min:alter_min+Ntmp-1)=tempwave(1:Ntmp)
	Norig=Norig_new
!
!	call map(tempwave,tempflux,Ntmp,modwave,modflux,Nmod)
!	origflux(alter_min:alter_max) = modflux(1:Nmod)
!
! Deallocate all memory
!
999	continue
	deallocate (Hwave)
	deallocate (Hfreq)
	deallocate (Hflux)
	deallocate (tempwave)
	deallocate (tempflux)
	deallocate (modwave)
	deallocate (modflux)
!
	return
	end

C-----------------------------------------------------------------------------
C
C subroutine determines the minimum delta lambda found in the model spectra
C

	subroutine get_model_res(wave,n,model_res)

	implicit none

	integer*4 n
	real*8 wave(n)
	real*8 model_res

	integer*4 i

	model_res = wave(n) - wave(1)
	do i=2,n
	   if ( (wave(i)-wave(i-1)) .lt. model_res) then
	      model_res = wave(i)-wave(i-1)
	    end if
	end do

	return
	end
	

C-----------------------------------------------------------------------------
C
C subroutine to make points linearly spaced in wavelength to a given
C delta lambda via linear interpolation
C

	subroutine linearize(new_wave,new_flux,n_new,
	1                    old_wave,old_flux,n_old,dlam)

	implicit none

	integer*4 n_new
	real*8 new_wave(n_new)
	real*8 New_flux(n_new)
!
	integer*4 n_old
	real*8 old_wave(n_old)
	real*8 old_flux(n_old)
!
	real*8 dlam 	! even spacing in wavelength

	integer*4 nnew
	integer*4 i
	real*8 dlamover2

	dlamover2 = dlam/2.0d0

	nnew = (old_wave(n_old)-old_wave(1))/dlam  ! this will shorten wavelength 
                                          ! array by a fraction of dlam

	if (nnew .gt. n_new) then
	   write(*,*)'Error in linearize, n_newX = ',n_new
	   write(*,*)'too small for this resolution'
	   stop 
	endif
	n_new=nnew	
!
	new_wave(1) = old_wave(1) + dlamover2
	do i=2,n_new
	   new_wave(i) = new_wave(i-1) + dlam
	end do

	call map2(old_wave,old_flux,n_old,
	1            new_wave,new_flux,n_new,dlam)

	return
	end

C-----------------------------------------------------------------------------
C subroutine to build wavelength array from frequency array (if forward=true)
C or build frequency array from wavelength array (if forward=false)
C 
C
	subroutine nu_to_lambda(freq,wave,n,nutolambda)

	implicit none
	include 'constants.inc'

	integer*4 n
	real*8 freq(n)
	real*8 wave(n)
	logical nutolambda ! true=nu->lambda, false=lambda->nu

	integer*4 i

	if (nutolambda) then
	 do i=1,n
	   if (freq(i) .eq. zero) then
	      write(*,*)'zero frequency found at i=',i
	   endif
	   wave(i) = (c_kms*10D-3)/freq(i) ! wavelength in Angstroms
	 end do
	else
	 do i=1,n
          if (wave(i) .eq. zero) then
	    write(*,*)'zero wavelength found at i=',i
	  endif
	  freq(i) = (c_kms*10D-3)/wave(i) ! wavelength in Angstroms
         end do
        end if

	return
	end

C-----------------------------------------------------------------------------
C subroutine to convert lambda array to ln(lambda) array if forward=true
C or convert ln(lambda) array to lambda array if forward=false
C
	subroutine lambda_to_lnlambda(wave,n,ltolnl)

	implicit none
	include 'constants.inc'

	integer*4 n
	real*8 wave(n)
	logical ltolnl  ! if T, lambda->lnlambda, if F, nlambda->lambda.

	integer*4 i

	if (ltolnl) then
	 do i=1,n
	  wave(i) = dlog(wave(i))
	 end do
	else 
	 do i=1,n
	  wave(i) = dexp(wave(i))
	 end do
	end if 

	return
	end

C-----------------------------------------------------------------------------
C subroutine to remove duplicate data points (poinst with same frequency
C
	subroutine remove_dups(freq,flux,n)

	implicit none
	include 'constants.inc'
	include 'parameters.inc'

	integer*4 n
	real*8 freq(n)
	real*8 flux(n)

	integer*4 i,j
	real*8 prev

	prev = zero
	j=0
	do i=1,n
	   if (prev .ne. freq(i)) then
	      j = j+1
	      freq(j) = freq(i)
	      flux(j) = flux(i)
	   else
	      write(*,*)'Duplicate data point ignored at nu = ',freq(i)
	   endif
	   prev = freq(i)
	end do
	n = j
	
	return
	end

C-----------------------------------------------------------------------------
C
C subroutine to initialize the H absorption array
C
	subroutine Habsinit(wave,flux,freq,modwave,n)

	implicit none
	include 'constants.inc'

	integer*4 n
	real*8 modwave(n)
	real*8 wave(n),flux(n),freq(n)

	integer*4 i

	wave(1:n) = modwave(1:n)
	freq(1:n) = (c_kms*10D-3)/wave(1:n)
	flux(1:n) = one

	return
	end

C-----------------------------------------------------------------------------
C subroutine to map array A onto array B using linear interpolation

	subroutine map(Awave,Aflux,AN,Bwave,Bflux,BN)

	implicit none
	include 'constants.inc'

	integer*4 AN,BN
	real*8 Awave(AN),Aflux(AN)
	real*8 Bwave(BN),Bflux(BN)

	integer*4 i,j

	i = 1
	j = 1

	do while (Bwave(j) .le. Awave(1))
	 j = j+1
	end do
	do while ((j .lt. BN) .and. (Bwave(j) .lt. Awave(AN)))
	 do while (.not.((Bwave(j).ge.Awave(i)).and.(Bwave(j).le.Awave(i+1))))
	  i = i+1
	 end do
	 Bflux(j) = Aflux(i) + (Aflux(i+1)-Aflux(i))*
     1   	  (Bwave(j) - Awave(i))/(Awave(i+1) - Awave(i))
	 j = j+1
	end do

	return 
	end


C-----------------------------------------------------------------------------
C subroutine to map array A onto array B using weighted averaging of fluxes.
C To determine the flux of a point in B, all points withing +/- dlam/2 are
C considered, as well as the the points longward and shortward of this
C range.  The fluxes are weighted by their wavelength distance (?) from
C the considered point.

	subroutine map2(Awave,Aflux,AN,Bwave,Bflux,BN,dlam)

	implicit none
	include 'constants.inc'

	integer*4 AN,BN
	real*8 Awave(AN),Aflux(AN)
	real*8 Bwave(BN),Bflux(BN)
	real*8 dlam,dlamover2

	integer*4 i,j,i_min,i_max
	integer*4 smaller_index,larger_index

	real*8 fluxsum
	real*8 dlamsum

	dlamover2=dlam/two
	i = 1
	j = 1

	do while (Bwave(j) .le. Awave(1))
	  j = j+1
	end do
!
	i_min=1
	i_max=1
	do while ((j .le. BN) .and. (Bwave(j) .lt. Awave(AN)))
C
C        when determining the flux of a point on the new grid (B), we
C        wish to consider all points within that points resolution range
C        (+/- dlam/2).
C        However, in the case where no old grid-points lie within this
C        range on one side or the other, we need to find the next nearest
C        point outside that range.  The new flux is the average
C        over these points, weighted by distance from the new point.
C        Warning: possible danger of coefficient blowing up if it is
C        very close to wavelength of point of interest (or same as.)
C
	 do while(awave(i_min) .lt. bwave(j)-dlamover2)
	   i_min=i_min+1
	 end do
	 if(awave(i_min) .gt. bwave(j))i_min=i_min-1
	 do while(awave(i_max) .lt. bwave(j)+dlamover2)
	   i_max=i_max+1
	 end do
	 if(awave(i_max) .lt. bwave(j))i_max=i_max+1
!
	 fluxsum = zero
	 dlamsum = zero    ! this is actulally the sum of dlam^-1
	 i = i_min
	 do while (i .le. i_max)
	   dlamsum = dlamsum + one/(abs(Bwave(j) - Awave(i)))
	   fluxsum = fluxsum + Aflux(i)/abs(Bwave(j) - Awave(i))
	   i = i+1
	 end do
	 Bflux(j) = fluxsum/dlamsum
	 j = j+1
	end do

	return 
	end

C-----------------------------------------------------------------------------
C function which returns the index of the first array element *larger* than
C or equal to the supplied number

	function larger_index(list,n,value)

	implicit none

	integer*4 larger_index,n
	real*8 list(n),value

	integer*4 i

	if (list(n) .lt. value) then
	   write(*,*)'Error - no element larger than or equal to ',value
	   ! some error handling function here
	endif
	i = 1
	do while (value .ge. list(i) .and. i .lt. n)
	   i = i+1
	end do
	larger_index = i

	return
	end

C-----------------------------------------------------------------------------
C function which returns the index of the last array element *smaller* or
C equal to the supplied number

	function smaller_index(list,n,value)

	implicit none

	integer*4 smaller_index,n
	real*8 list(n),value

	integer*4 i

	if (list(1) .gt. value) then
	   write(*,*)'Error - no element smaller than or equal to ',value
	   ! some error handling function here
	endif
	i = 1
	do while (list(i) .le. value .and. i .le. n)
	   i = i+1
	end do
	smaller_index = i-1

	return
	end

