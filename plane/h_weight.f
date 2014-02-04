c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine h_weight(nw,angle,temp)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Calculate integration weights for the flux, H.
c (see description in routine weights and blue derivation book)
c
c Written    8-95 DLM
c
c Altered    8-95 DLM Made more general.  Checks now depend on end points
c                     instead of 0->1.
c
c Altered 2/21/96 DLM Updated to F90 standard
c
c--------------------------------------------------------------------
c
      implicit none
c
c Number of angles: passed from calling routine
c
      integer :: nw
c
c Local and temparary variables
c
      integer :: j
c
      real*8 :: a1,a2,sum,check
      real*8, dimension(nw) :: temp
c
c Angle variables
c
      real*8, dimension(nw) :: angle
c
c--------------------------------------------------------------------
c
      temp(:)=0.0d0
c
      do j=1,nw-1
        a1=(angle(j)+angle(j+1))/4.0d0
        a2=(angle(j)*angle(j)+
     *      angle(j)*angle(j+1)+
     *      angle(j+1)*angle(j+1))/6.0d0
        temp(j)=temp(j)-angle(j+1)*a1+a2
        temp(j+1)=temp(j+1)+angle(j)*a1-a2
      enddo
c
c Check that 1/2*int[mu*dmu]=[mu_max^2-mu+min^2]/4
c
      sum=0.0d0
      do j=1,nw
        sum=sum+temp(j)
      enddo
c
      check=(angle(1)**2-angle(nw)**2)/4.0d0
c
      if(abs(check-sum).gt.1.0d-12)then
        print*,' H integration weights need to be normalized (1)'
        print*,' error=',abs(check-sum)
      endif
c
c Check that 1/2*int[mu^2*dmu]=[mu_max^3-mu_min^3]/6
c
      sum=0.0d0
      do j=1,nw
        sum=sum+temp(j)*angle(j)
      enddo
c
      check=(angle(1)**3-angle(nw)**3)/6.0d0
c
      if(abs(check-sum).gt.1.0d-12)then
        print*,' H integration weights need to be normalized (2)'
        print*,' error=',abs(check-sum)
      endif
c
      return
      end
