c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine runge_kutta(y,dydx,n,x,h,beta,dbetadr,gamma,deriv)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Numerical Recipies routine to calculate the 4th order Runge-Kutta
c solution of a set of ordinary differential equations
c
c Altered 2/21/96 DLM Updated to F90 standard
c
c--------------------------------------------------------------------
c
      implicit none
c
      external deriv
c
c Number of equations to solve: passed from calling routine
c
      integer :: n
c
c Local variables
c
      integer :: i
      integer, parameter :: nmax=2
      real*8, dimension(nmax) :: dym,dyt,yt
      real*8 h6,hh,xh
c
c  Value to be determine
c
      real*8,  dimension(n) :: dydx,y
      real*8 h,x,beta,dbetadr,gamma
c
      hh=h*0.5d0
      h6=h/6.0d0
      xh=x+hh
      x=x+h
      do i=1,n
        yt(i)=y(i)+hh*dydx(i)
      enddo
      call deriv(xh,yt,dyt,beta,dbetadr,gamma)
      do i=1,n
        yt(i)=y(i)+hh*dyt(i)
      enddo
      call deriv(xh,yt,dym,beta,dbetadr,gamma)
      do i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
      enddo
      call deriv(x,yt,dyt,beta,dbetadr,gamma)
      do i=1,n
        y(i)=y(i)+h6*(dydx(i)+dyt(i)+2.0d0*dym(i))
      enddo
c
      return
      end
c
