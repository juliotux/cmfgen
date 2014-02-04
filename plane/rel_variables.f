c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine rel_variables(nd,chi,eta,nu_dnu,b,
     *       I_prev,chi_tau,source_prime)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Determine the needed variables for the relativistic formal solution.
c The method of Short Characteristic from Olson & Kunsasz (87) JQSRT
c and Hauschildt (92) JQSRT is used.
c
c The general form of the transfer equation is:
c
c   d(I)/d(tau) + I = S
c
c For relativistic tranfer:
c
c   d(I)/(chi_tau*ds) + I = source_prime
c
c The variables calculated and returned are:
c
c   source_prime=(eta+alpha*Iprev)/(alpha+chi_prime)
c   chi_tau=alpha+chi_prime
c
c   where:
c     chi_prime = chi+3*b_*
c     alpha = nu_dnu*b_* = freq(k)*/(freq(k-1)-freq(k))*b_*
c     b_* = gamma*((1-mu_*^2)*beta/r+gamma^2*mu_**(mu_*+beta)dbetadr)
c
c         where: _* is _p or _m for ray in plus or minus direction
c
c The variables passed to this routine are:
c      chi(nd) = total opacity (continuum, line, scattering)
c      eta(nd) = total emissivity
c        b(nd) = advection and abberation terms
c   I_prev(nd) = Intensity at previous frequency
c       nu_dnu = freq(k)*/(freq(k-1)-freq(k))
c
c written  5/23/97 DLM
c
c--------------------------------------------------------------------
c
      implicit none
c
c Grid size variables
c
      integer :: nd
c
c Opacity and emissivity variables
c
      real*8, dimension(nd) :: chi,eta
c
c Intensity variable
c
      real*8, dimension(nd) :: I_prev
c
c Frequency variable
c
      real*8 nu_dnu
c
c Advection and abberation terms
c
      real*8, dimension(nd) :: b
c
c Transfer variables
c
      real*8, dimension(nd) :: chi_tau
      real*8, dimension(nd) :: source_prime
c
c--------------------------------------------------------------------
c
c Use Fortran 90 array operations
c
c            chi_tau = alpha+chi_prime
c
      chi_tau=(nu_dnu+3.0d0)*b+chi
c
c            source_prime = (eta+alpha*I_prev)/(alpha+chi_prime)
c
      source_prime=(eta+nu_dnu*b*I_prev)/chi_tau
c
      return
      end
c
