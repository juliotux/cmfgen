c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      function wr_pwd()
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Function to return the present working directory (pwd) on a
c UNIX machine.  The pwd string is truncated after the user's
c name, if it is contained in the pwd.
c
c On VMS systems the string will return blank.
c
c  created  1/13/96  DLM
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      implicit none
c
      integer i
c
      character(len=120) wr_pwd,pwd,user
c
       call getenv('PWD',pwd)
       call getenv('USER',user)
       i=index(pwd,trim(user))
       if(i.ne.0)then
         wr_pwd=pwd(i+9:)
       else
         wr_pwd=pwd
       endif
!
      wr_pwd=' '
c
      return
c
      stop
      end
c
