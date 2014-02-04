c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      function wr_fmt_string(length) result(fmt_string)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Function to create a format string to needed to read a string
c from an internal string.  SGI would do this with the default "(a)"
c format, but VMS would give access violation if string being
c read was longer than string being read from.
c
c  create 1/13/96  DLM
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      implicit none
c
      integer :: length
c
      character(len=120) :: fmt_string
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      if(length.lt.10)then
        write(fmt_string,"('(a',i1,')')")length
      elseif(length.lt.100)then
        write(fmt_string,"('(a',i2,')')")length
      else
        write(fmt_string,"('(a',i3,')')")length
      endif
c
      return
      end
c
