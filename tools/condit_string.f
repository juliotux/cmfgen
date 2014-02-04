c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      function condit_string(answer,nchar) result(string)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Subroutine to condition an input string.  On entry STRING is
c parameter string input by the user and NCHAR is 0.  On exit,
c STRING will have been conditioned to the format required by the
c calling routine and NCHAR will be the number of characters in
c STRING that are used.
c
c Conditioning consists of removing extra spaces (one space between
c each parameter), converting all ","s to spaces, and making all
c ASCII letters uppercase.
c
c written  11/19/96  DLM  Modeled after FUNCTION CONDIT(A,N) that
c                         worked only under VMS.
c
c altered   2/27/96  DLM  Changed name to condit_string
c
c                         Changed to a function instead of a subroutine
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      implicit none
c
      integer :: nchar,i,j
c
      character(len=*) :: answer
      character(len=120) :: string
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      string=adjustl(answer)
c     
      nchar=len_trim(string)
c
c Remove all double spaces between parameters
c
      i=index(string,'  ')
      spaces: do
        i=index(string,'  ')
        if(i.gt.nchar)exit
        string(i:)=string(i+1:)//' '
        nchar=nchar-1
      end do spaces
c
c Convert all ","s to spaces
c
      commas: do
        i=index(string,',')
        if(i.eq.0)exit
        string(i:i)=' '
      enddo commas
c
c Change all letters to uppercase
c
      upper: do i=1,nchar
        j=ichar(string(i:i))
        select case (j)
          case (97:122)                 !a-z
            string(i:i)=char(j-32)      !A-Z
        end select
      enddo upper
c
      return
      end
c
