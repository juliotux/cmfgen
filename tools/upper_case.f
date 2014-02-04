c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      function upper_case(string) result(all_caps)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Function to change all lower case letters in the variable
c string to upper case.
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      implicit none
c
      integer :: i,j,nchar
      character(len=*) :: string
      character(len=120) :: all_caps
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Change all letters to uppercase
c
      all_caps="  "
      nchar=len_trim(string)
      upper: do i=1,nchar
        j=ichar(string(i:i))
        select case (j)
          case (97:122)                   !a-z
            all_caps(i:i)=char(j-32)      !A-Z
          case default
            all_caps(i:i)=char(j)
        end select
      enddo upper
c
      return
      end
c
