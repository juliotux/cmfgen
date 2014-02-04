c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      module mod_wr_string
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      interface wr_string
        module procedure wr_real,
     *                   wr_real_dim,
     *                   wr_integer,
     *                   wr_integer_dim,
     *                   wr_logical,
     *                   wr_character
      end interface
c
      contains
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      function wr_real(input) result(default)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c This function write the real variable "input" to the string
c "default".  It automatically decides which format type is best
c (F or E) and trims zeros in decimal places
c
c written  11/24/96  DLM
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      implicit none
c
      integer :: num_char,i
c
      integer, parameter :: num_e_decimal=6
      integer, parameter :: num_f_decimal=6
c
      real*8 :: input,t1
c
      character(len=120) :: default
      character(len=20) :: fmt
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Compute number of digits to be written
c
      default=' '
      if((abs(input).ge.1.0d+06).or.(abs(input).le.1.0d-3))then
c
c E format
c
        num_char=6+num_e_decimal
        if(input.lt.0.0)num_char=num_char+1       !minus sign
        write(fmt,'(a5,i2.2,a1,i2.2,a1)')
     *       '(1p,e',num_char,'.',num_e_decimal,')'
        write(default,fmt)input
c
c Remove zeros in decimal place
c
        i=index(default,'E')-1
        do while (default(i:i).eq.'0')
          if(default(i-1:i-1).ne.'.')default(i:)=default(i+1:)
          i=i-1
        end do
c
      else
c
c F format
c
        t1=abs(input)
        num_char=2
        if(t1.gt.1.0)then
          num_char=num_char+log10(t1)
        end if
        num_char=num_char+num_f_decimal
        if(input.lt.0.0)num_char=num_char+1       !minus sign
c
        if(num_char.lt.10)then
          write(fmt,"('(f',i1,'.',i1,')')")num_char,num_f_decimal
        else
          write(fmt,"('(f',i2,'.',i1,')')")num_char,num_f_decimal
        endif
        write(default,fmt)input
c
c Remove zeros in decimal place
c
        i=len_trim(default)
        do while (default(i:i).eq.'0')
          if(default(i-1:i-1).ne.'.')default(i:i)=' '
          i=i-1
        end do
c
      endif      
c
      return
      end function wr_real
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      function wr_real_dim(input,dim) result(default)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c This function write the real variable "input" to the string
c "default".  It automatically decides which format type is best
c (F or E) and trims zeros in decimal places
c
c written  11/24/96  DLM
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      implicit none
c
      integer :: num_char,i,j,dim,position
c
      integer, parameter :: num_e_decimal=6
      integer, parameter :: num_f_decimal=6
c
      real*8, dimension(dim) :: input
      real*8 :: t1
c
      character(len=120) :: default
      character(len=20) :: fmt
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Compute number of digits to be written
c
      default=' '
      position=1
c
c Loop over dimension of array "input"
c
      do j=1,dim
c
        if((abs(input(j)).ge.1.0d+06).or.(abs(input(j)).le.1.0d-3))then
c
c E format
c
          num_char=6+num_e_decimal
          if(input(j).lt.0)num_char=num_char+1     !minus sign
          write(fmt,'(a5,i2.2,a1,i2.2,a1)')
     *         '(1p,e',num_char,'.',num_e_decimal,')'
          write(default(position:),fmt)input(j)
c
c Remove zeros in decimal place
c
          i=position+index(default(position:),'E')-2
          do while (default(i:i).eq.'0')
            if(default(i-1:i-1).ne.'.')default(i:)=default(i+1:)
            i=i-1
          end do
c
        else
c
c F format
c
          t1=abs(input(j))
          num_char=2
          if(t1.gt.1.0)then
            num_char=num_char+log10(t1)
          end if
          num_char=num_char+num_f_decimal
          if(input(j).lt.0)num_char=num_char+1 !minus sign
c
          if(num_char.lt.10)then
            write(fmt,"('(f',i1,'.',i1,')')")num_char,num_f_decimal
          else
            write(fmt,"('(f',i2,'.',i1,')')")num_char,num_f_decimal
          endif
          write(default(position:),fmt)input(j)
c
c Remove zeros in decimal place
c
          i=len_trim(default)
          do while (default(i:i).eq.'0')
            if(default(i-1:i-1).ne.'.')default(i:i)=' '
            i=i-1
          end do
c
        endif
c
        if(j.lt.dim)then
          position=len_trim(default)+1
          default(position:position)=","
          position=position+1
        endif
c
      enddo
c
      return
      end function wr_real_dim
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      function wr_integer(input) result(default)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c This function writes the integer variable "input" to the string
c "default".
c
c written  11/24/96  DLM
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      implicit none
c
      integer :: num_char
c
      integer :: input
c
      real*8 :: t1
c
      character(len=120) :: default
      character(len=20) :: fmt
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Compute number of digits to be written
c
      default=' '
c
      t1=abs(float(input))
      num_char=1
      if(t1.gt.1.0)then
        num_char=num_char+log10(t1)
      end if
      if(input.lt.0)num_char=num_char+1    !minus sign
c
      if(num_char.lt.10)then
        write(fmt,"('(i',i1,')')")num_char
      else
        write(fmt,"('(i',i2,')')")num_char
      endif
      write(default,fmt)input
c
      return
      end function wr_integer
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      function wr_integer_dim(input,dim) result(default)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c This function write the integer variable "input" to the string
c "default".
c
c written  11/24/96  DLM
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      implicit none
c
      integer :: num_char,j,dim,position
c
      integer, dimension(dim) :: input
c
      real*8 :: t1
c
      character(len=120) :: default
      character(len=20) :: fmt
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Compute number of digits to be written
c
      default=' '
      position=1
c
c Loop over dimension of array "input"
c
      do j=1,dim
c
        t1=abs(float(input(j)))
        num_char=1
        if(t1.gt.1.0)then
          num_char=num_char+log10(t1)
        end if
        num_char=num_char
        if(input(j).lt.0)num_char=num_char+1     !minus sign
c
        if(num_char.lt.10)then
          write(fmt,"('(i',i1,')')")num_char
        else
          write(fmt,"('(i',i2,')')")num_char
        endif
        write(default(position:),fmt)input(j)
c
        if(j.lt.dim)then
          position=len_trim(default)+1
          default(position:position)=","
          position=position+1
        endif
c
      enddo
c
      return
      end function wr_integer_dim
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      function wr_logical(input) result(default)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c This function writes the logical variable "input" to the string
c "default".
c
c written  11/24/96  DLM
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      implicit none
c
      logical :: input
c
      character(len=120) :: default
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      if(input)then
        default='T'
      else
        default='F'
      endif
c
      return
      end function wr_logical
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      function wr_character(input) result(default)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c This function writes the string variable "input" to the string
c "default".  This routine is included only for completeness.
c It would easier to set "default=input" in the calling routine.
c
c written  11/24/96  DLM
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      implicit none
c
      character(len=*) :: input
c
      character(len=120) :: default
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      default=input
c
      return
      end function wr_character
c
      end module mod_wr_string
c
