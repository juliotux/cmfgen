c
c    **** MODULE USR_HIDDEN ***
c
c This module creates the generic subroutine called USR_HIDDEN.
c When a call is made to USR_HIDDEN, the interface determines
c which subroutine to call by matching the types and sizes of
c variables passed in the call.
c
c
c The possible subroutines are:
c
c   RD_HIDDEN_REAL           => read one hidden real variable
c   RD_HIDDEN_REAL_DIM       => read an array of hidden real variable
c   RD_HIDDEN_INTEGER        => read one hidden integer variable
c   RD_HIDDEN_INTEGER_DIM    => read an array of hidden integer variables
c   RD_HIDDEN_LOGICAL        => read one hidden logical variable
c   RD_HIDDEN_STRING         => read one hidden string variable
c
c The general call to USR_HIDDEN is for one input variable:
c
c   CALL USR_HIDDEN(VARIABLE,'VAR_NAME','DEFAULT','DESCRIPTION')
c
c and for an array of real or integer variables:
c
c   CALL USR_HIDDEN(VARIABLE,DIMENSION,REQUIRED,'VAR_NAME',
c  *    'DEFAULT','DESCRIPTION')
c
c where:
c
c     VARIABLE             => variable or array to be read
c
c     DIMENSION (optional) => dimension of array
c
c     REQUIRED  (optional) => required number of array variables to read
c
c     'VAR_NAME'           => string or string variable used to flag an
c                             option in .sve file
c     'DEFAULT'            => string or string variable containing defaults
c                             values for variable
c     'DESCRIPTION'        => string or string variable containing a description
c                             of the variable
c
c    **** MODULE USR_HIDDEN ***
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c written  3/4/97  DLM  Modeled after usr_option
c
c altered  3/17/97 DLM  Removed main_option from list of passed 
c                         variables.  Now store hidden options in
c                         variable string.
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      module mod_usr_hidden
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      interface usr_hidden
        module procedure rd_hidden_real,
     *                   rd_hidden_real_dim,
     *                   rd_hidden_integer,
     *                   rd_hidden_integer_dim,
     *                   rd_hidden_logical,
     *                   rd_hidden_string
      end interface
c
      contains
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine rd_hidden_real(input,var_name,default,description)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c This routine reads a real value into the hidden variable
c input. The user is never prompted for a hidden variable and the
c default value passed to this subroutine is always used unless the
c hidden option is entered with the main option. For example:
c   rd_obs (scale=1.2)
c The value for the hidden sub-option "scale" is set to 1.2.  This value
c for the hidden option (in this case scale) is stored in the variable
c string by the routine pre_ and post_sve_file.  Hidden
c options are not written in .sve files unless their values are changed
c
c written   3/4/96  DLM  Modified copy of usr_option
c
c altered  3/17/97  DLM  Removed main_option from list of passed 
c                          variables.  Now store hidden options in
c                          variable string.
c                        Now read and write var_name as all capital
c                          letters.  Cannot reassign var_name because
c                          a specific string is passed to usr_option
c                          instead of a character variable.  Thus must
c                          use another variable, all_caps, to create
c                          an all upper case string (SGI will allow
c                          change to var_name, but VMS will not!).
c                        Use "{}" instead of "[]" to delimit sub-options
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      use mod_sve_file
      implicit none
c
      integer :: l,o_start,o_end
c
      real*8 :: input
c
      character(len=*) :: var_name,description,default
      character(len=120) :: answer,sve_string,option,all_caps
c
      logical :: present
c
      external upper_case
      character(len=120) :: upper_case
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Write and search for variable names in all capital letters
c
      all_caps=upper_case(var_name)
      answer=" "
c
c Determine if .sve file is being read.  If so, read hidden option
c into answer.
c
      call pre_sve_file(answer,all_caps,present)
c
c Check string for hidden option.
c
      option="{"//trim(all_caps)//"="
      l=index(string,trim(option))
      if(l.gt.0)then
c
c If hidden option found then put into answer
c
        o_start=l+len_trim(option)
        o_end=o_start+index(string(o_start:),'}')-2
        answer=string(o_start:o_end)
c        if(.not.sve_read)
c     *       write(*,"(8x,a,a,'   **hidden option changed**')")
c     *       trim(option(2:)),trim(answer)
c
      endif
c
c Read input either from answer or default
c
      if(len_trim(answer).ne.0)then
        read(answer,*,err=300)input
        sve_string=answer
c
      elseif(len_trim(default).ne.0)then
        read(default,*)input
        sve_string=default
      endif
c
c Write option to current .sve file., and ouput to terminal.
C
      if(.not. sve_read)then
        write(*,"(1X,A,' =>',A)")trim(all_caps),trim(sve_string)
      end if
      call wr_sve_file(sve_string,all_caps)
c
      return
c
c If given value for hidden option cannot be read, print message
c and read default value.
c
 300  print*,' error in rd_hidden_real'
      print*,' not able to read hidden value'
      print*,trim(option),trim(answer)
      print*
      print*,' default value used'
      read(default,*)input
      sve_string=default
c
      return
      end subroutine rd_hidden_real
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine rd_hidden_real_dim(input,dim,required,var_name,
     *     default,description)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c This routine reads an array of real values into the hidden variable
c input. The user is never prompted for a hidden variable and the
c default value passed to this subroutine is always used unless the
c hidden option is entered with the main option. For example:
c   rd_obs (scales=1.1,1.2,1.3)
c The value for the hidden sub-option "scales" is set to 1.1,1.2,1.3.
c These value for the hidden option (in this case scale) are stored in
c the variable string by the routine pre_ and post_sve_file. Hidden
c options are not written in .sve files unless their values are changed
c
c written   3/4/96  DLM  Modified copy of usr_option
c
c altered  3/17/97  DLM  Removed main_option from list of passed 
c                          variables.  Now store hidden options in
c                          variable string.
c                        Now read and write var_name as all capital
c                          letters.  Cannot reassign var_name because
c                          a specific string is passed to usr_option
c                          instead of a character variable.  Thus must
c                          use another variable, all_caps, to create
c                          an all upper case string (SGI will allow
c                          change to var_name, but VMS will not!).
c                        Use "{}" instead of "[]" to delimit sub-options
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      use mod_wr_string
      use mod_sve_file
      implicit none
c
      integer :: dim,required
      integer :: i,l,o_start,o_end,found
      integer :: m,start,stop,step,count
c
      real*8, dimension(dim) :: input
      real*8 t1
c
      character(len=*) :: var_name,description,default
      character(len=120) :: answer,sve_string,option,number,all_caps
c
      logical :: present
c
      external upper_case
      character(len=120) :: upper_case
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Write and search for variable names in all capital letters
c
      all_caps=upper_case(var_name)
      answer=" "
c
c Determine if .sve file is being read.  If so, read hidden option
c into answer.
c
      call pre_sve_file(answer,all_caps,present)
c
c Check main option for hidden option.
c
      option="{"//trim(all_caps)//"="
      l=index(string,trim(option))
      if(l.gt.0)then
c
c If hidden option found then put into answer
c
        o_start=l+len_trim(option)
        o_end=o_start+index(string(o_start:),'}')-2
        answer=string(o_start:o_end)
c        if(.not.sve_read)
c     *       write(*,"(8x,a,a,'   **hidden option changed**')")
c     *       trim(option(2:)),trim(answer)
c
      endif
c
c Read input() either from answer or default
c
      l=len_trim(answer)
      if(l.ne.0)then
        l=index(answer,':')
        if(l.gt.0)then
          read(answer(:l-1),*)start
          m=index(answer(l+1:),':')
          if(m.gt.0)then
            read(answer(l+1:l+m-1),*)stop
            read(answer(l+m+1:),*)step
          else
            read(answer(l+1:),*)stop
            step=1
          endif
          count=0
          do l=start,stop,step
            t1=float(l)
            number=wr_string(t1)
            if(l.eq.start)then
              answer=number
            else
              answer=trim(answer)//","//trim(number)
            endif
            count=count+1
            if(count.eq.dim)exit
          enddo
          write(*,"(8x,a,'=',a)")trim(var_name),trim(answer)
        endif
        l=len_trim(answer)
        input(:)=0.0
        read(answer(:l),*,end=100,err=300)(input(i),i=1,dim)
 100    do i=dim,1,-1
          if(input(i).ne.0.0)exit
        enddo
        found=i
        if(found.lt.required)goto 300
        sve_string=answer
c
      elseif(len_trim(default).ne.0)then
        read(default,*,end=400)(input(i),i=1,dim)
 400    sve_string=default
      endif
c
c Write option to current .sve file., and ouput to terminal.
C
      if(.not. sve_read)then
        write(*,"(1X,A,' =>',A)")trim(all_caps),trim(sve_string)
      end if
      call wr_sve_file(sve_string,all_caps)
c
      return
c
c If given value for hidden option cannot be read, print message
c and read default value.
c
 300  print*,' error in rd_hidden_real_dim'
      write(*,"(1x,i2,' values required, only found ',i2)")
     *     required,found
      print*,' default value used'
      read(answer,*,end=500)(input(i),i=1,dim)
 500  sve_string=default
c
      return
      end subroutine rd_hidden_real_dim
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine rd_hidden_integer(input,var_name,default,description)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c This routine reads an integer value into the hidden variable
c input. The user is never prompted for a hidden variable and the
c default value passed to this subroutine is always used unless the
c hidden option is entered with the main option. For example:
c   rd_obs (level=40)
c The value for the hidden sub-option "level" is set to 40. This value
c for the hidden option (in this case level) is stored in the variable
c string by the routine pre_ and post_sve_file.  Hidden
c options are not written in .sve files unless their values are changed
c
c written   3/4/96  DLM  Modified copy of usr_option
c
c altered  3/17/97  DLM  Removed main_option from list of passed 
c                          variables.  Now store hidden options in
c                          variable string.
c                        Now read and write var_name as all capital
c                          letters.  Cannot reassign var_name because
c                          a specific string is passed to usr_option
c                          instead of a character variable.  Thus must
c                          use another variable, all_caps, to create
c                          an all upper case string (SGI will allow
c                          change to var_name, but VMS will not!).
c                        Use "{}" instead of "[]" to delimit sub-options
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      use mod_sve_file
      implicit none
c
      integer :: l,o_start,o_end
c
      integer :: input
c
      character(len=*) :: var_name,description,default
      character(len=120) :: answer,sve_string,option,all_caps
c
      logical :: present
c
      external upper_case
      character(len=120) :: upper_case
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Write and search for variable names in all capital letters
c
      all_caps=upper_case(var_name)
      answer=" "
c
c Determine if .sve file is being read.  If so, read hidden option
c into answer.
c
      call pre_sve_file(answer,all_caps,present)
c
c Check main option for hidden option.
c
      option="{"//var_name//"="
      l=index(string,trim(option))
      if(l.gt.0)then
c
c If hidden option found then put into answer
c
        o_start=l+len_trim(option)
        o_end=o_start+index(string(o_start:),'}')-2
        answer=string(o_start:o_end)
c        if(.not.sve_read)
c     *       write(*,"(8x,a,a,'   **hidden option changed**')")
c     *       trim(option(2:)),trim(answer)
c
      endif
c
c Read input either from answer or default
c
      if(len_trim(answer).ne.0)then
        read(answer,*,err=300)input
        sve_string=answer
c
      elseif(len_trim(default).ne.0)then
        read(default,*)input
        sve_string=default
      endif
c
c Write option to current .sve file., and ouput to terminal.
C
      if(.not. sve_read)then
        write(*,"(1X,A,' =>',A)")trim(all_caps),trim(sve_string)
      end if
      call wr_sve_file(sve_string,all_caps)
C
      return
c
c If given value for hidden option cannot be read, print message
c and read default value.
c
 300  print*,' error in rd_hidden_integer'
      print*,' not able to read hidden value'
      print*,trim(option),trim(answer)
      print*
      print*,' default value used'
      read(default,*)input
      sve_string=default
c
      return
      end subroutine rd_hidden_integer
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine rd_hidden_integer_dim(input,dim,required,var_name,
     *     default,description)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c This routine reads an array of integer values into the hidden variable
c input. The user is never prompted for a hidden variable and the
c default value passed to this subroutine is always used unless the
c hidden option is entered with the main option. For example:
c   rd_obs (levels=1,2,3)
c The value for the hidden sub-option "levels" is set to 1,2,3.  Thes values
c for the hidden option (in this case levels) are stored in the variable
c string by the routine pre_ and post_sve_file.  Hidden
c options are not written in .sve files unless their values are changed
c
c written   3/4/96  DLM  Modified copy of usr_option
c
c altered  3/17/97  DLM  Removed main_option from list of passed 
c                          variables.  Now store hidden options in
c                          variable string.
c                        Now read and write var_name as all capital
c                          letters.  Cannot reassign var_name because
c                          a specific string is passed to usr_option
c                          instead of a character variable.  Thus must
c                          use another variable, all_caps, to create
c                          an all upper case string (SGI will allow
c                          change to var_name, but VMS will not!).
c                        Use "{}" instead of "[]" to delimit sub-options
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      use mod_wr_string
      use mod_sve_file
      implicit none
c
      integer :: dim,required
      integer :: i,l,o_start,o_end,found
      integer :: m,start,stop,step,count
c
      integer, dimension(dim) :: input
c
      character(len=*) :: var_name,description,default
      character(len=120) :: answer,sve_string,option,number,all_caps
c
      logical :: present
c
      external upper_case
      character(len=120) :: upper_case
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Write and search for variable names in all capital letters
c
      all_caps=upper_case(var_name)
      answer=" "
c
c Determine if .sve file is being read.  If so, read hidden option
c into answer.
c
      call pre_sve_file(answer,all_caps,present)
c
c Check main option for hidden option.
c
      option="{"//trim(all_caps)//"="
      l=index(string,trim(option))
      if(l.gt.0)then
c
c If hidden option found then put into answer
c
        o_start=l+len_trim(option)
        o_end=o_start+index(string(o_start:),'}')-2
        answer=string(o_start:o_end)
c        if(.not.sve_read)
c     *       write(*,"(8x,a,a,'   **hidden option changed**')")
c     *       trim(option(2:)),trim(answer)
c
      endif
c
c Read input either from answer or default
c
      l=len_trim(answer)
      if(l.ne.0)then
        l=index(answer,':')
        if(l.gt.0)then
          read(answer(:l-1),*)start
          m=index(answer(l+1:),':')
          if(m.gt.0)then
            read(answer(l+1:l+m-1),*)stop
            read(answer(l+m+1:),*)step
          else
            read(answer(l+1:),*)stop
            step=1
          endif
          count=0
          do l=start,stop,step
            number=wr_string(l)
            if(l.eq.start)then
              answer=number
            else
              answer=trim(answer)//","//trim(number)
            endif
            count=count+1
            if(count.eq.dim)exit
          enddo
          write(*,"(8x,a,'=',a)")trim(var_name),trim(answer)
        endif
        l=len_trim(answer)
        input(:)=huge(l)
        read(answer(:l),*,end=100,err=300)(input(i),i=1,dim)
        found=0
 100    do i=dim,1,-1
          if(input(i) .ne. huge(l))then
            found=i
	    exit
	  else
	    input(i)=0
	  end if
        enddo
        if(found.lt.required)goto 300
        sve_string=answer
c
      elseif(len_trim(default).ne.0)then
        read(default,*,end=400)(input(i),i=1,dim)
 400    sve_string=default
      endif
c
c Write option to current .sve file., and ouput to terminal.
C
      if(.not. sve_read)then
        write(*,"(1X,A,' =>',A)")trim(all_caps),trim(sve_string)
      end if
      call wr_sve_file(sve_string,all_caps)
c
      return
c
c If given value for hidden option cannot be read, print message
c and read default value.
c
 300  print*,' error in rd_hidden_integer_dim'
      write(*,"(1x,i2,' values required, only found ',i2)")
     *     required,found
      print*,' default value used'
      read(default,*,end=500)(input(i),i=1,dim)
 500  sve_string=default
c
      return
      end subroutine rd_hidden_integer_dim
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine rd_hidden_logical(input,var_name,default,description)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c This routine reads a logical value into the hidden variable
c input. The user is never prompted for a hidden variable and the
c default value passed to this subroutine is always used unless the
c hidden option is entered with the main option. For example:
c   rd_obs (over=T)
c The value for the hidden sub-option "over" is set to TRUE.  This value
c for the hidden option (in this case over) is stored in the variable
c string by the routine pre_ and post_sve_file.  Hidden
c options are not written in .sve files unless their values are changed
c
c  written  3/4/96  DLM  Modified copy of usr_option
c
c altered  3/17/97  DLM  Removed main_option from list of passed 
c                          variables.  Now store hidden options in
c                          variable string.
c                        Now read and write var_name as all capital
c                          letters.  Cannot reassign var_name because
c                          a specific string is passed to usr_option
c                          instead of a character variable.  Thus must
c                          use another variable, all_caps, to create
c                          an all upper case string (SGI will allow
c                          change to var_name, but VMS will not!).
c                        Use "{}" instead of "[]" to delimit sub-options
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      use mod_sve_file
      implicit none
c
      integer :: l,o_start,o_end
c
      logical :: input
      character(len=*) :: var_name,description,default
      character(len=120) :: answer,sve_string,option,all_caps
c
      logical :: present
c
      external upper_case
      character(len=120) :: upper_case
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Write and search for variable names in all capital letters
c
      all_caps=upper_case(var_name)
      answer=" "
c
c Determine if .sve file is being read.  If so, read hidden option
c into answer.
c
      call pre_sve_file(answer,all_caps,present)
c
c Check string for hidden option.
c
      option="{"//trim(all_caps)//"="
      l=index(string,trim(option))
      if(l.gt.0)then
c
c If hidden option found then put into answer
c
        o_start=l+len_trim(option)
        o_end=o_start+index(string(o_start:),'}')-2
        answer=string(o_start:o_end)
c        if(.not.sve_read)
c     *       write(*,"(8x,a,a,'   **hidden option changed**')")
c     *       trim(option(2:)),trim(answer)
c
      endif
c
c Read input either from answer or default
c
      if(len_trim(answer).ne.0)then
        read(answer,*,err=300)input
        sve_string=answer
c
      elseif(len_trim(default).ne.0)then
        read(default,*)input
        sve_string=default
      endif
c
c Write option to current .sve file., and ouput to terminal.
C
      if(.not. sve_read)then
        write(*,"(1X,A,' =>',A)")trim(all_caps),trim(sve_string)
      end if
      call wr_sve_file(sve_string,all_caps)
c
      return
c
c If given value for hidden option cannot be read, print message
c and read default value.
c
 300  print*,' error in rd_hidden_logical'
      print*,' not able to read hidden value'
      print*,trim(option),trim(answer)
      print*
      print*,' default value used'
      read(default,*)input
      sve_string=default
c
      return
      end subroutine rd_hidden_logical
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine rd_hidden_string(input,var_name,default,description)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c This routine reads a string value into the hidden variable
c input. The user is never prompted for a hidden variable and the
c default value passed to this subroutine is always used unless the
c hidden option is entered with the main option. For example:
c   rd_obs (name=junk)
c The value for the hidden sub-option "name" is set to "junk".  This value
c for the hidden option (in this case junk) is stored in the variable
c string by the routine pre_ and post_sve_file.  Hidden
c options are not written in .sve files unless their values are changed
c
c written  3/4/96  DLM  Modified copy of usr_option
c
c altered  3/17/97  DLM  Removed main_option from list of passed 
c                          variables.  Now store hidden options in
c                          variable string.
c                        Now read and write var_name as all capital
c                          letters.  Cannot reassign var_name because
c                          a specific string is passed to usr_option
c                          instead of a character variable.  Thus must
c                          use another variable, all_caps, to create
c                          an all upper case string (SGI will allow
c                          change to var_name, but VMS will not!).
c                        Use "{}" instead of "[]" to delimit sub-options
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      use mod_sve_file
      implicit none
c
      integer :: l,o_start,o_end
c
      character(len=*) :: input
      character(len=*) :: var_name,description,default
      character(len=120) :: answer,sve_string,option
      character(len=120) :: fmt_string,all_caps
c
      logical :: present
c
      external wr_fmt_string
      character(len=120) :: wr_fmt_string
c
      external upper_case
      character(len=120) :: upper_case
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Write and search for variable names in all capital letters
c
      all_caps=upper_case(var_name)
      answer=" "
c
c Determine if .sve file is being read.  If so, read hidden option
c into answer.
c
      call pre_sve_file(answer,all_caps,present)
c
c Check main option for hidden option.
c
      option="{"//trim(all_caps)//"="
      l=index(string,trim(option))
      if(l.gt.0)then
c
c If hidden option found then put into answer
c
        o_start=l+len_trim(option)
        o_end=o_start+index(string(o_start:),'}')-2
        answer=string(o_start:o_end)
c        if(.not.sve_read)
c     *       write(*,"(8x,a,a,'   **hidden option changed**')")
c     *       trim(option(2:)),trim(answer)
c
      endif
c
c Read input either from answer or default
c
      if(len_trim(answer).ne.0)then
!        fmt_string=wr_fmt_string(l)
!        read(answer,fmt_string,err=300)input
	input=trim(answer)
        sve_string=answer
c
      elseif(len_trim(default).ne.0)then
	input=trim(default)
!       read(default,*)input
        sve_string=default
      endif
c
c Write option to current .sve file., and ouput to terminal.
C
      if(.not. sve_read)then
        write(*,"(1X,A,' =>',A)")trim(all_caps),trim(sve_string)
      end if
      call wr_sve_file(sve_string,all_caps)
c
      return
c
c If given value for hidden option cannot be read, print message
c and read default value.
c
 300  print*,' error in rd_hidden_string'
      print*,' not able to read hidden value'
      print*,trim(option),trim(answer)
      print*
      print*,' default value used'
      read(default,*)input
      sve_string=default
c
      return
      end subroutine rd_hidden_string
c
      end module mod_usr_hidden
c
