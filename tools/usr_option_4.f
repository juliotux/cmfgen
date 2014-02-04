c
c    **** MODULE USR_OPTION_4 ***
c
c This module creates the generic subroutine called USR_OPTION.
c When a call is made to USR_OPTION, the interface determines
c which subroutine to call by matching the types and sizes of
c variables passed in the call.
c
c DISPGEN and subroutines require double precision (*8) where PGPLOT
c routines require single precision (*4).  Thus, this module is compiled
c without the -r8 compiler flag.  These (*4) routines are
c archived in the library usr_option/libusr.a.  This
c then requires all the DISPGEN routines to only USE MOD_USR_OPTION.  PGPLOT,
c however, must also USE MOD_USR_OPTION_4.
c
c The possible subroutines are:
c
c   RD_REAL_4         => read one real *4 variable
c   RD_REAL_DIM_4     => read an array of real *4 variable
c
c The general call to USR_OPTION is for one input variable:
c
c   CALL USR_OPTION(VARIABLE,'VAR_NAME','DEFAULT','DESCRIPTION')
c
c and for an array of real or integer variables:
c
c   CALL USR_OPTION(VARIABLE,DIMENSION,REQUIRED,'VAR_NAME',
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
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      module mod_usr_option_4
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      interface usr_option
        module procedure rd_real_4,
     *                   rd_real_dim_4
      end interface
c
      contains
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine rd_real_4(input,var_name,default,description)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c This routine reads a real*4 value into the variable input.  Normally
c the user is prompted for an answer, but if a previously written
c .sve file is being read, the option is read from file and the user
c is not prompted.
c
c  written  11/23/96  DLM
c
c  altered   1/13/97  DLM  Added external function WR_FMT_STRING
c                          for VMS's benefit
c                          Changed (a,... to (1x,a for VMS's benefit
c
c  altered   3/17/97  DLM  Now read and write var_name as all capital
c                          letters.  Cannot reassign var_name because
c                          a specific string is passed to usr_option
c                          instead of a character variable.  Thus must
c                          use another variable, all_caps, to create
c                          an all upper case string (SGI will allow
c                          change to var_name, but VMS will not!).
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      implicit none
c
      integer :: i
c
      real :: input
c
      character(len=*) :: var_name,description,default
      character(len=120) :: answer,sve_string,all_caps
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
c
c Determine if .sve file is being read.  If so, read next option,
c either a main option or a sub option.
c
      call pre_sve_file(answer,all_caps,present)
c
c If option read from .sve file, do not prompt usr.
c
      if(present)goto 400
c
c Prompt usr for input.
c
      write(*,"(3x,a)")trim(description)
 300  if(len_trim(default).ne.0)then
        write(*,"(1x,a,' [',a,']=>')",advance="no")
     *       var_name,trim(default)
      else
        write(*,"(1x,a,' =>')",advance="no")var_name
      endif
      read(*,"(a)")answer
c
c If main option, 'close' or 'file=' flag then take approiate action.
c
      call post_sve_file(answer,all_caps,present)
c
c If 'close' then prompt user for input again.
c
      if(.not.present)goto 300
c
c Read answer, from user or from .sve file, into input variable.
c
 400  if(len_trim(answer).ne.0)then
        read(answer,*,err=300)input
        sve_string=answer
      elseif(len_trim(default).ne.0)then
        read(default,*)input
        sve_string=default
      else
        write(*,"(' no default given; you must enter a value')")
        goto 300
      endif
c
c Write option to current .sve file.
c
      call wr_sve_file(sve_string,all_caps)
c
      return
c
      stop
      end subroutine rd_real_4
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine rd_real_dim_4(input,dim,required,
     *     var_name,default,description)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c This routine reads an array of real*4 values into the variable array input.
c Normally the user is prompted for an answer, but if a previously written
c .sve file is being read, the option is read from file and the user
c is not prompted.  The passed varaibles 'dim' and 'required' indicates
c the size of the array and the minimum number of variables the user must
c input respectively.
c
c  written  11/23/96  DLM
c
c  altered   3/4/97   DLM  Added do loops as inputs
c
c  altered   3/17/97  DLM  Now read and write var_name as all capital
c                          letters.  Cannot reassign var_name because
c                          a specific string is passed to usr_option
c                          instead of a character variable.  Thus must
c                          use another variable, all_caps, to create
c                          an all upper case string (SGI will allow
c                          change to var_name, but VMS will not!).
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      use mod_wr_string_4
      implicit none
c
      integer :: dim,required,found
      integer :: i,l,m,start,stop,step,count
c
      real, dimension(dim) :: input
c
      character(len=*) :: var_name,description,default
      character(len=120) :: answer,sve_string,number,all_caps
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
c
      call pre_sve_file(answer,all_caps,present)
      if(present)goto 400
      write(*,"(3x,a)")trim(description)
 300  if(len_trim(default).ne.0)then
        write(*,"(1x,a,' [',a,']=>')",advance="no")
     *       var_name,trim(default)
      else
        write(*,"(1x,a,' =>')",advance="no")var_name
      endif
      read(*,"(a)")answer
      call post_sve_file(answer,all_caps,present)
      if(.not.present)goto 300
 400  if(len_trim(answer).ne.0)then
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
            number=wr_string(float(l))
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
        if(found.lt.required)then
          write(*,"(1x,i2,' values required, only found ',i2)")
     *         required,found
          goto 300
        endif
        sve_string=answer
      else
        if(len_trim(default).ne.0)then
          read(default,*,end=200)(input(i),i=1,dim)
 200      found=i-1
          if(found.lt.dim)then
            do i=found+1,dim
              input(i)=0.0
            enddo
          endif
          sve_string=default
        else
          write(*,
     *         "(' no default given; you must enter a value')")
          goto 300
        endif
      endif
c
      call wr_sve_file(sve_string,all_caps)
c
      return
c
      stop
      end subroutine rd_real_dim_4
c
      end module mod_usr_option_4
c
