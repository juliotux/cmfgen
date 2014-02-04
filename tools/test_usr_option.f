c
      program test_usr_option
c
c Routine to test usr_option as it is used in DISPGEN
c
c Written  2/26/96  DLM
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      use mod_usr_option
      use mod_usr_hidden
c      use mod_wr_string
      implicit none
c
      real*8 :: a,g
      integer, dimension(10) :: c,f
      real*8, dimension(10) :: d,e
      integer :: b,h
      logical :: i,logic
c
      character(len=120) :: main,st_hid,str
c
 100  call sve_file('reset')
c
      call usr_option(main,'test1','','main option')
c
      if(main(1:2).eq."ex")then
        goto200
      elseif(main(1:4).eq."box=")then
        call wr_box_file(main)
      else
        call usr_option(a,'real','10.','sub-option 1')
        call usr_option(b,'int','20','sub-option 2')
        call usr_hidden(g,'hid_real','20.','hidden real')
        call usr_hidden(h,'hid_int','20','hidden int')
        call usr_option(c,10,3,'i_levels','1,2,3,4,5','levels')
        call usr_option(d,10,3,'r_levels','1.,2.,3.','levels')
        call usr_hidden(e,10,3,'real_hid','20.,30.,40.',
     *       'real_hid')
        call usr_hidden(f,10,3,'int_hid','20,30,40',
     *       'int_hid')
        call usr_option(logic,'logic','False','logical sub-option')
        call usr_option(str,'str','Hello','string sub-option')
        call usr_hidden(i,'log_hid','FALSE',
     *       'log_hid')
        call usr_hidden(st_hid,'st_hid','JUNK IT',
     *       'st_hid')
      endif
c
      goto 100
c
 200  continue
c
      stop
      end
c
