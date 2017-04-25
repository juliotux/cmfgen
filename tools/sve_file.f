c
c       **** SUBROUTINE SVE_FILE ****
c
c This subroutine is used by the USR_OPTION routines to read and write
c .sve files.  There are 6 subroutine, 4 that are called only by the
c usr_option routine and 2 (sve_file and wr_box_file) which are 
c called by the routine that calls usr_option.
c
c This set of subroutines was written to be as self contained as possible.
c To use USR_OPTION simply include a USE MOD_USR_OPTION as the first line
c of any routine and link to the libtools library.
c
c------------------------------------------------------------------------------
c
c    Definitions:
c       Main option => This is the top of the input
c                      structure.  For example, if_carb,
c                      dc_feiv.  The main option is written as
c                      the first entry on a new line in the .sve
c                      file without any flags.  Each main
c                      option may or may not have sub options.
c    
c        Sub option => These are options that query for specific
c                      information about the given main option.
c                      They are written in the .sve files on the
c                      same line as its main option, deliminted with
c                      [var_name=sub option]
c    
c     Hidden option => These are options where the default value passed
c                      to usr_hidden is always used unless the user
c                      indicates otherwise by including the hidden option
c                      on the main option input line (ie rd_obs (over=T)).
c
c          sve file => This is a .sve file with a single main option
c                      and it's subsequent sub options.  The name of the
c                      .sve file the name of the main option with the suffix
c                      .sve (ie if_iron.sve) unless a sve= is appended to
c                      the main option and them the name is "file".sve
c                      (ie "if_iron sve=junk" => junk.sve).
c
c          box file => This is a file with a list of .sve files to process.
c                      Each line contains the name of the .sve file (ie junk.sve)
c                      This .box file must be created by the program that
c                      calls usr_option.
c
c          var_name => This is the generic name sent in the USR_OPTION
c                      call that will be used to denote the current option.
c                      var_name does not have to be the name of the varaible
c                      being read and can contain any "normal" characters
c                      (ie commas, periods etc) except quotes and double
c                      quotes.
c
c         .filename => Indicates usr_option should read the previously written
c                      file filename.sve and execute the main option and
c                      sub-options listed there.  A sub-option can be changed
c                      including the option to change and the new value along
c                      with the filename (ie .filename (over=T))
c
c         #filename => Indicated .box file should be read which contains a list
c                      of .sve files to execute. The list of .sve files can
c                      contain hidden options and changes to option in the
c                      individual .sve files
c
c------------------------------------------------------------------------------
c
c  Subroutines:
c  POST_SVE_FILE: This subroutine is called after the user inputs an
c                 answer.  It only acts if a main option is being
c                 input and if a .sve or .box file is to be read.
c                 It also handles changes to .sve files.
c  
c   PRE_SVE_FILE: This subroutine acts before the user inputs an answer
c                 If a .sve file is being read, the current option (main or sub)
c                 is read from the .sve file and the user is not prompted.  If the
c                 option is not found, the user is then prompted for an 
c                 answer.
c  
c    WR_SVE_FILE: If a .sve file is being written, this routine writes the
c                 last answer to the current .sve file.  It determines if
c                 the option is main or sub and writes it in the correct
c                 format.
c  
c       SVE_FILE: This subroutine is called to set several flags.  The only
c                 passed string currently used is 'RESET'.  This resets
c                 several flags and if a .box file is being read, the next
c                 .sve file is opened and read.
c
c  OPEN_SVE_FILE: Subroutine that determines the name of the .sve or .box file and opens
c                 it.  If it is sent just the name of a main option, a new file named
c                 "main_option.sve" is opened.  If sve=filename is apended to the
c                 main option then a new file named "filename.sve" is opened.  If
c                 .filename is sent then an old file named "filename.sve" is opened.  Finally,
c                 if box=filename is passed then a .box file is opened or if #filename
c                 then an old file named "filename.box" is opened
c
c    WR_BOX_FILE: Subroutine that creates a .box file.
c
c------------------------------------------------------------------------------
c
c   Creation of .sve file:
c     Creation of .sve files is carried out automatically.  Each main
c       option creates a file names "main_option.sve" which will contain one line
c       with the main option listed first and all subsequent sub-options listed
c       as [var_name=sub-option].
c     If sve="filename" is appended to the main option, the name of the save file will
c       be the name give (ie "if_carb sve=junk" creates junk.sve containing the
c       main option if_carb and its sub-options).
c
c   Read options from .sve file:
c     To read a .sve file, preceed the name of the .sve file with a "."
c     (ex =>.if_carb  ; reads if_carb.sve file).  Changes are made to sub-options
c     by including (sub-option=new_value) on the line with .filename.
c
c   Read options from .box file:
c     To read a .box file, preceed the name of the .sve file with a "#"
c     (ex =>#if_carb  ; reads if_carb.box file).  Changes to .sve files
c     must be done when the .box fie is created.
c
c   Creation of .box file:
c     Append box="filename" to the main option answer to create a .box file called
c     filename.box.  The user will be prompted for the names of .sve files to include
c     in the .box file.  The .sve files must exist in order to put them in a .box
c     file.  A carriage return signals that the .box is complete.  Changes to the
c     sub-options in the .sve file can be included (ie filename (levs=1,3,5,7))
c
c   Editing of .sve file:
c     .sve and .box files are editable with a text editor.  Be sure to follow the format
c     for main and sub options described above.  If a sub-option is not present in a
c     .sve file, the user is prompted for the options.  In this way, the .sve file can be
c     edited and a specific sub-option removed to force the user to input the value (ie
c     a smoothing parameter or scale factor).
c
c------------------------------------------------------------------------------
c
c To use USR_OPTION, a call to SVE_FILE('reset') is required before each main
c option.  Then next call to USR_OPTION is then the main option and all
c subsequent calls to USR_OPTION and USR_HIDDEN are sup-options, until a call to 
c SVE_FILE('reset'). This way all sub-options are written in the same save file.  
c If .box files are to be written, an IF branch to call WR_SVE_FILE(MAIN) must be
c created.
c 
c------------------------------------------------------------------------------
c
c General coding to use .sve  and .box files
c
c 100  call sve_file('reset')
c      ....
c c       main option
c      call usr_option(main,'var_name','default','description')
c      ....
c      if(main.eq.option)then
c c         sub options     
c        call usr_option(variable,'var_name','default','description')
c        call usr_option(variable,'var_name','default','description')
c        call usr_hidden(variable,'var_name',main_option,'default','description')
c        ...
c      elseif(main.eq.???)then
c c         sub options     
c      ...
c      elseif(main.eq.'box=')then
c c         write .box file
c        call wr_box_file(main)
c      elseif(main.eq.'exit')then
c        goto 200
c      endif
c      ....
c      goto 100
c      ....
c 200  continue
c      ....
c
c See usr_option.f for more information on USR_OPTION
c
c       **** SUBROUTINE SVE_FILE ****
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c written 12/9/96  DLM
c
c altered  2/26/96 DLM Changed entry points to modules.
c                      Removed extended sve files and replaced with
c                      box files.
c
c*****************************************************************
c
      module mod_sve_file
c
c*****************************************************************
c
      integer, parameter :: unit_sve=64
      integer, parameter :: unit_box=65
c
      character(len=500) string
c
      logical :: sve_read,sve_write,main_option,box_read
c
      end module mod_sve_file
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine pre_sve_file(answer,var_name,present)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c This subroutine checks to see if a .sve file is being read.  It is
c the first action taken in the usr_option subroutines and check to
c see if the current option should be read from file or if the usr
c should be prompted.
c 
c 1) If a .sve file is being read, then the next option, main or sub, 
c    is read and passed back to usr_option in answer and present
c    is set to .TRUE.  
c
c 2) If a .sve file is not being read or the option described by "var_name"
c    is not found, present is set to .FALSE. and usr_option prompts the user
c    for the option.
c
c written  12/9/96  DLM
c
c altered  2/27/97  DLM  Removed extended .sve files and include .box
c                          file.  Not many changes here, but see routine
c                          post_sve_file.
c                        Name changed to post_sve_file
c
c altered  3/17/97  DLM  Use "{}" instead of "[]" to delimit sub-options
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      use mod_sve_file
      implicit none
c
      character(len=*) :: answer,var_name
c
      logical :: present
c
      integer :: l,k,o_end,o_start
c
      character(len=120) :: option
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      if(sve_read)then
c
c Operate only if options available in .sve file
c
        if(main_option)then
c
c If main option then read first string
c
          l=index(string,' ')-1
          k=index(string,'(')-1
          if((k.gt.0).and.(k.lt.l))l=k
          answer=string(:l)
          write(*,"(6x,'Main option=',a)")trim(answer)
          present=.TRUE.
c
        else
c
c If sub-option then search for "var_name"
c
          option="{"//trim(var_name)//"="
          l=index(string,trim(option))
          if(l.gt.0)then
c
c If found then read value into answer which will be read in usr_option
c
            o_start=l+len_trim(option)
            o_end=o_start+index(string(o_start:),'}')-2
            if(o_end.lt.o_start)then
              print*,' apparent error in option format'
              print*,' could not find `}`'
              present=.FALSE.
            else
              answer=string(o_start:o_end)
              present=.true.
              write(*,"(8x,a,a)")trim(option(2:)),trim(answer)
            endif
          else
c
c If not found then usr_option will prompt user
c
            present=.false.
c
          endif
c
        endif
c
      else
c
        present=.false.
c
      endif
c
      return
      end subroutine pre_sve_file
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine post_sve_file(answer,var_name,present)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c This subroutine checks the user given answer to see if a .box or .sve
c file should be opened and/or read.  It is called by usr_option after
c the user has entered an option.  It is not called if a .sve file is
c being read.
c
c 1) If "answer" starts with "E", "EX" or "GR" then return to usr_option.
c
c 2) If "answer" is blank then return to usr_option, which will read
c    default answer.
c
c 3) If "answer" begins with a "." then the .sve file is opened, options
c    (sub and main) are read from the .sve file and main option is read
c    into "answer".  The sub-options are read into "answer" by subsequent
c    call to usr_option via the subroutine pre_sve_file.
c
c 4) If "answer" begins with a "#" then the .box file is opened, the name
c    of the first .sve file is read, the .sve file is opened, options
c    (sub and main) are read from the .sve file and main option is read
c    into "answer".  The sub-options are read into "answer" by subsequent
c    call to usr_option via the subroutine pre_sve_file.
c
c 5) If sve=filename is found in "answer" then .sve file is opened with the name
c    of the option passed in "answer" (ie "filename.sve").  Options are written 
c    to the .sve file by the subroutine wr_sve_file called by usr_option.
c
c 6) If box=filename is found in "answer" then return to usr_option.  Create of
c    .box file must be handled by the routine that called usr_option via a call
c    to wr_box_file.
c
c 7) Else .sve file is opened with the name of the option passed in "answer"
c    (ie "option.sve").  Options are written to the .sve file by the
c    subroutine wr_sve_file called by usr_option.  If hidden options are
c    included with the main option, the new values are stored in the variable
c    "string" and later read by usr_hidden.
c
c written  12/9/96  DLM
c
c altered  2/27/97  DLM  Removed extended .sve files and include .box file.
c                          Required the additions that subroutine open_sve_file
c                          would open both .sve or .box files.
c                        Name changed to post_sve_file.
c
c altered  3/17/97  DLM  Use "{}" instead of "[]" to delimit sub-options.
c                        Use upper_case instead of condit_sting to create
c                          all capital letter strings.
c                        Changed .sve files to have one option per line.
c                          Required all lines to be read in and concatinated
c                          into variable string.  String is still checked 
c                          elsewhere for sub-options and hidden options.
c                        Hidden options included with a main options are now
c                          written to variable string.  Previously the main option
c                          had to be passed in the call to usr_option.  Now
c                          the calls to usr_option and _hidden are the same.
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      use mod_sve_file
      implicit none
c
      character(len=*) :: answer,var_name
c
      logical :: present
c
      integer :: s,b,l,k,i,o_start,o_end,num_hidden
      integer :: ios
c
      character(len=120) :: all_caps,option,next_option
c
      integer, parameter :: max_hidden=10
      character(len=120), dimension(max_hidden) :: new_var,new_value
c
      external upper_case
      character(len=120) :: upper_case      !sets all letters to upper case
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      present=.TRUE.
c
      all_caps=upper_case(answer)
c
      if(main_option)then
        s=index(all_caps,'SVE=')
        b=index(all_caps,'BOX=')
c
        if(all_caps(1:2).eq.'EX')then
c
c 1) Return to usr_option to execute this option
c
        elseif(len_trim(answer).eq.0)then
c
c 2) Return to usr_option read in default value
c
        elseif(answer(1:1).eq.'.')then
c
c 3) Open previously written .sve file and read first option
c
          call open_sve_file(answer,present)
          if(.not.present)then
            print*,' error trying to open .sve file'
            print*,' file does not exist'
            return
          endif
c
c    Check for hidden options or changes to sub-options in .sve file given
c    by user in "answer"
c
          k=1
          do i=1,max_hidden
            o_start=k+index(answer(k:),'(')
            if(o_start.le.k)then
c
c    No more hidden options or changes
c                         
              num_hidden=i-1
              exit
c
c    Determine hidden option name or change
c
            else
              o_end=o_start+index(answer(o_start:),'=')-1
              if(o_end.le.o_start)then
                print*,' apparent error in option format'
                print*,' could not find `=` in input string'
                num_hidden=i-1
                exit
              endif
              new_var(i)=upper_case(answer(o_start:o_end))
c
c    Determine new value for hidden option or change
c
              o_start=o_end+1
              o_end=o_start+index(answer(o_start:),')')-2
              if(o_end.lt.o_start)then
                print*,' apparent error in option format'
                print*,' could not find `)` in input string'
                num_hidden=i-1
                exit
              endif
              new_value(i)=answer(o_start:o_end)
              k=o_end+2
            endif
          enddo
c
c   Read options from .sve file, which consists of one line for each
c   option, into the variable "string".  String is used to store the
c   values of the sub-options, which are read later by usr_option via
c   pre_sve_file.
c
          read(unit_sve,"(a)",iostat=ios)string
          if(ios.ne.0)then
            print*,' error reading .sve file ',string
            close(unit_sve)
            close(unit_box)
            box_read=.FALSE.
            sve_read=.FALSE.
          else
            sve_read=.TRUE.
          endif
          read(unit_sve,"(a)",iostat=ios)next_option
          do while(ios.eq.0)
            l=len_trim(string)+2
            k=len_trim(next_option)
            string(l:)=next_option(:k)
            read(unit_sve,"(a)",iostat=ios)next_option
          enddo
	  rewind(unit_sve)    					!Jdh 21-Mar-1997
c
c   Read main_option from string into answer which will
c   be read by usr_option.
c
          l=index(string,' ')-1
          k=index(string,'(')-1
          if((k.gt.0).and.(k.lt.l))l=k
          answer=string(:l)
          write(*,"(6x,'Main option=',a)")trim(answer)
          sve_read=.TRUE.
c
c    If there are hidden options or changes put new_values into string
c
          if(num_hidden.gt.0)then
            do i=1,num_hidden
              option="{"//trim(new_var(i))
              o_start=index(string,trim(option))
c
c    If hidden option already in string or change then replace value
c
              if(o_start.gt.0)then
                o_start=o_start+len_trim(option)-1
                o_end=o_start+index(string(o_start:),'}')-1
                if(o_end.lt.o_start)then
                  print*,' apparent error in option format'
                  print*,' could not find `}` in .sve file'
                  print*,string
                  exit
                endif
                string=string(:o_start)//trim(new_value(i))//
     *               string(o_end:)
c
c    Else put hidden option at the end of string
c
              else
                o_start=len_trim(string)+1
                string=string(:o_start)//trim(option)//
     *               trim(new_value(i))//"}"
              endif
            enddo
          endif
c
        elseif(answer(1:1).eq.'#')then
c
c 4) Open previously written .box file, read name of first .sve file,
c    open this .sve file and read main option
c
          call open_sve_file(answer,present)
          if(.not.present)then
            print*,' error trying to open .box file'
            print*,' file does not exist'
            return
          endif
          box_read=.TRUE.
c
c    Read name of .sve file from .box file
c
          read(unit_box,"(a)")string
c
c    Open previously written .sve file
c
          if(string(1:1).ne.'.')string="."//string
          call open_sve_file(string,present)
          if(.not.present)then
            print*,' error trying to open .sve file'
            print*,' file does not exist'
            present=.FALSE.
            return
          endif
c
c    Check for hidden options or changes to sub-options in .sve file
c
          k=1
          do i=1,max_hidden
            o_start=k+index(string(k:),'(')
            if(o_start.le.k)then
c
c    No more hidden options or changes
c
              num_hidden=i-1
              exit
c
c    Determine hidden option name or change
c
            else
              o_end=o_start+index(string(o_start:),'=')-1
              if(o_end.le.o_start)then
                print*,' apparent error in option format'
                print*,' could not find `=` in box string'
                print*,' '//trim(string)
                num_hidden=i-1
                exit
              endif
              new_var(i)=upper_case(answer(o_start:o_end))
c
c    Determine new value for hidden option or change
c
              o_start=o_end+1
              o_end=o_start+index(string(o_start:),')')-2
              if(o_end.lt.o_start)then
                print*,' apparent error in option format'
                print*,' could not find `)` in box string'
                print*,' '//trim(string)
                num_hidden=i-1
                exit
              endif
              new_value(i)=string(o_start:o_end)
              k=o_end+2
            endif
          enddo
c
c   Read options from .sve file, which consists of one line for each
c   option, into the variable "string".  String is used to store the
c   values of the sub-options, which are read later by usr_option via
c   pre_sve_file.
c
          read(unit_sve,"(a)",iostat=ios)string
          if(ios.ne.0)then
            print*,' error reading .sve file ',string
            close(unit_sve)
            close(unit_box)
            box_read=.FALSE.
            sve_read=.FALSE.
          else
            sve_read=.TRUE.
          endif
          read(unit_sve,"(a)",iostat=ios)next_option
          do while(ios.eq.0)
            l=len_trim(string)+2
            k=len_trim(next_option)
            string(l:)=next_option(:k)
            read(unit_sve,"(a)",iostat=ios)next_option
          enddo
	  rewind(unit_sve)    					!Jdh 19-Feb-2001
c
c   Read main_option from string into answer which is will
c   be read by usr_option.
c
          l=index(string,' ')-1
          k=index(string,'(')-1
          if((k.gt.0).and.(k.lt.l))l=k
          answer=string(:l)
          write(*,"(6x,'Main option=',a)")trim(answer)
          sve_read=.TRUE.
c
c   If there are hidden options or changes put new_values into string
c
          if(num_hidden.gt.0)then
            do i=1,num_hidden
              option="{"//trim(new_var(i))
              o_start=index(string,trim(option))
c
c   If hidden option already in string or change then replace value
c
              if(o_start.gt.0)then
                o_start=o_start+len_trim(option)-1
                o_end=o_start+index(string(o_start:),'}')-1
                string=string(:o_start)//trim(new_value(i))//
     *               string(o_end:)
c
c   Else put hidden option at the end of string
c
              else
                o_start=len_trim(string)+1
                string=string(:o_start)//trim(option)//
     *               trim(new_value(i))//"}"
              endif
            enddo
          endif
c
        elseif(s.gt.0)then
c
c 5) If sve=filename found then open new file named "filename.sve"
c
          call open_sve_file(answer,present)
          sve_read=.FALSE.
          present=.TRUE.
c
        elseif(b.gt.0)then
c
c 6) If box= found then return to usr_option.  Creation of .box files handled
c    by routine that call usr_option.
c
        elseif(b.eq.0)then
c
c 7) If sve= and box= not given then open file "option".sve
c
          call open_sve_file(answer,present)
          sve_read=.FALSE.
          present=.TRUE.
c
c    If hidden options included with main_option then write them to variable
c    string.  They will be read later by usr_hidden.
c
          k=1
          l=1
          do i=1,max_hidden
            o_start=k+index(answer(k:),'(')
            if(o_start.le.k)then
c
c    No more hidden options or changes
c
              num_hidden=i-1
              exit
c
c    Put hidden option into string.
c
            else
              o_end=o_start+index(answer(o_start:),'=')-1
              if(o_end.le.o_start)then
                print*,' apparent error in option format'
                print*,' could not find `=` in box string'
                print*,' '//trim(answer)
                num_hidden=i-1
                exit
              endif
              string(l:)="{"//upper_case(answer(o_start:o_end))
              l=len_trim(string)+1
              o_start=o_end+1
              o_end=o_start+index(answer(o_start:),')')-2  
              if(o_end .lt. o_start)then			!le to lt jdh
                print*,' apparent error in option format'
                print*,' could not find `)` in box string'
                print*,' '//trim(answer)
                num_hidden=i-1
                exit
              endif
              string(l:)=answer(o_start:o_end)//"}"
              k=o_end+1
              l=len_trim(string)+2
            endif
          enddo
c
        else
c
c 8) Should never get here!!!
c
         print*,' what to do, what to do'
          stop
c
        endif
      endif
c
      return
      end subroutine post_sve_file
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine sve_file(answer)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c This subroutine checks and resets all the sve_file flags.  If a .box
c file is being read, then the next .sve file name in .box file is read,
c the .sve file is opened and the first line is read into "string".  The
c options are read from "string" by pre_sve_file.
c
c written  12/9/96  DLM
c
c altered  2/27/97  DLM  Removed extended .sve files and include .box file.
c
c altered  3/17/97  DLM  Use "{}" instead of "[]" to delimit sub-options.
c                        Use upper_case instead of condit_sting to create
c                        all capital letter strings.
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      use mod_sve_file
      implicit none
c
      integer :: ios,l,i,k,num_hidden,o_start,o_end
c
      logical :: present
c
      character(len=*) :: answer
      character(len=120) :: all_caps,option,next_option
c
      integer, parameter :: max_hidden=10
      character(len=120), dimension(max_hidden) :: new_var,new_value
c
      external upper_case
      character(len=120) :: upper_case      !sets all letters to upper case
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      all_caps=upper_case(answer)
c
      if(trim(all_caps).eq.'RESET')then
        close(unit_sve)
        sve_read=.FALSE.
        main_option=.TRUE.
        if(box_read)then
          read(unit_box,"(a)",iostat=ios)string
          if(ios.ne.0)then
            close(unit_box)
            box_read=.FALSE.
            sve_read=.FALSE.
          else
            call open_sve_file(string,present)
            if(.not.present)then
              print*,' .sve file not found =>'//trim(string)
              print*,' .box files closed'
              close(unit_box)
              sve_read=.FALSE.
              box_read=.FALSE.
            endif
c
c    Check for hidden options or changes to sub-options in .sve file
c
            k=1
            do i=1,max_hidden
              o_start=k+index(string(k:),'(')
              if(o_start.le.k)then
c
c    No more hidden options or changes
c
                num_hidden=i-1
                exit
c
c    Determine hidden option name or change
c
              else
                o_end=o_start+index(string(o_start:),'=')-1
                if(o_end.le.o_start)then
                  print*,' apparent error in option format'
                  print*,' could not find `=` in box string'
                  num_hidden=i-1
                  exit
                endif
                new_var(i)=upper_case(string(o_start:o_end))
c
c    Determine new value for hidden option or change
c
                o_start=o_end+1
                o_end=o_start+index(string(o_start:),')')-2
                if(o_end.lt.o_start)then
                  print*,' apparent error in option format'
                  print*,' could not find `)` in box string'
                  print*,' '//trim(string)
                  num_hidden=i-1
                  exit
                endif
                new_value(i)=string(o_start:o_end)
                k=o_end+2
              endif
            enddo
c
c   Read main_option from .sve file and put it into default, which will
c   be read by usr_option
c
c   Read options from .sve file, which consists of one line for each
c   option, into the variable "string".  String is used to store the
c   values of the sub-options, which are read later by usr_option via
c   pre_sve_file.
c
          read(unit_sve,"(a)",iostat=ios)string
          if(ios.ne.0)then
            print*,' error reading .sve file ',string
            close(unit_sve)
            close(unit_box)
            box_read=.FALSE.
            sve_read=.FALSE.
          else
            sve_read=.TRUE.
          endif
          read(unit_sve,"(a)",iostat=ios)next_option
          do while(ios.eq.0)
            l=len_trim(string)+2
            k=len_trim(next_option)
            string(l:)=next_option(:k)
            read(unit_sve,"(a)",iostat=ios)next_option
          enddo
	  rewind(unit_sve)    					!Jdh 19-Feb-2001
c
c   If there are hidden options or changes put new_values into string
c
            if(num_hidden.gt.0)then
              do i=1,num_hidden
                option="{"//trim(new_var(i))
                o_start=index(string,trim(option))
c
c   If hidden option already in string or change then replace value
c
                if(o_start.gt.0)then
                  o_start=o_start+len_trim(option)-1
                  o_end=o_start+index(string(o_start:),'}')-1
                  string=string(:o_start)//trim(new_value(i))//
     *                 string(o_end:)
c
c   Else put hidden option at the end of string
c
                else
                  o_start=len_trim(string)+1
                  string=string(:o_start)//trim(option)//
     *                 trim(new_value(i))//"}"
                endif
              enddo
            endif
c
          endif
	else
	  string=' '
        endif
      else
        print*,'unknown string passed to sve_file =>',answer
        stop
      endif
c
      return
      end subroutine sve_file
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine wr_sve_file(answer,var_name)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c This subroutine is used by the usr_option routines to write
c the option given in the current .sve file.
c
c written  12/9/96  DLM
c
c altered  2/27/97  DLM  Included the use of all_caps.
c                        Name changed to wr_sve_file.
c
c altered  3/17/97  DLM  Use "{}" instead of "[]" to delimit sub-options.
c                        Use upper_case instead of condit_sting to create
c                        all capital letter strings.
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      use mod_sve_file
      implicit none
c
      character(len=*) :: answer,var_name
c
      integer :: o_end,l
c
      character(len=120) :: fmt_string,option,all_caps
c
      external upper_case
      character(len=120) :: upper_case
c
      external wr_fmt_string
      character(len=120) :: wr_fmt_string
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      all_caps=upper_case(answer)
c
C      if(    (all_caps(1:2).eq.'EX')
C     *   .or.(all_caps(1:1).eq.'X')
C     *   .or.(all_caps(1:4).eq.'BOX='))then
      if( all_caps(1:4).eq.'BOX=')then
        return
C      elseif(.not.sve_read)then
       else
        o_end=index(answer,' ')-1
        l=index(answer,'(')-1
        if((l.gt.0).and.(l.lt.o_end))o_end=l
        if(main_option)then
          option=answer(:o_end)
        else
          option="{"//trim(var_name)//"="//answer(:o_end)//"}"
          o_end=len_trim(option)
        endif
        fmt_string=wr_fmt_string(o_end)
        write(unit_sve,fmt_string)
     *       option(:o_end)
      endif        
      main_option=.FALSE.
c
      return
      end subroutine wr_sve_file
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine open_sve_file(answer,present)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Subroutine to determine the name of a .sve or .box file
c and opens it.
c
c written  12/9/96  DLM
c
c altered  2/27/97  DLM  Removed extended .sve files and include .box file.
c                        Name changed to open_sve_file
c
c altered  3/17/97  DLM  Change open() to generic gen_asci_open
c                        Use upper_case instead of condit_string to create all
c                        capital letter string
c
c altered 27/01/04  DJH  Installed FILE_SHOULD_EXIT variable so that routine
c                        does not create a sve file when using the .option.
c                        Designed to prevent the creation of empty .sve files.
c 
c altered 08/06/15  DJH  .box file no longer created if it does not exist.
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      use mod_sve_file
      implicit none
c
      integer :: o_start,o_end,l,s,b,unit,ios
c
      character(len=*) :: answer
      character(len=120) :: filename,ending,all_caps
c
      logical :: file_should_exist
      logical :: present
c
      external upper_case
      character(len=120) :: upper_case
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      o_start=1
      o_end=index(answer,' ')-1
      l=index(answer,'(')-1
      if((l.gt.0).and.(l.lt.o_end))o_end=l
      all_caps=upper_case(answer)
      s=index(all_caps,'SVE=')
      b=index(all_caps,'BOX=')
c
c Determine start and end of filename to open and
c the ending (.sve or .box)
c
      file_should_exist=.false.
      if(answer(1:1).eq.'.')then
        o_start=2
        ending=".sve"
        unit=unit_sve
        file_should_exist=.true.
      elseif(answer(1:1).eq.'#')then
        o_start=2
        ending=".box"
        unit=unit_box
        file_should_exist=.true.
      elseif(s.ge.1)then
        o_start=s+4
        o_end=o_start+index(answer(o_start:),' ')-2
        l=o_start+index(answer(o_start:),'(')-2
        if((l.gt.o_start-1).and.(l.lt.o_end))o_end=l
        l=o_end-o_start+1
        if(l.gt.9)then
          o_end=o_start+8
          print*,' sve file name too long, truncated to ',
     *         answer(o_start:o_end)
        endif
        ending=".sve"
        unit=unit_sve
      elseif(b.ge.1)then
        o_start=b+4
        o_end=o_start+index(answer(o_start:),' ')-2
        l=o_start+index(answer(o_start:),'(')-2
        if((l.gt.o_start-1).and.(l.lt.o_end))o_end=l
        l=o_end-o_start+1
        if(l.gt.9)then
          o_end=o_start+8
          print*,' box file name too long, truncated to ',
     *         answer(o_start:o_end)
        endif
        ending=".box"
        unit=unit_box
      else
        ending=".sve"
        unit=unit_sve
      endif
      if(o_end-o_start+1.gt.3)then
        if(index(all_caps(o_start:o_end),'.SVE').gt.0)then
          filename=answer(o_start:o_end)
        elseif(index(all_caps(o_start:o_end),'.BOX').gt.0)then
          filename=answer(o_start:o_end)
        else
          filename=answer(o_start:o_end)//trim(ending)
        endif
      else
        filename=answer(o_start:o_end)//trim(ending)
      endif
!
      inquire(file=filename,exist=present)
      if(file_should_exist .and. .not. present)then
        return
      end if
c
      call gen_asci_open(unit,filename,'unknown',' ',' ',0,ios)
      if(ios.ne.0)
     *     print*,' error opening file ',trim(filename)
      return
c
      end subroutine open_sve_file
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine wr_box_file(filename)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Subroutine to create a .box file.
c
c written  2/27/97  DLM
c
c altered  3/17/97  DLM Use upper_case instead of condit_string to create
c                       all capital letter string
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      use mod_sve_file
      implicit none
c
      integer :: l,m
c
      character(len=*) :: filename
      character(len=120) :: answer,all_caps,svename
c
      logical :: present
c
      external upper_case
      character(len=120) :: upper_case
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      call open_sve_file(filename,present)
c
      l=1
      do while(l.gt.0)
c
        write(*,"(1x,'Enter .sve file name [cr to end]=>')",
     *       advance="no")
        read(*,"(a)")answer
        l=len_trim(answer)
        if(l.gt.0)then
          l=index(answer,' ')-1
          m=index(answer,'(')-1
          if((m.gt.0).and.(m.lt.l))l=m
          all_caps=upper_case(answer)
          m=index(all_caps(:l),'.SVE')
          if(m.le.0)svename=answer(:l)//".sve"
          if(svename(1:1).eq.'.')svename=svename(2:)
          inquire(file=svename,exist=present)
          if(.not.present)then
            print*,' file not found: ',trim(svename)
          else
            write(unit_box,"(a)")trim(answer)
          endif
        endif
c
      enddo
c
      close(unit_box)
c
      return
c
      end subroutine wr_box_file
c
