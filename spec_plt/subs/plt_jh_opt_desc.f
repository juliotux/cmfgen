	SUBROUTINE PLT_JH_OPT_DESC
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
	CHARACTER(LEN=1) KEY
!
	WRITE(6,'(X,A)')' '
	WRITE(6,'(X,A)')' '
	WRITE(6,'(X,A)')'General listing of options in PLT_JH. Warning: Code originally developed'
	WRITE(6,'(X,A)')'to plot J but later extended to plot H, Eta, and Chi. Be carefull with'
	WRITE(6,'(X,A)')'labeling and units. Another possible trap is that EDDFACTOR is ouput'
	WRITE(6,'(X,A)')'on the FINE grid, where as H, Eta and Chi are ouput on the regular grid.' 
	WRITE(6,'(X,A)')'If comparing J, H etc, it may be better to use the EDDFACTOR file from the'
	WRITE(6,'(X,A)')'CMFGEN computation, or alternatively, do not add extra points in the CMF_FLUX'
	WRITE(6,'(X,A)')'computation.'

	WRITE(6,'(X,A)')'Options are ordered under subject. Associated with each option are requested' 
	WRITE(6,'(X,A)')'inputs. Some inputs are not prompted for, and can only be changed from their'
	WRITE(6,'(X,A)')'default values by specifying them in the call. eg.'
	WRITE(6,'(X,A)')' '
	WRITE(6,'(X,A)')'       RD_MOD(OVER=T)'
	WRITE(6,'(X,A)')' '
	WRITE(6,'(X,A)')'Such inputs are placed in [ ].'
	WRITE(6,'(X,A)')' '

	WRITE(6,'(X,A)')'SVE and BOX file:'
	WRITE(6,'(8X,A)')'Default is to write file MAIN_OPT_STR.sve'
	WRITE(6,'(8X,A)')'Append sve=filename to write a new .sve file (no brackets)'
	WRITE(6,'(8X,A)')'Type box=filename to write a .box file containing several .sve files'
	WRITE(6,'(8X,A)')'Type .filename to read .sve file'
	WRITE(6,'(8X,A)')'Type #filename to read .box file'
!
        WRITE(6,'(A)')' '
	WRITE(6,'(3A)')RED_PEN,'Hit return to continue',DEF_PEN
        READ(5,'(A)')KEY
	WRITE(6,'(3A)')RED_PEN,'Input:',DEF_PEN
        WRITE(6,'(A)')'    RD_MOD  	Read in new data into MODEL B. Any previous model B data is lost.'
        WRITE(6,'(A)')'                     Data may be J, H, Eta, or CHI.'
	WRITE(6,'(A)')'    MOD_A	Switch to model A as the default. B data is not destroyed.'
	WRITE(6,'(A)')'    MOD_B	Switch to model B as the default. A data is not destroyed.'
!
	WRITE(6,'(A)')' '
	WRITE(6,'(3A)'),RED_PEN,'Output',DEF_PEN
	WRITE(6,'(A)')'   EXTJ          Extends J (from A model) in R, creating a new data file J_DATA'
	WRITE(6,'(A)')'                    and its corresponding file J_DATA_INFO. These can be used to'
	WRITE(6,'(A)')'                    help facilitate convergence when extending model to much'
 	WRITE(6,'(A)')'                    larger radii or when adding extra data points.'
	WRITE(6,'(A)')'   WR_ID         Writes out file names of files read.'
!
        WRITE(6,'(A)')' '
	WRITE(6,'(3A)')RED_PEN,'Hit return to continue',DEF_PEN
        READ(5,'(A)')KEY
	WRITE(6,'(3A)')RED_PEN,'X Axis Options',DEF_PEN
	WRITE(6,'(A)')'    LX		Switch between logarithmic and linear X-AXIS'
	WRITE(6,'(A)')'    LOGX        (i.e. do opposite to current setting)'
	WRITE(6,'(A)')'    LINX'
	WRITE(6,'(A)')'    XU		Change units for X-axis. Options are:'
    	WRITE(6,'(A)')'          Hz'
	WRITE(6,'(A)')'          um'
	WRITE(6,'(A)')'          Ang'
	WRITE(6,'(A)')'          keV'
	WRITE(6,'(A)')'          eV'
	WRITE(6,'(A)')'          km/s'
!
        WRITE(6,'(A)')' '
	WRITE(6,'(3A)')RED_PEN,'Y Axis Options',DEF_PEN
	WRITE(6,'(A)')'    LY		Switch between logarithmic and linear Y-AXIS'
	WRITE(6,'(A)')'    LOGY         (i.e. do opposite to current setting)'
	WRITE(6,'(A)')'    LINY'
	WRITE(6,'(A)')'    YU	        Change Y plot unit. Options are:'
	WRITE(6,'(A)')'          NAT         Default - avoids stuffing up ETA etc.'
	WRITE(6,'(A)')'          Flam        ergs/cm^2/s/Ang)'
	WRITE(6,'(A)')'          FNU         Jy'
	WRITE(6,'(A)')'          NU_FNU      ergs/cm^2/s'
!
        WRITE(6,'(A)')' '
	WRITE(6,'(3A)')RED_PEN,'Hit return to continue',DEF_PEN
        READ(5,'(A)')KEY
	WRITE(6,'(3A)')RED_PEN,'Main options',DEF_PEN
	WRITE(6,'(A)')'    JD           Plot variable at a given depth as a function of Lambda'
	WRITE(6,'(A)')'                    DEPTH   -- Depth index'
	WRITE(6,'(A)')'    RSQJD        Plots RSQJ as a function of Lambda.'
	WRITE(6,'(A)')'                    DEPTH   -- Depth index'
	WRITE(6,'(A)')'                    [SCALE] -- Scaling factor for plot.'
	WRITE(6,'(A)')'                    [ZEROV] -- Shifts plot to observers frame (velocity shift only)'
	WRITE(6,'(A)')'    R3J          Plots Int r^3.J dv as a function of depth.'
	WRITE(6,'(A)')'    JNU          Plots J at a given wavelength'
	WRITE(6,'(A)')'                    LAMBDA  -- Wavelength in Angstroms (-ve assumed to be Hz).'
	WRITE(6,'(A)')'                    [SCALE] -- Scaling factor for plot.'
	WRITE(6,'(A)')'                    [RSQJ]  -- Plot r^2.J?'
	WRITE(6,'(A)')'    SPHJ         Dirty option to modify J for sphericity effects.'
	WRITE(6,'(A)')'    ES           Convolves (J) with electron scattering redistribution function.'
	WRITE(6,'(A)')'    EJ           Plots radiation energy density; computes E(J) summed over all depths'
	WRITE(6,'(A)')'    BB           Plots blackbody spectrum at a given temperature.'
	WRITE(6,'(A)')'    CF           Plots the integrated flux (H) as a function of depth'
	WRITE(6,'(A)')'    CFD          Plots the integrated flux (H) as a function of frequencey at a given depth.  '
	WRITE(6,'(A)')'    PHOT         Normalized plot showing source of ionizations for 1/v^3 cross-section,'
	WRITE(6,'(A)')'    dBdR         Compute 1/3 dB/dR for comparision with CHI x FLUX (check on diffusion approximation)'
	WRITE(6,'(A)')'    INT          Compute Int J dv & Int dJ/dlnv dv, and plot as a function of depth (diagnostic)'
	WRITE(6,'(A)')'    DNU.         Plots c.dNU/NU as a function on NU (diagnostic) '
!
	RETURN
	END
