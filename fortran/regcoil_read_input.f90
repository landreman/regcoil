subroutine regcoil_read_input

  implicit none

  integer :: numargs
  character(len=200) :: inputFilename
  integer :: ios

  ! getcarg is in LIBSTELL
  call getcarg(1, inputFilename, numargs)

  if (numargs<1) then
     stop "One argument is required: the input namelist file, which must be named regcoil_in.XXXXX"
  end if
  if (numargs>1) then
     print *,"WARNING: Arguments after the first will be ignored."
  end if

  call regcoil_read_input_file(trim(inputFilename), ios)
  if (ios /= 0) stop

end subroutine regcoil_read_input


subroutine regcoil_read_input_file(inputFilename, ios)

  use regcoil_variables

  implicit none

  character(len=*), intent(in) :: inputFilename
  integer, intent(out) :: ios
  integer :: fileUnit, didFileAccessWork
  character(len=200) :: linebuf
  integer :: slashpos
  character(len=200) :: basename

  ios = 0

  ! Prefer legacy naming when the basename starts with regcoil_in.
  slashpos = index(inputFilename, '/', back=.true.)
  if (slashpos == 0) then
     basename = trim(inputFilename)
  else
     basename = trim(inputFilename(slashpos+1:))
  end if
  if (basename(1:11) == "regcoil_in.") then
     output_filename = "regcoil_out" // trim(basename(11:)) // ".nc"
  else
     output_filename = "regcoil_out.nc"
  end if

  fileUnit=11
  open(unit=fileUnit, file=inputFilename, action="read", status="old", iostat=didFileAccessWork)
  if (didFileAccessWork /= 0) then
     print *,"Error opening input file ", trim(inputFilename)
     ios = 1
     return
  end if

  read(fileUnit, nml=regcoil_nml, iostat=didFileAccessWork)
  if (didFileAccessWork /= 0) then
     print *,"Error!  I was able to open the file ", trim(inputFilename), &
            " but not read data from the regcoil_nml namelist in it."
     backspace(fileUnit)
     read(fileUnit, fmt="(A)") linebuf
     print *, "Invalid line in namelist: ", trim(linebuf)
     print *, "(Error may be earlier)"
     if (didFileAccessWork==-1) then
        print *,"Make sure there is a carriage return after the / at the end of the namelist!"
     end if
     close(unit = fileUnit)
     ios = 2
     return
  end if
  if (verbose) print *,"Successfully read parameters from regcoil_nml namelist in ", trim(inputFilename), "."
  close(unit = fileUnit)

  if (verbose) then
     print *,"Resolution parameters:"
     print "(a,i5)","   ntheta_plasma  =",ntheta_plasma
     print "(a,i5)","   ntheta_coil    =",ntheta_coil
     print "(a,i5)","   nzeta_plasma   =",nzeta_plasma
     print "(a,i5)","   nzeta_coil     =",nzeta_coil
     print "(a,i5)","   mpol_potential =",mpol_potential
     print "(a,i5)","   ntor_potential =",ntor_potential

     select case (symmetry_option)
     case (1)
        print *,"Symmetry: sin(m*theta - n*zeta) modes only"
     case (2)
        print *,"Symmetry: cos(m*theta - n*zeta) modes only"
     case (3)
        print *,"Symmetry: both sin(m*theta - n*zeta) and cos(m*theta - n*zeta) modes"
     case default
        print *,"Error! Invalid setting for symmetry_option: ",symmetry_option
        ios = 3
        return
     end select
  end if

end subroutine regcoil_read_input_file
