subroutine regcoil_read_input

  use regcoil_variables

  implicit none

  integer :: numargs
  character(len=200) :: inputFilename
  integer :: fileUnit, didFileAccessWork, i
  integer, parameter :: uninitialized = -9999

  ! getcarg is in LIBSTELL
  call getcarg(1, inputFilename, numargs)

  if (numargs<1) then
     stop "One argument is required: the input namelist file, which must be named regcoil_in.XXXXX"
  end if
  if (numargs>1) then
     print *,"WARNING: Arguments after the first will be ignored."
  end if
  if (inputFilename(1:11) .ne. "regcoil_in.") then
     stop "Input file must be named regcoil_in.XXX for some extension XXX"
  end if

  output_filename = "regcoil_out" // trim(inputFilename(11:)) // ".nc"

  fileUnit=11
  open(unit=fileUnit, file=inputFilename, action="read", status="old", iostat=didFileAccessWork)
  if (didFileAccessWork /= 0) then
     print *,"Error opening input file ", trim(inputFilename)
     stop
  else
     read(fileUnit, nml=regcoil_nml, iostat=didFileAccessWork)
     if (didFileAccessWork /= 0) then
        print *,"Error!  I was able to open the file ", trim(inputFilename), &
               " but not read data from the regcoil_nml namelist in it."
        if (didFileAccessWork==-1) then
           print *,"Make sure there is a carriage return after the / at the end of the namelist!"
        end if
        stop
     end if
     if (verbose) print *,"Successfully read parameters from regcoil_nml namelist in ", trim(inputFilename), "."
  end if
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
        stop
     end select
  end if

end subroutine regcoil_read_input
