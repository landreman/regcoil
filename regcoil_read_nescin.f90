subroutine regcoil_read_nescin(nescin_filename, r, drdtheta, drdzeta, d2rdtheta2, d2rdthetadzeta, d2rdzeta2, ntheta, nzetal, theta, zetal, compute_2nd_derivs)

  use regcoil_variables, only: nfp, xm, xn, mnmax, rmnc_global => rmnc, zmns_global => zmns, rmns_global => rmns, zmnc_global => zmnc
  use safe_open_mod
  use stel_constants
  use stel_kinds
  
  implicit none
  
  character(*) :: nescin_filename
  integer, intent(in) :: ntheta, nzetal
  real(dp), dimension(3,ntheta,nzetal) :: r, drdtheta, drdzeta
  real(dp), dimension(3,ntheta,nzetal) :: d2rdtheta2, d2rdthetadzeta, d2rdzeta2
  real(dp), dimension(ntheta)  :: theta
  real(dp), dimension(nzetal) :: zetal
  logical :: compute_2nd_derivs
  
  integer :: iunit = 7, itheta, izeta, iflag
  integer :: m, n, ntotal, k, mr, nr, istat
  real(dp) :: rmnc, zmns, rmns, zmnc
  real(dp) :: angle, sinangle, cosangle, dsinangledtheta, dsinangledzeta, dcosangledtheta, dcosangledzeta
  real(dp) :: angle2, sinangle2, cosangle2, dsinangle2dzeta, dcosangle2dzeta
  real(dp) :: d2sinangle2dzeta2, d2cosangle2dzeta2
  real(dp) :: d2sinangledtheta2, d2sinangledthetadzeta, d2sinangledzeta2
  real(dp) :: d2cosangledtheta2, d2cosangledthetadzeta, d2cosangledzeta2

  character(300) :: myline
  character(*), parameter :: matchString = "------ Current Surface"

  r=0
  drdtheta=0
  drdzeta=0

  if (compute_2nd_derivs) then
     d2rdtheta2 = 0
     d2rdthetadzeta = 0
     d2rdzeta2 = 0
  end if

  call safe_open(iunit, istat, trim(nescin_filename), 'old', 'formatted')
  if (istat .ne. 0) then
     stop 'Error opening nescin file'
  endif


  ! Skip down to the line in the file that begins with matchString
  do
     read (iunit,"(a)") myline
     if (myline(:len(matchString)) == matchString) then
        exit
     end if
  end do
  read (iunit, *)

  read (iunit, *) ntotal
  print *,"  Reading",ntotal,"Fourier modes from nescin"


  read (iunit, *)
  read (iunit, *)
  do k = 1, ntotal
     read (iunit, *) m, n, rmnc, zmns, rmns, zmnc

     do izeta = 1,nzetal
        angle2 = zetal(izeta)
        sinangle2 = sin(angle2)
        cosangle2 = cos(angle2)
        dsinangle2dzeta = cosangle2
        dcosangle2dzeta = -sinangle2
        d2sinangle2dzeta2 = -sinangle2
        d2cosangle2dzeta2 = -cosangle2
        do itheta = 1,ntheta
           angle = m*theta(itheta) + n*nfp*zetal(izeta)
           sinangle = sin(angle)
           cosangle = cos(angle)
           dsinangledtheta = cosangle*m
           dcosangledtheta = -sinangle*m
           dsinangledzeta = cosangle*n*nfp
           dcosangledzeta = -sinangle*n*nfp
           d2sinangledtheta2  = -m*m*sinangle
           d2sinangledthetadzeta = -m*n*nfp*sinangle
           d2sinangledzeta2  = -n*n*nfp*nfp*sinangle
           d2cosangledtheta2  = -m*m*cosangle
           d2cosangledthetadzeta = -m*n*nfp*cosangle
           d2cosangledzeta2  = -n*n*nfp*nfp*cosangle
           
           r(1,itheta,izeta) = r(1,itheta,izeta) + rmnc * cosangle * cosangle2 + rmns * sinangle * cosangle2
           r(2,itheta,izeta) = r(2,itheta,izeta) + rmnc * cosangle * sinangle2 + rmns * sinangle * sinangle2
           r(3,itheta,izeta) = r(3,itheta,izeta) + zmns * sinangle             + zmnc * cosangle
           
           drdtheta(1,itheta,izeta) = drdtheta(1,itheta,izeta) + rmnc * dcosangledtheta * cosangle2 + rmns * dsinangledtheta * cosangle2
           drdtheta(2,itheta,izeta) = drdtheta(2,itheta,izeta) + rmnc * dcosangledtheta * sinangle2 + rmns * dsinangledtheta * sinangle2
           drdtheta(3,itheta,izeta) = drdtheta(3,itheta,izeta) + zmns * dsinangledtheta + zmnc * dcosangledtheta
           
           drdzeta(1,itheta,izeta) = drdzeta(1,itheta,izeta) + rmnc * (dcosangledzeta * cosangle2 + cosangle * dcosangle2dzeta) &
                + rmns * (dsinangledzeta * cosangle2 + sinangle * dcosangle2dzeta)
           drdzeta(2,itheta,izeta) = drdzeta(2,itheta,izeta) + rmnc * (dcosangledzeta * sinangle2 + cosangle * dsinangle2dzeta) &
                + rmns * (dsinangledzeta * sinangle2 + sinangle * dsinangle2dzeta)
           drdzeta(3,itheta,izeta) = drdzeta(3,itheta,izeta) + zmns * dsinangledzeta + zmnc * dcosangledzeta

           if (compute_2nd_derivs) then
              d2rdtheta2(1,itheta,izeta) = d2rdtheta2(1,itheta,izeta) + rmnc * d2cosangledtheta2 * cosangle2 + rmns * d2sinangledtheta2 * cosangle2
              d2rdtheta2(2,itheta,izeta) = d2rdtheta2(2,itheta,izeta) + rmnc * d2cosangledtheta2 * sinangle2 + rmns * d2sinangledtheta2 * sinangle2
              d2rdtheta2(3,itheta,izeta) = d2rdtheta2(3,itheta,izeta) + zmns * d2sinangledtheta2 + zmnc * d2cosangledtheta2

              d2rdthetadzeta(1,itheta,izeta) = d2rdthetadzeta(1,itheta,izeta) + rmnc * (d2cosangledthetadzeta * cosangle2 + dcosangledtheta * dcosangle2dzeta) &
                   + rmns * (d2sinangledthetadzeta * cosangle2 + dsinangledtheta * dcosangle2dzeta)
              d2rdthetadzeta(2,itheta,izeta) = d2rdthetadzeta(2,itheta,izeta) + rmnc * (d2cosangledthetadzeta * sinangle2 + dcosangledtheta * dsinangle2dzeta) &
                   + rmns * (d2sinangledthetadzeta * sinangle2 + dsinangledtheta * dsinangle2dzeta)
              d2rdthetadzeta(3,itheta,izeta) = d2rdthetadzeta(3,itheta,izeta) + zmns * d2sinangledthetadzeta + zmnc * d2cosangledthetadzeta
           
              d2rdzeta2(1,itheta,izeta) = d2rdzeta2(1,itheta,izeta) + rmnc * (d2cosangledzeta2 * cosangle2 + dcosangledzeta * dcosangle2dzeta &
                   + dcosangledzeta * dcosangle2dzeta + cosangle * d2cosangle2dzeta2) &
                   + rmns * (d2sinangledzeta2 * cosangle2 + dsinangledzeta * dcosangle2dzeta &
                   + dsinangledzeta * dcosangle2dzeta + sinangle * d2cosangle2dzeta2)
              d2rdzeta2(2,itheta,izeta) = d2rdzeta2(2,itheta,izeta) + rmnc * (d2cosangledzeta2 * sinangle2 + dcosangledzeta * dsinangle2dzeta &
                   + dcosangledzeta * dsinangle2dzeta + cosangle * d2sinangle2dzeta2) &
                   + rmns * (d2sinangledzeta2 * sinangle2 + dsinangledzeta * dsinangle2dzeta &
                   + dsinangledzeta * dsinangle2dzeta + sinangle * d2sinangle2dzeta2)
              d2rdzeta2(3,itheta,izeta) = d2rdzeta2(3,itheta,izeta) + zmns * d2sinangledzeta2 + zmnc * d2cosangledzeta2
           end if
        end do
     end do
  end do
  
  
  close(iunit)

end subroutine  regcoil_read_nescin
