subroutine read_nescin(nescin_filename, r, drdu, drdv, d2rdu2, d2rdudv, d2rdv2, nu, nvl, u, vl, compute_2nd_derivs)

  use global_variables, only: nfp, geometry_option_outer, xm, xn, mnmax, rmnc_global => rmnc, zmns_global => zmns, rmns_global => rmns, zmnc_global => zmnc
  use safe_open_mod
  use stel_constants
  use stel_kinds
  
  implicit none
  
  character(*) :: nescin_filename
  integer, intent(in) :: nu, nvl
  real(dp), dimension(3,nu,nvl) :: r, drdu, drdv
  real(dp), dimension(3,nu,nvl) :: d2rdu2, d2rdudv, d2rdv2
  real(dp), dimension(nu)  :: u
  real(dp), dimension(nvl) :: vl
  logical :: compute_2nd_derivs
  
  integer :: iunit = 7, iu, iv, iflag
  integer :: m, n, ntotal, k, mr, nr, istat
  real(dp) :: rmnc, zmns, rmns, zmnc
  real(dp) :: angle, sinangle, cosangle, dsinangledu, dsinangledv, dcosangledu, dcosangledv
  real(dp) :: angle2, sinangle2, cosangle2, dsinangle2dv, dcosangle2dv
  real(dp) :: d2sinangle2dv2, d2cosangle2dv2
  real(dp) :: d2sinangledu2, d2sinangledudv, d2sinangledv2
  real(dp) :: d2cosangledu2, d2cosangledudv, d2cosangledv2

  character(300) :: myline
  character(*), parameter :: matchString = "------ Current Surface"

  r=0
  drdu=0
  drdv=0

  if (compute_2nd_derivs) then
     d2rdu2 = 0
     d2rdudv = 0
     d2rdv2 = 0
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


  if (geometry_option_outer==4) then
     ! Clear arrays associated with the plasma surface for offsetting,
     ! and replace them with the nescin values.
     deallocate(xm,xn,rmnc_global,zmns_global)
     if (allocated(rmns_global)) then
        deallocate(rmns_global)
     end if
     if (allocated(zmnc_global)) then
        deallocate(zmnc_global)
     end if
     mnmax = ntotal
     allocate(xm(mnmax),stat=iflag)
     if (iflag .ne. 0) stop "Allocation error! read_nescin 1"
     allocate(xn(mnmax),stat=iflag)
     if (iflag .ne. 0) stop "Allocation error! read_nescin 2"
     allocate(rmnc_global(mnmax),stat=iflag)
     if (iflag .ne. 0) stop "Allocation error! read_nescin 3"
     allocate(zmns_global(mnmax),stat=iflag)
     if (iflag .ne. 0) stop "Allocation error! read_nescin 4"
     allocate(rmns_global(mnmax),stat=iflag)
     if (iflag .ne. 0) stop "Allocation error! read_nescin 5"
     allocate(zmnc_global(mnmax),stat=iflag)
     if (iflag .ne. 0) stop "Allocation error! read_nescin 6"
  end if

  read (iunit, *)
  read (iunit, *)
  do k = 1, ntotal
     read (iunit, *) m, n, rmnc, zmns, rmns, zmnc

     if (geometry_option_outer==4) then
        ! Set arrays associated with offsetting surfaces
        xm(k) = m
        xn(k) = -n*nfp
        rmnc_global(k) = rmnc
        zmns_global(k) = zmns
        rmns_global(k) = rmns
        zmnc_global(k) = zmnc
     end if

     do iv = 1,nvl
        angle2 = twopi*vl(iv)/nfp
        sinangle2 = sin(angle2)
        cosangle2 = cos(angle2)
        dsinangle2dv = cosangle2*twopi/nfp
        dcosangle2dv = -sinangle2*twopi/nfp
        d2sinangle2dv2 = -twopi*twopi/(nfp*nfp)*sinangle2
        d2cosangle2dv2 = -twopi*twopi/(nfp*nfp)*cosangle2
        do iu = 1,nu
           angle = twopi*(m*u(iu) + n*vl(iv))
           sinangle = sin(angle)
           cosangle = cos(angle)
           dsinangledu = cosangle*twopi*m
           dcosangledu = -sinangle*twopi*m
           dsinangledv = cosangle*twopi*n
           dcosangledv = -sinangle*twopi*n
           d2sinangledu2  = -twopi*twopi*m*m*sinangle
           d2sinangledudv = -twopi*twopi*m*n*sinangle
           d2sinangledv2  = -twopi*twopi*n*n*sinangle
           d2cosangledu2  = -twopi*twopi*m*m*cosangle
           d2cosangledudv = -twopi*twopi*m*n*cosangle
           d2cosangledv2  = -twopi*twopi*n*n*cosangle
           
           r(1,iu,iv) = r(1,iu,iv) + rmnc * cosangle * cosangle2 + rmns * sinangle * cosangle2
           r(2,iu,iv) = r(2,iu,iv) + rmnc * cosangle * sinangle2 + rmns * sinangle * sinangle2
           r(3,iu,iv) = r(3,iu,iv) + zmns * sinangle             + zmnc * cosangle
           
           drdu(1,iu,iv) = drdu(1,iu,iv) + rmnc * dcosangledu * cosangle2 + rmns * dsinangledu * cosangle2
           drdu(2,iu,iv) = drdu(2,iu,iv) + rmnc * dcosangledu * sinangle2 + rmns * dsinangledu * sinangle2
           drdu(3,iu,iv) = drdu(3,iu,iv) + zmns * dsinangledu + zmnc * dcosangledu
           
           drdv(1,iu,iv) = drdv(1,iu,iv) + rmnc * (dcosangledv * cosangle2 + cosangle * dcosangle2dv) &
                + rmns * (dsinangledv * cosangle2 + sinangle * dcosangle2dv)
           drdv(2,iu,iv) = drdv(2,iu,iv) + rmnc * (dcosangledv * sinangle2 + cosangle * dsinangle2dv) &
                + rmns * (dsinangledv * sinangle2 + sinangle * dsinangle2dv)
           drdv(3,iu,iv) = drdv(3,iu,iv) + zmns * dsinangledv + zmnc * dcosangledv

           if (compute_2nd_derivs) then
              d2rdu2(1,iu,iv) = d2rdu2(1,iu,iv) + rmnc * d2cosangledu2 * cosangle2 + rmns * d2sinangledu2 * cosangle2
              d2rdu2(2,iu,iv) = d2rdu2(2,iu,iv) + rmnc * d2cosangledu2 * sinangle2 + rmns * d2sinangledu2 * sinangle2
              d2rdu2(3,iu,iv) = d2rdu2(3,iu,iv) + zmns * d2sinangledu2 + zmnc * d2cosangledu2

              d2rdudv(1,iu,iv) = d2rdudv(1,iu,iv) + rmnc * (d2cosangledudv * cosangle2 + dcosangledu * dcosangle2dv) &
                   + rmns * (d2sinangledudv * cosangle2 + dsinangledu * dcosangle2dv)
              d2rdudv(2,iu,iv) = d2rdudv(2,iu,iv) + rmnc * (d2cosangledudv * sinangle2 + dcosangledu * dsinangle2dv) &
                   + rmns * (d2sinangledudv * sinangle2 + dsinangledu * dsinangle2dv)
              d2rdudv(3,iu,iv) = d2rdudv(3,iu,iv) + zmns * d2sinangledudv + zmnc * d2cosangledudv
           
              d2rdv2(1,iu,iv) = d2rdv2(1,iu,iv) + rmnc * (d2cosangledv2 * cosangle2 + dcosangledv * dcosangle2dv &
                   + dcosangledv * dcosangle2dv + cosangle * d2cosangle2dv2) &
                   + rmns * (d2sinangledv2 * cosangle2 + dsinangledv * dcosangle2dv &
                   + dsinangledv * dcosangle2dv + sinangle * d2cosangle2dv2)
              d2rdv2(2,iu,iv) = d2rdv2(2,iu,iv) + rmnc * (d2cosangledv2 * sinangle2 + dcosangledv * dsinangle2dv &
                   + dcosangledv * dsinangle2dv + cosangle * d2sinangle2dv2) &
                   + rmns * (d2sinangledv2 * sinangle2 + dsinangledv * dsinangle2dv &
                   + dsinangledv * dsinangle2dv + sinangle * d2sinangle2dv2)
              d2rdv2(3,iu,iv) = d2rdv2(3,iu,iv) + zmns * d2sinangledv2 + zmnc * d2cosangledv2
           end if
        end do
     end do
  end do
  
  
  close(iunit)

end subroutine  read_nescin
