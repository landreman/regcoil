module read_efit_mod

  use stel_kinds

  implicit none

  private

  public :: read_efit

  integer :: nwb, nhb
  real(dp) :: this_theta, psiN_desired, R_mag, Z_mag
  real(dp), dimension(:),   allocatable :: efit_R, efit_Z
  real(dp), dimension(:,:), allocatable :: efit_psi

contains

  subroutine read_efit(filename, psiN_desired_in, mnmax, rmnc, zmns, rmns, zmnc)
    
    ! Most of this subroutine is adapted from GS2's geometry module (geo/eeq.f90)
    
    ! Inputs: 
    !   filename: EFIT file to read
    !   psiN_desired
    use splines
    use stel_constants
    use stel_kinds
    
    implicit none
    
    character(*) :: filename
    integer, intent(in) :: mnmax
    real(dp), intent(in) :: psiN_desired_in
    real(dp), dimension(mnmax) :: rmnc, zmns, rmns, zmnc
    
    character(10) ::  char
    integer :: i, j, nw, nh, ierr
    integer :: big
    real(dp) :: xdum
    real(dp) :: rwid, zhei, rcentr, rleft, psi_0, psi_a, bcentr
    real(dp), dimension(:),   allocatable :: psi_bar, fp, qsf, pressure, beta, dummy, spsi_bar, sefit_R, sefit_Z
    real(dp), dimension(:,:), allocatable :: sefit_psi
    real(dp), dimension(:,:,:), allocatable :: dpm, dtm
    integer :: ndum, nbbbs  ! nbbbs = number of points for LCFS in efit file.
    real(dp), allocatable, dimension (:) :: rbbbs, zbbbs, thetab_atan, r_bound !boundary of plasma
    
    ! Variables for spline fit of psi(R,Z):
    real(dp), dimension(:), allocatable :: zp, temp, zx1, zxm, zy1, zyn
    real(dp) :: zxy11, zxym1, zxy1n, zxymn
    
    type (spline):: LCFS_spline
    integer :: ntheta = 20
    real(dp), dimension(:), allocatable :: theta_atan, theta_general, R_surface, Z_surface, r_LCFS
    
    ! Variables related to the root finding:
    integer :: fzeroFlag
    real(dp) :: rootSolve_abserr, rootSolve_relerr, rootSolve_min, rootSolve_max
    
    real(dp) :: r_minor
    real(dp) :: rtemp, ztemp, rtemps, ztempc, cosmu, sinmu, dnorm
    integer :: m
    real(dp) :: offset1, offset2
    real(dp), dimension(:), allocatable :: thetab_atan_3, r_bound_3
    
    psiN_desired = psiN_desired_in

    ! nw, nh = original number of grid points in width and height of efit file. ('small' grid)
    ! big = some factor by which the EFIT (R,Z) grid is refined.
    big = 20
    ! nwb, nhb = number of grid points in width and height for refined grid. ('big' grid)
    ! efit_psi = psi on 'big' grid
    ! sefit_psi = psi on 'small' (original) grid
    
    
    print *,"  Reading EFIT file ",trim(filename)
    
    
    if (psiN_desired > 1) then
       stop "Error! psiN_desired must be <= 1."
    end if
    
    if (psiN_desired <= 0) then
       stop "Error! psiN_desired must be more than 0."
    end if
    
    
    ! Initialize a theta grid.
    ! (Factor 10 here is arbitrary. Really it just needs to be large enough to resolve mnmax modes without aliasing.)
    ntheta = mnmax * 10
    allocate(theta_general(ntheta))
    allocate(theta_atan(ntheta))
    do i=1,ntheta
       theta_general(i) = i-1
    end do
    theta_general = theta_general*twopi/ntheta - pi
    
    offset1 = 2*pi*0.71
    offset2 = 2*pi*0.3
    theta_atan = theta_general - 0.6*sin((theta_general-offset1)/2)*cos((theta_general-offset1)/2) &
         / (0.3 + sin((theta_general-offset1)/2)**2) &
         - 0.32*sin((theta_general-offset2)/2)*cos((theta_general-offset2)/2) &
         / (0.2 + sin((theta_general-offset2)/2)**2)
    
    ! Beginning of code taken from gs2/geo/eeq.f90
    
    
    
    open(unit=5,file=trim(filename),status='old',form='formatted')
    
    ! Read the data
    
    read(5,1000) char, char, char, char, char, i, nw, nh
1000 format(5(a10),i2,i4,i4)
    
    !if (verbosity>1) write (*,*) 'EFIT dimensions are: ', nw, nh
    
    nwb = nw * big
    nhb = nh * big
    
    ! Next line replaced with the allocates by MJL
    !call alloc_module_arrays(nwb, nwb, nhb, nw, nh)
    !subroutine alloc_module_arrays(np, nw, nh, nws, nhs, ntime)
    allocate(psi_bar(nwb), fp(nwb), qsf(nwb), pressure(nwb), beta(nwb))
    allocate(dummy(nw), efit_R(nwb), efit_Z(nhb))
    allocate(spsi_bar(nw), sefit_R(nw), sefit_Z(nh))
    allocate(efit_psi(nwb, nhb), sefit_psi(nw, nh))
    !allocate(dpm(nwb, nhb, 2), dtm(nwb, nhb, 2))
    
    
    read(5,2020) rwid, zhei, rcentr, rleft, xdum      
    read(5,2020) R_mag, Z_mag, psi_0, psi_a, bcentr
2020 format (5e16.9)
    
    
    !
    ! pbar is defined by
    ! pbar = (psi-psi_0)/(psi_a-psi_0)
    ! fp and q are functions of pbar
    !
    do i=1,2
       read(5,2020)xdum,xdum,xdum,xdum,xdum
    enddo
    
    do i=1,nw
       spsi_bar(i) = float(i-1) / float(nw-1)
       sefit_R(i) = rleft +rwid * float(i-1) / float(nw-1)
    enddo
    
    do j=1,nh
       sefit_Z(j) = ((float(j-1) / float(nh-1) ) - 0.5)*zhei
    enddo
    
    do i=1,nwb
       psi_bar(i) = float(i-1) / float(nwb-1)
       efit_R(i) = rleft +rwid * float(i-1) / float(nwb-1)
    enddo
    
    do j=1,nhb
       efit_Z(j) = ((float(j-1) / float(nhb-1) ) - 0.5)*zhei
    enddo
    
    ! Read T = R * B_toroidal
    read(5,2020) (dummy(j),   j = 1, nw)
    !MJL call inter_cspl(nw, spsi_bar, dummy, nwb, psi_bar, fp)
    ! Read pressure
    read(5,2020) (dummy(j), j = 1, nw)
    !MJL call inter_cspl(nw, spsi_bar, dummy, nwb, psi_bar, pressure)
    ! Read T * (dT/dpsi)?
    read(5,2020) (dummy(j),     j = 1, nw)
    ! Read dp/dpsi?
    read(5,2020) (dummy(j),     j = 1 ,nw)
    ! Read psi(R,Z):
    read(5,2020) ((sefit_psi(i,j) , i = 1, nw) , j = 1, nh)
    
    ! Read q(psi)
    read(5,2020) (dummy(j) ,   j = 1, nw)
    
    ! Read # of points describing the LCFS:
    read(5,2022) nbbbs, ndum
2022 format (2i5)      
    
    allocate(rbbbs(nbbbs), zbbbs(nbbbs), thetab_atan(nbbbs), r_bound(nbbbs))
    ! Read in location of the LCFS:
    read(5,2020) (rbbbs(i), zbbbs(i) , i = 1, nbbbs)
    
    ! End of info in the EFIT file.
    close(5)
    
    
    ! get r_boundary(theta)
    thetab_atan = atan2 ((zbbbs-Z_mag), (rbbbs-R_mag))
    r_bound = sqrt( (rbbbs - R_mag)**2 + (zbbbs - Z_mag)**2 )
    
    call sort(thetab_atan, r_bound, zbbbs, rbbbs, nbbbs)
    
    
    ! Allow for duplicated points near +- pi:
    
    if(thetab_atan(1) == thetab_atan(2)) then
       thetab_atan(1) = thetab_atan(1) + twopi
       call sort(thetab_atan, r_bound, zbbbs, rbbbs, nbbbs)
    endif
    
    if(thetab_atan(nbbbs-1) == thetab_atan(nbbbs)) then
       thetab_atan(nbbbs) = thetab_atan(nbbbs) - twopi
       call sort(thetab_atan, r_bound, zbbbs, rbbbs, nbbbs)
    endif
    
    ! It isn't likely that a duplicate point would exist near theta = 0, 
    ! so I am not allowing this possibility for now.
    
    do i=1,nbbbs-1
       if(thetab_atan(i+1) == thetab_atan(i)) then
          !          write(*,*) 'Duplicates near theta = 0 not allowed.'
          !          write(*,*) i, i+1, ' Stopping.'
          !          stop
          !
          ! put in kluge for duplicate points, which happens near theta=0:
          thetab_atan(i+1) = thetab_atan(i+1)+1.e-8
       endif
    enddo
    
    
!!$  print *,"psi_0:",psi_0
!!$  print *,"psi_a:",psi_a
!!$  print *,"sefit_R:"
!!$  print *,sefit_R
!!$  print *,"sefit_Z:"
!!$  print *,sefit_Z
    
    ! Normalize and shift psi so it goes from 0 on axis to 1 at the LCFS.
    sefit_psi = 1 - (sefit_psi-psi_a)/(psi_0-psi_a)
    
!!$  print *,"max(sefit_psi):",maxval(sefit_psi)
!!$  print *,"min(sefit_psi):",minval(sefit_psi)
    
    
    
    ! Use splines to interpolate psi(R,Z) from the original coarse grid to the fine ('big') grid:
    allocate(zp(3*nw*nh), temp(nw+2*nh))
    allocate(zx1(nh), zxm(nh), zy1(nw), zyn(nw))
    
    call fitp_surf1(nw, nh, sefit_R, sefit_Z, sefit_psi, &
         nw, zx1, zxm, zy1, zyn, zxy11, zxym1, zxy1n, zxymn, &
         255, zp, temp, 1.0_dp, ierr)
    
    do j = 1, nhb
       do i = 1, nwb
          efit_psi(i,j) = fitp_surf2(efit_R(i), efit_Z(j), nw, nh, &
               sefit_R, sefit_Z, sefit_psi, nw, zp, 1.0_dp)
       enddo
    enddo
    
    deallocate(zp, temp)
    
!!$  print *,"max(efit_psi):",maxval(efit_psi)
!!$  print *,"min(efit_psi):",minval(efit_psi)
    
    
    ! Next, use spline interpolation to get r_LCFS on the theta grid we want.
    ! Since 'thetab_atan' does not run from exactly -pi to pi, there will be a bit of extrapolation
    ! at the ends of the grid. This could be eliminated if it is a problem.
    allocate(thetab_atan_3(nbbbs*3))
    allocate(r_bound_3(nbbbs*3))
    thetab_atan_3(1:nbbbs) = thetab_atan-2*pi
    thetab_atan_3((nbbbs+1):(2*nbbbs)) = thetab_atan
    thetab_atan_3((2*nbbbs+1):(3*nbbbs)) = thetab_atan+2*pi
    r_bound_3(1:nbbbs) = r_bound
    r_bound_3((nbbbs+1):(2*nbbbs)) = r_bound
    r_bound_3((2*nbbbs+1):(3*nbbbs)) = r_bound
    
    !call new_spline(nbbbs, thetab_atan, r_bound, LCFS_spline)
    call new_spline(nbbbs*3, thetab_atan_3, r_bound_3, LCFS_spline)
    allocate(r_LCFS(ntheta))
    do i = 1,ntheta
       r_LCFS(i) = splint(theta_atan(i), LCFS_spline)
    end do
    
!!$  print *," "
!!$  print *,"thetab_atan:",thetab_atan
!!$  print *," "
!!$  print *,"theta:",theta
!!$  print *," "
!!$  print *,"r_LCFS:",r_LCFS
!!$  print *," "

    ! Now solve the root-finding problem at each theta
    allocate(R_surface(ntheta), Z_surface(ntheta))
    
    rootSolve_abserr = 1.0e-10_dp
    rootSolve_relerr = 1.0e-10_dp
    do i = 1,ntheta
       print *,"i=",i,"of",ntheta
       rootSolve_min = 0.0_dp
       rootSolve_max = r_LCFS(i)
       this_theta = theta_atan(i)
       
       if (abs(psiN_desired-1) < 1d-7) then
          r_minor = r_LCFS(i)
          if (i==1) then
             print *,"Using EFIT LCFS directly"
          end if
       else
          if (i==1) then
             print *,"Computing a flux surface inside the EFIT LCFS."
          end if
          
          call fzero(fzero_residual, rootSolve_min, rootSolve_max, rootSolve_max * psiN_desired, &
               rootSolve_relerr, rootSolve_abserr, fzeroFlag)
          ! Note: fzero returns its answer in rootSolve_min
          r_minor = rootSolve_min
          if (fzeroFlag == 4) then
             stop "ERROR: fzero returned error 4: no sign change in residual"
          else if (fzeroFlag > 2) then
             print *,"WARNING: fzero returned an error code:",fzeroFlag
          end if
       end if
       
       R_surface(i) = Rpos(r_minor, this_theta)
       !Z_surface(i) = Zpos(r_minor, this_theta) - Z_mag  ! Shift vertically so magnetic axis is now at Z=0
       Z_surface(i) = Zpos(r_minor, this_theta) 
    end do
    
    print *,"Begin Fourier coefficients of plasma surface ----------"
    ! Begin code adapted from BNORM to Fourier transform data on the theta grid to sin/cos coefficients
    do m = 0, mnmax-1
       dnorm = (1.0_dp)/ntheta
       if (m.ne.0) dnorm = 2*dnorm
       rtemp = 0
       ztemp = 0
       rtemps = 0
       ztempc = 0
       do i = 1, ntheta
          cosmu = cos(m*theta_general(i))*dnorm
          sinmu = sin(m*theta_general(i))*dnorm
          rtemp  = rtemp  + R_surface(i) * cosmu
          ztemp  = ztemp  + Z_surface(i) * sinmu
          rtemps = rtemps + R_surface(i) * sinmu
          ztempc = ztempc + Z_surface(i) * cosmu
       end do
       rmnc(m+1) = rtemp
       zmns(m+1) = ztemp
       rmns(m+1) = rtemps
       zmnc(m+1) = ztempc
       ! Example of format for VMEC input:
       ! RBC(-12, 0) = -5.1637E-05 ZBS(-12, 0) = -9.4960E-05
       print "(a,i3.3,a,e20.12,a,i3.3,a,e20.12)","RBC(0,",m,")=",rmnc(m+1)," ZBS(0,",m,")=",zmns(m+1)
       print "(a,i3.3,a,e20.12,a,i3.3,a,e20.12)","RBS(0,",m,")=",rmns(m+1)," ZBC(0,",m,")=",zmnc(m+1)
    end do
    print *,"End Fourier coefficients of plasma surface ----------"
    
    print *,"Begin Fourier coefficients of plasma surface ----------"
    print *,'m, rmnc, rmns, zmnc, zmns'
    do m = 0, mnmax-1
       print "(i3.3,4e20.12)",m,rmnc(m+1),rmns(m+1),zmnc(m+1),zmns(m+1)
    end do
    print *,"End Fourier coefficients of plasma surface ----------"
    
!!$  print *,"r_surface:",r_surface
!!$  print *," "
    print *,"  Succeeded reading and processing EFIT file."
        
  end subroutine read_efit
    
! --------------------------------------------------------------------------------

  function fzero_residual(r)
    implicit none
    
    real(dp), intent (in) :: r
    real(dp) :: f, fzero_residual
    
    call eqitem(r, this_theta, efit_psi, f)
    fzero_residual = f - psiN_desired
  end function fzero_residual

! --------------------------------------------------------------------------------

  subroutine eqitem(r, thetin, f, fstar)
    ! Given a quantity f(:,:) on the 'big' (i.e. refined) (R,Z) grid,
    ! linearly interpolate to evaluate the quantity at a given (r,theta)
    ! where theta is the angle defined by tan(theta) = (Z-Z_axis)/(R-R_axis).
    
    implicit none
    
    real(dp), intent (in) :: r, thetin, f(:,:)
    real(dp), intent (out) :: fstar
    integer :: i, j, istar, jstar
    real(dp) :: st, dt, sr, dr
    real(dp) :: r_pos, z_pos
    
    r_pos = Rpos(r, thetin)
    z_pos = Zpos(r, thetin)
    
    ! find point on R mesh
    
    if(r_pos >= efit_R(nwb) .or. r_pos <= efit_R(1)) then
       write(*,*) 'No evaluation of eqitem allowed outside'
       write(*,*) 'or on edge of R domain'
       write(*,*) r, thetin, efit_R(1), efit_R(nwb), r_pos
       stop
    endif
    
    ! ensure point is on Z mesh
    
    if(z_pos >= efit_Z(nhb) .or. z_pos <= efit_Z(1)) then
       write(*,*) 'No evaluation of eqitem allowed outside'
       write(*,*) 'or on edge of Z domain'
       write(*,*) r, thetin, efit_Z(1), efit_Z(nhb), z_pos
       stop
    endif
    
    istar=0
    do i=2,nwb
       if(istar /= 0) cycle
       if(r_pos < efit_R(i)) then
          dr = r_pos - efit_R(i-1)
          sr = efit_R(i) - r_pos
          istar=i-1
       endif
    enddo
    
    ! Now do Z direction
    
    jstar=0
    do j=1,nhb
       if(jstar /= 0) cycle
       if(z_pos < efit_Z(j)) then
          dt = z_pos - efit_Z(j-1)
          st = efit_Z(j) - z_pos
          jstar=j-1
       endif
    enddo
    
    ! use opposite area stencil to interpolate
    
    fstar=f(istar    , jstar    ) * sr * st &
         +f(istar + 1, jstar    ) * dr * st &
         +f(istar    , jstar + 1) * sr * dt &
         +f(istar + 1, jstar + 1) * dr * dt
    fstar = fstar &
         /abs(efit_R(istar+1)-efit_R(istar)) &
         /(efit_Z(jstar+1)-efit_Z(jstar))
    
  end subroutine eqitem

! -------------------------------------------------------------------------------- 

 function Zpos (r, theta)
    implicit none
    real(dp), intent (in) :: r, theta
    real(dp) :: Zpos
    
    Zpos = Z_mag + r * sin(theta)
  end function Zpos
  
! --------------------------------------------------------------------------------

  function Rpos (r, theta)
    implicit none
    real(dp), intent (in) :: r, theta
    real(dp) :: Rpos
    
    Rpos = R_mag + r * cos(theta)
  end function Rpos
        
  
end module read_efit_mod

! --------------------------------------------------------------------------------
  
subroutine sort(a, b, c, d, N)
    
  use stel_kinds
  
  implicit none
  
  integer, intent(in) :: N
  real(dp), dimension(N), intent(in out) :: a, b, c, d
  real(dp) :: tmp
  integer :: i, jmax
  logical :: sorted
  
  jmax = size(a)
  
  do 
     sorted = .true.
     do i=1,jmax-1
        if(a(i+1) < a(i)) then
           tmp=a(i); a(i)=a(i+1); a(i+1)=tmp             
           tmp=b(i); b(i)=b(i+1); b(i+1)=tmp
           tmp=c(i); c(i)=c(i+1); c(i+1)=tmp
           tmp=d(i); d(i)=d(i+1); d(i+1)=tmp
           sorted = .false.
        endif
     enddo
     if (sorted) exit
  enddo
  
end subroutine sort
