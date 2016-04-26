subroutine compute_Kmn

  ! Notation in this subroutine refers to the appendix of Merkel 1987 Varenna paper.
  ! Although Merkel J Comp Phys 66, 83 (1986) appears similar, it is filled with mathematical errors!
  ! These errors appear to be absent in the Varenna version.

  use global_variables, only: nu_middle, nv_middle, drdu_middle, drdv_middle, d2rdu2_middle, d2rdudv_middle, d2rdv2_middle, normal_middle, &
       Merkel_Kmn, xm_outer, xn_outer, mnmax_outer, mpol_outer, ntor_outer
  use stel_kinds

  implicit none

  integer :: iu, iv, index, factorial_max, l, max_index, m, n, lmax, imn
  integer :: tic, toc, countrate, iflag
  real(dp) :: a, b, c, AA, BB, CC
  real(dp) :: aPlus, aMinus, d
  real(dp) :: AAPlus, AAMinus, DD
  real(dp) :: k1Plus, k1Minus, k2Plus, k2Minus
  real(dp), dimension(:), allocatable :: log_factorial
  real(dp), dimension(:), allocatable :: TPlus, TMinus, SPlus, SMinus
  real(dp), dimension(:,:), allocatable :: cSPlus, cSMinus
  real(dp), dimension(:,:,:), allocatable :: factor

  call system_clock(tic,countrate)
  print *,"Beginning computation of Kmn."

  ! Initialize array of n!
  factorial_max = mpol_outer+ntor_outer
  !allocate(factorial(0:factorial_max),stat=iflag)
  allocate(log_factorial(factorial_max+1),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  log_factorial = 0
  do index = 2,factorial_max
     log_factorial(index+1) = log(1.0_dp*index) + log_factorial(index)
  end do

  max_index = max(mpol_outer+ntor_outer,1)
  allocate(TPlus(max_index+1),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(TMinus(max_index+1),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(SPlus(max_index+1),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(SMinus(max_index+1),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(cSPlus(mpol_outer+1,ntor_outer+1),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(cSMinus(mpol_outer+1,ntor_outer+1),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(factor(max(mpol_outer,ntor_outer)+1,mpol_outer+1,ntor_outer+1),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  factor = 0
  do m = 0,mpol_outer
     do n = 0,ntor_outer              
        lmax = (m+n-abs(m-n))/2
        if (lmax+1 > max(mpol_outer,ntor_outer)+1) then 
           stop "Out of bounds!"
        end if
        do l = 0,lmax
           factor(l+1,m+1,n+1) = 0.5*((-1)**(l+(abs(m-n)-m+n)/2)) * exp( log_factorial(1+(m+n+abs(m-n))/2 + l) &
                - (log_factorial(1+(m+n-abs(m-n))/2-l) + log_factorial(1+abs(m-n)+l) + log_factorial(1+l)))

           !factor(l+1,m+1,n+1) = 0.5*((-1)**(l+(abs(m-n)-m+n)/2)) * factorial(1+(m+n+abs(m-n))/2 + l) &
           !     / (factorial(1+(m+n-abs(m-n))/2-l) * factorial(1+abs(m-n)+l) * factorial(1+l))

           !factor(l+1,m+1,n+1) = 0.5*((-1)**(l+(abs(m-n)-m+n)/2)) * factorial((m+n+abs(m-n))/2 + l) &
           !     / (factorial((m+n-abs(m-n))/2-l) * factorial(abs(m-n)+l) * factorial(l))
        end do
     end do
  end do

  do iu = 1,nu_middle
     do iv = 1,nv_middle
        ! Since Fortran is case-insensitive, Merkel's uppercase characters are doubled, e.g. AA

        a = drdu_middle(1,iu,iv)*drdu_middle(1,iu,iv) + drdu_middle(2,iu,iv)*drdu_middle(2,iu,iv) + drdu_middle(3,iu,iv)*drdu_middle(3,iu,iv)
        b = drdu_middle(1,iu,iv)*drdv_middle(1,iu,iv) + drdu_middle(2,iu,iv)*drdv_middle(2,iu,iv) + drdu_middle(3,iu,iv)*drdv_middle(3,iu,iv)
        c = drdv_middle(1,iu,iv)*drdv_middle(1,iu,iv) + drdv_middle(2,iu,iv)*drdv_middle(2,iu,iv) + drdv_middle(3,iu,iv)*drdv_middle(3,iu,iv)
        
        AA = 0.5*( d2rdu2_middle(1,iu,iv)*normal_middle(1,iu,iv) &
             + d2rdu2_middle(2,iu,iv)*normal_middle(2,iu,iv) &
             + d2rdu2_middle(3,iu,iv)*normal_middle(3,iu,iv) )
        
        BB = 0.5*( d2rdudv_middle(1,iu,iv)*normal_middle(1,iu,iv) &
             + d2rdudv_middle(2,iu,iv)*normal_middle(2,iu,iv) &
             + d2rdudv_middle(3,iu,iv)*normal_middle(3,iu,iv) )
        
        CC = 0.5*( d2rdv2_middle(1,iu,iv)*normal_middle(1,iu,iv) &
             + d2rdv2_middle(2,iu,iv)*normal_middle(2,iu,iv) &
             + d2rdv2_middle(3,iu,iv)*normal_middle(3,iu,iv) )

        aPlus  = a + 2*b + c;
        aMinus = a - 2*b + c;
        d = c - a;
        AAPlus  = AA + 2*BB + CC;
        AAMinus = AA - 2*BB + CC;
        DD = CC-AA;
        k1Plus  = (a+c)*  BB  - (AA+CC)*  b;
        k1Minus = (a+c)*(-BB) - (AA+CC)*(-b);
        k2Plus  = CC*(a+b) -   BB *(c-a) - AA*(c+b);
        k2Minus = CC*(a-b) - (-BB)*(c-a) - AA*(c-b);

        TPlus(1)  = 1/sqrt(a+2*b+c) * log((sqrt(c*(a+2*b+c)) + c + b) / (sqrt(a*(a+2*b+c)) - a - b));
        TMinus(1) = 1/sqrt(a-2*b+c) * log((sqrt(c*(a-2*b+c)) + c - b) / (sqrt(a*(a-2*b+c)) - a + b));
        TPlus(2)  = 1/(a+2*b+c) * (2*(sqrt(c)-sqrt(a)) - (c-a)*TPlus(1));
        TMinus(2) = 1/(a-2*b+c) * (2*(sqrt(c)-sqrt(a)) - (c-a)*TMinus(1));
                    
        do l = 2,max_index
           TPlus(l+1)  = 1/(l*(a+2*b+c))*(2*(sqrt(c) + ((-1)**l)*sqrt(a)) - (2*l-1)*(c-a)*TPlus(l-1+1)  - (l-1)*(a-2*b+c)*TPlus(l-2+1));
           TMinus(l+1) = 1/(l*(a-2*b+c))*(2*(sqrt(c) + ((-1)**l)*sqrt(a)) - (2*l-1)*(c-a)*TMinus(l-1+1) - (l-1)*(a+2*b+c)*TMinus(l-2+1));
        end do
                    
        SPlus(1) = 1/((aPlus ** (1.5_dp)) * (aPlus*aMinus-d*d)) * (0 &
             + 1/sqrt(aMinus+aPlus-2*d)*sqrt(aPlus)*(AAMinus*aPlus*(aPlus-d)-aMinus*(AAPlus*d+aPlus*AAPlus-2*aPlus*DD)+2*d*(AAPlus*d-aPlus*DD)) &
             + 1/sqrt(aMinus+aPlus+2*d)*sqrt(aPlus)*(AAMinus*aPlus*(aPlus+d)+aMinus*(AAPlus*d-aPlus*AAPlus-2*aPlus*DD)+2*d*(AAPlus*d-aPlus*DD)) &
             + AAPlus*(d*d-aPlus*aMinus) * log( (d-aPlus+sqrt(aPlus*(aMinus+aPlus-2*d))) / (d+aPlus+sqrt(aPlus*(aMinus+aPlus+2*d))) ))
                    
        SMinus(1) = 1/((aMinus ** (1.5_dp)) * (aMinus*aPlus-d*d)) * (0 &
             + 1/sqrt(aPlus+aMinus-2*d)*sqrt(aMinus)*(AAPlus*aMinus*(aMinus-d)-aPlus*(AAMinus*d+aMinus*AAMinus-2*aMinus*DD)+2*d*(AAMinus*d-aMinus*DD)) &
             + 1/sqrt(aPlus+aMinus+2*d)*sqrt(aMinus)*(AAPlus*aMinus*(aMinus+d)+aPlus*(AAMinus*d-aMinus*AAMinus-2*aMinus*DD)+2*d*(AAMinus*d-aMinus*DD)) &
             + AAMinus*(d*d-aMinus*aPlus) * log( (d-aMinus+sqrt(aMinus*(aPlus+aMinus-2*d))) / (d+aMinus+sqrt(aMinus*(aPlus+aMinus+2*d))) ))
                    
        do l = 1,max_index
           SPlus(l+1) = ( ((AA+2*BB+CC)*(a*c-b*b) + l*(k1Plus*(a+2*b+c)+k2Plus*(c-a)))*TPlus(l+1) &
                +l*(k1Plus*(c-a)+k2Plus*(a-2*b+c))*TPlus(l-1+1) &
                -((c+b)*k1Plus+(c-b)*k2Plus)/sqrt(c) - ((-1)**l)*((a+b)*k1Plus-(a-b)*k2Plus)/sqrt(a)) &
                /((a+2*b+c)*(a*c-b*b))
                        
           SMinus(l+1) = ( ((AA-2*BB+CC)*(a*c-b*b) + l*(k1Minus*(a-2*b+c)+k2Minus*(c-a)))*TMinus(l+1) &
                +l*(k1Minus*(c-a)+k2Minus*(a+2*b+c))*TMinus(l-1+1) &
                -((c-b)*k1Minus+(c+b)*k2Minus)/sqrt(c) - ((-1)**l)*((a-b)*k1Minus-(a+b)*k2Minus)/sqrt(a)) &
                /((a-2*b+c)*(a*c-b*b))
        end do

        cSPlus  = 0
        cSMinus = 0
        do m = 0,mpol_outer
           do n = 0,ntor_outer              
              lmax = (m+n-abs(m-n))/2
              !do l = 0,lmax
              do l = lmax,0,-1
                 index = abs(m-n)+2*l+1
                 cSPlus( m+1,n+1) = cSPlus( m+1,n+1) + factor(l+1,m+1,n+1) *  SPlus(index)
                 cSMinus(m+1,n+1) = cSMinus(m+1,n+1) + factor(l+1,m+1,n+1) * SMinus(index)
              end do
           end do
        end do
                    
        index = (iv-1)*nu_middle + iu
        do imn = 1,mnmax_outer
           m = xm_outer(imn);
           n = xn_outer(imn);
           if (n<0) then
              if (m .ne. 0) then
                 ! m>=1, n<0, so flip cPlus and cMinus
                 Merkel_Kmn(index,imn) = cSMinus(m+1,-n+1) + cSMinus(m-1+1,-n+1) + cSMinus(m+1,-n-1+1) + cSMinus(m-1+1,-n-1+1);
              else
                 ! m=0, n<0, so flip cPlus and cMinus
                 Merkel_Kmn(index,imn) = cSMinus(m+1,-n+1) + cSMinus(m+1,-n-1+1) + cSPlus(m+1,-n+1) + cSPlus(m+1,-n-1+1);
              end if
           else
              ! n >= 0
              if (m .ne. 0) then 
                 if (n .ne. 0) then
                    ! m>=1, n>=1
                    Merkel_Kmn(index,imn) = cSPlus(m+1,n+1) + cSPlus(m-1+1,n+1) + cSPlus(m+1,n-1+1) + cSPlus(m-1+1,n-1+1);
                 else
                    ! m>=1, n=0
                    Merkel_Kmn(index,imn) = cSPlus(m+1,n+1) + cSPlus(m-1+1,n+1) + cSMinus(m+1,n+1) + cSMinus(m-1+1,n+1);
                 end if
              else
                 if (n .ne. 0) then
                    ! m=0, n>=1
                    Merkel_Kmn(index,imn) = cSPlus(m+1,n+1) + cSPlus(m+1,n-1+1) + cSMinus(m+1,n+1) + cSMinus(m+1,n-1+1);
                 else
                    ! m=0, n=0
                    Merkel_Kmn(index,imn) = cSPlus(m+1,n+1) + cSPlus(m+1,n+1) + cSMinus(m+1,n+1) + cSMinus(m+1,n+1);
                 end if
              end if
           end if
        end do
     end do
  end do

  deallocate(TPlus,TMinus,SPlus,SMinus,cSPlus,cSMinus,log_factorial,factor)

  call system_clock(toc)
  print *,"Done computing Kmn. Took ",real(toc-tic)/countrate," sec."

end subroutine compute_Kmn


