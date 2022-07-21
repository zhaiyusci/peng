!***********************************************************************
! The main program is used here show how to call the potential energy 
! function.  The users should modify it.
!***********************************************************************
! program main
  ! implicit none
  ! real(8)  :: R, V, dVdR, d2VdR2

  ! R = 1.d0
  ! write(6,888) "R(A)", "V(cm**-1)", "dVdR(cm**-1/A)", "d2VdR2(cm**-1/A**2)"
  ! 888 format(4a25)
  ! do while (R .le. 10.d0)
    ! call V_MLR_Ar_Ar(R, V, dVdR, d2VdR2)
    ! write(6,999) R, V, dVdR, d2VdR2
    ! R = R+0.5d0
  ! enddo
  ! 999 format(4f25.4)
! end program


!***********************************************************************
subroutine V_MLR_Ar_Ar(R, V, dVdR, d2VdR2)
  ! 
  !  Morse/Long-Range potential energy curve for Ar-Ar system.
  !  Raw data from:
  !    Myatt, Dham, Chandrasekhar, McCourt, and Le Roy
  !    Mol. Phys. 116, 1598 (2018)
  !
  !* On input:
  !    real(8) R:       (Angstrom)
  !      nuclear distance where the potential energy is calculated.
  !* On output:
  !    real(8) V:       (cm**-1)
  !      potential energy calculated at R (in cm-1).
  !    real(8) dVdR:    (cm**-1/Angstrom)
  !      first-order and
  !    real(8) d2VdR2V: (cm**-1/Angstrom**2)
  !      second-order derivatives of V at R.
  !
  implicit none
  real(8), intent(in)  :: R
  real(8), intent(out) :: V, dVdR, d2VdR2
  !*********************************************************************
  ! Parameters...
  real(8), parameter   :: De = 99.49, Re = 3.766, Rref = 5.3
  integer, parameter   :: p = 6 , q = 3
  integer, parameter   :: N_C = 3 
  integer, parameter   :: m(N_C) = (/ 6, 8, 10 /)
  real(8), parameter   :: C(N_C) = &
    (/ 3.105D+05, 2.219D+06, 1.899D+07 /)
  integer, parameter   :: beta_max = 4
  real(8), parameter   :: beta(0:beta_max) = &
    (/0.04743, 0.0877, 0.112, 0.57, 0.33/)
  real(8), parameter   :: rho = 1.1, s=-0.5
  !*********************************************************************
  real(8)              :: beta_inf = 0.0d0
  real(8)              :: A, dAdR, d2AdR2
  real(8)              :: uLR, duLRdR, d2uLRdR2
  real(8)              :: uLRe
  real(8)              :: EE, dEEdR, d2EEdR2
  real(8)              :: E, dEdR, d2EdR2
  real(8)              :: beta_MLR, dbetadR, d2betadR2
  real(8)              :: yp, dypdR, d2ypdR2

  call get_y(Re, p, R, yp, dypdR, d2ypdR2)

  call get_uLR(Re, uLRe, duLRdR, d2uLRdR2)
  beta_inf = dlog(2.0*De/uLRe)

  call get_uLR(R, uLR, duLRdR, d2uLRdR2)

  call get_beta(R, beta_MLR, dbetadR, d2betadR2)

  E = -beta_MLR*yp
  dEdR = -dbetadR*yp-beta_MLR*dypdR
  d2EdR2 = -d2betadR2*yp-2*dbetadR*dypdR-beta_MLR*d2ypdR2

  EE = dexp(E)
  dEEdR = EE*dEdR
  d2EEdR2 = dEEdR*dEdR+EE*d2EdR2

  A = 1.d0-uLR/uLRe*EE
  dAdR = -duLRdR/uLRe*EE-uLR/uLRe*dEEDR
  d2AdR2 = -d2uLRdR2/uLRe*EE-2*duLRdR/uLRe*dEEDR-uLR/uLRe*d2EEdR2

  V = De*A**2-De
  dVdR = 2*De*A*dAdR
  d2VdR2 = 2*De*(dAdR**2+A*d2AdR2)
  return 

contains

  subroutine get_uLR(R, uLR, duLRdR, d2uLRdR2) 
    implicit none
    real(8), intent(in)  :: R
    real(8), intent(out) :: uLR, duLRdR, d2uLRdR2
    integer              :: i
    real(8)              :: D, dDdR, d2DdR2
    uLR = 0.d0
    duLRdR = 0.d0
    d2uLRdR2 = 0.d0
    do i = 1, N_C
      call get_D(m(i), R, D, dDdR, d2DdR2)
      uLR = uLR+D*C(i)/(R**m(i))
      duLRdR = duLRdR+dDdR*C(i)/R**m(i)-m(i)*D*C(i)/R**(m(i)+1)
      d2uLRdR2 = d2uLRdR2+d2DdR2*C(i)/R**m(i) &
        -2*dDdR*m(i)*C(i)/R**(m(i)+1) &
        +D*m(i)*(m(i)+1)*C(i)/R**(m(i)+2)
    end do
    return
  end subroutine

  subroutine get_D(m, R, D, dDdR, D2DdR2) 
    !
    ! Generalised Douketis-Scoles-type damping functions
    !
    ! Ref: Le Roy, Haugen, Tao, and Li,
    !      Mol. Phys. 109, 1598 (2011)
    !      doi:10.1080/00268976.2010.527304
    !
    implicit none
    integer, intent(in)  :: m
    real(8), intent(in)  :: R
    real(8), intent(out) :: D, dDdR, d2DdR2
    integer              :: i
    real(8)              :: FF, dFFdR, d2FFdR2
    real(8)              :: F, dFdR, d2FdR2
    ! These three lines are taken from reference above
    real(8), parameter   :: bds(6) = (/4.99d0, 4.53d0, 3.95d0, 3.69d0, 3.30d0, 2.50d0/)
    real(8), parameter   :: cds(6) = (/0.34d0, 0.36d0, 0.39d0, 0.405d0, 0.423d0, 0.468d0/)
    real(8), parameter   :: sds(6) = (/2.0d0, 1.0d0, 0.0d0, -0.5d0, -1.0d0, -2.0d0/)
    real(8)              :: b = 0.0d0, c = 0.0d0
    ! real(8)              :: res = 0.0d0
    logical              :: finds = .true.
    save b, c, finds

    if(finds) then
      do i = 1, 6
        if (abs(sds(i) - s) .lt. 1.0e-6) then
          b = bds(i)
          c = cds(i)
          finds = .false.
          exit
        end if
      end do
    endif

    F = -b*rho*r/(m*1.d0)-c*(rho*r)*(rho*r)/dsqrt(m*1.d0)
    dFdR = -b*rho/(m*1.d0)-2.d0*c*rho**2*r/dsqrt(m*1.d0)
    d2FdR2 = -2.d0*c*rho**2/dsqrt(m*1.d0)

    FF = dexp(F)
    dFFdR = FF*dFdR
    d2FFdR2 = dFFdR*dFdR+FF*d2FdR2

    D = (1.d0-FF)**(m*1.d0+s)
    dDdR = -(m*1.d0+s)*(1.d0-FF)**(m*1.d0+s-1.d0)*dFFdR
    d2DdR2 = -(m*1.d0+s) &
      *((m*1.d0+s-1.d0)*(1.d0-FF)**(m*1.d0+s-2.d0)*(-dFFdR**2) &
      +(1.d0-FF)**(m*1.d0+s-1.d0)*d2FFdR2)
    return 
  end subroutine

  subroutine get_y(Rref, q, R, y, dydR, d2ydR2) 
    real(8), intent(in)  :: Rref, R
    integer, intent(in)  :: q
    real(8), intent(out) :: y, dydR, d2ydR2
    real(8)              :: plus

    plus = (R**q+Rref**q)

    y = (R**q-Rref**q)/plus
    dydR = 2.d0*q*Rref**q*R**(q-1)/plus**2
    d2ydR2 = 2.d0*q*Rref**q/plus**4 &
      *(plus**2*(q-1)*R**(q-2)-plus*2.d0*q*R**(2*q-2))
    return 
  end subroutine

  subroutine get_beta(R, beta_MLR, dbetadR, d2betadR2) 
    implicit none
    real(8), intent(in)  :: R
    ! real(8) :: RR
    real(8), intent(out) :: beta_MLR, dbetadR, d2betadR2
    integer              :: i
    real(8)              :: S, dSdR, d2SdR2
    real(8)              :: yp, dypdR, d2ypdR2
    real(8)              :: yq, dyqdR, d2yqdR2

    call get_y(Rref, p, R, yp, dypdR, d2ypdR2)
    call get_y(Rref, q, R, yq, dyqdR, d2yqdR2)

    S = 0.d0
    dSdR = 0.d0
    d2SdR2 = 0.d0
    do i = 0, beta_max
      S = S+beta(i)*yq**i;
      if (i.gt.0) then
        dSdR = dSdR+beta(i)*i*yq**(i-1)*dyqdR;
      endif
      if (i.eq.1) then
        d2SdR2 = d2SdR2+beta(i)*i* &
          (yq**(i-1)*d2yqdR2)
      elseif (i.gt.1) then
        d2SdR2 = d2SdR2+beta(i)*i* &
          ((i-1)*yq**(i-2)*dyqdR**2+yq**(i-1)*d2yqdR2)
      endif
    end do
    beta_MLR = beta_inf*yp+(1.d0-yp)*S
    dbetadR = beta_inf*dypdR-dypdR*S+(1.d0-yp)*dSdR
    d2betadR2 = beta_inf*d2ypdR2-d2ypdR2*S &
      -2.d0*dypdR*dSdR+(1.d0-yp)*d2SdR2
    return 
  end subroutine
endsubroutine

