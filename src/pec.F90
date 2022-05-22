
subroutine getfeatures(pot, rmin, eps, sgm)
  implicit none
  double precision :: R1,R2,V2,sgm,Rmin,eps,dVdR,Vm,Rm,inc
  external pot

  R1 = 1.0d0
  R2 = 1.0d0
  call pot(R2, V2, dVdR)
  do while (dabs(dvdr).gt.1.0d-10)
    inc=-dvdr*0.001d0
    ! In case we go too fast... 0.5 ang is safe in most cases.
    if (dabs(inc).gt.0.5d0) inc = inc/dabs(inc)*0.5d0
    R2=R2+inc
    call pot(R2, V2, dVdR)
  enddo
  Rmin=R2
  eps=-V2
  do while (R2-R1.gt.1.0d-10)
    Rm=(R1+R2)*0.5d0
    call pot(Rm, Vm, dVdR)
    if (Vm.gt.0.0d0) R1=Rm
    if (Vm.le.0.0d0) R2=Rm
    ! write(*,999) R1, R2, R2-R1, Rm
  enddo
  sgm=Rm
end subroutine

!***********************************************************************
double precision function vf(r)
  implicit none
  double precision, intent(in) :: r
  double precision :: v, dv, rr, rmin, eps, sgm
  common /feat/ rmin, eps, sgm
  rr=r*sgm
  call pot(rr,v,dv)
  vf=v/eps
  return
end function

double precision function vd(r)
  implicit none
  double precision, intent(in) :: r
  double precision :: v, dv, rr, rmin, eps, sgm
  common /feat/ rmin, eps, sgm
  rr=r*sgm
  call pot(rr,v,dv)
  vd=dv/eps*sgm
  return
end function

