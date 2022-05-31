subroutine pot(r,v,dv)
  implicit none
  double precision, intent(in) :: r
  double precision, intent(out) :: v, dv
  double precision :: d2v
  double precision :: r_6
  r_6=r**(-6)
  v = (4 * (r_6 * r_6 - r_6))
  dv = (4 * ((-12) * r_6 * r_6 / r - (-6) * r_6 / r))
  return
end subroutine

