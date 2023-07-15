subroutine lasso_c(rows,cols,x,y,t,w,wrms) bind(c,name="lasso_c")
  use iso_c_binding, only: c_int, c_double
  use wrappers, only: lasso_for_acps
  integer(c_int), intent(in), value :: rows, cols
  real(c_double), intent(in) :: x(rows,cols), y(rows)
  real(c_double), intent(in), value :: t
  real(c_double), intent(inout) :: w(cols)
  real(c_double), intent(inout) :: wrms

  w = lasso_for_acps(rows,cols,x,y,t)
  wrms = sqrt(sum((y - matmul(x,w))**2))

end subroutine lasso_c

