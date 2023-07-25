subroutine lasso_c(rows,cols,x,y,t,w,wrms,maxcoefp) bind(c,name="lasso_c")
  use iso_c_binding, only: c_int, c_double, c_ptr, c_associated, c_f_pointer
  use wrappers, only: lasso_for_acps
  integer(c_int), intent(in), value :: rows, cols
  real(c_double), intent(in) :: x(rows,cols), y(rows)
  real(c_double), intent(in), value :: t
  real(c_double), intent(inout) :: w(cols)
  real(c_double), intent(inout) :: wrms
  type(c_ptr), intent(in), value :: maxcoefp

  real(c_double), pointer :: maxcoef(:)

  if (c_associated(maxcoefp)) then
     call c_f_pointer(maxcoefp,maxcoef,shape=(/cols/))
     w = lasso_for_acps(rows,cols,x,y,t,maxcoef=maxcoef,verbose=.true.)
  else
     w = lasso_for_acps(rows,cols,x,y,t)
  end if
  wrms = sqrt(sum((y - matmul(x,w))**2))

end subroutine lasso_c

subroutine lasso_cv_c(rows,cols,x,y,t,w,wrms) bind(c,name="lasso_cv_c")
  use iso_c_binding, only: c_int, c_double
  use wrappers, only: lasso_cv
  integer(c_int), intent(in), value :: rows, cols
  real(c_double), intent(in) :: x(rows,cols), y(rows)
  real(c_double), intent(inout) :: t
  real(c_double), intent(inout) :: w(cols)
  real(c_double), intent(inout) :: wrms

  call lasso_cv(rows,cols,x,y,t,w,.true.)
  wrms = sqrt(sum((y - matmul(x,w))**2))

end subroutine lasso_cv_c
