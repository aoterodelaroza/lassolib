module wrappers
  implicit none

  private
  public :: lasso_for_acps

contains

  subroutine lasso_c(rows,cols,x,y,t,w) bind(c,name="lasso_c")
    use iso_c_binding, only: c_int, c_double
    integer(c_int), intent(in), value :: rows, cols
    real(c_double), intent(in) :: x(rows,cols), y(rows)
    real(c_double), intent(in), value :: t
    real(c_double), intent(inout) :: w(cols)

    w = lasso_for_acps(rows,cols,x,y,t)

  end subroutine lasso_c

  function lasso_for_acps(rows,cols,x,y,t,maxcoef,yadd,verbose) result(w)
    use lassofun, only: lasso
    integer, intent(in) :: rows, cols
    real*8, intent(in) :: x(rows,cols), y(rows)
    real*8, intent(in) :: t
    real*8, intent(in), optional :: maxcoef(cols)
    real*8, intent(in), optional :: yadd(rows,*)
    logical, intent(in), optional :: verbose
    real*8 :: w(cols)

    real*8, allocatable :: factor(:), xtilde(:,:)
    logical, allocatable :: idx(:)
    logical :: verbose_
    real*8 :: ratio, wrms
    integer :: i

    verbose_ = .false.
    if (present(verbose)) verbose_ = verbose

    if (present(yadd)) then
       write (*,*) "not implemented yet!"
       stop 1
    end if
    if (present(maxcoef)) then
       ! initialize
       allocate(factor(cols),xtilde(rows,cols))
       factor = 1d0
       xtilde = 0d0

       do while (.true.)
          do i = 1, cols
             xtilde(:,i) = x(:,i) * factor(i)
          end do
          w = lasso(xtilde,y,t)
          do i = 1, cols
             w(i) = w(i) * factor(i)
          end do
          wrms = sqrt(sum((y - matmul(x,w))**2))
          idx = (abs(w) > maxcoef)
          if (verbose_) then
             write (*,'("norm-1 = ",F12.4,", wrms = ",F12.4,", maxcoef: ok = ",I4," / violations = ",I4)') &
                sum(abs(w)),wrms,cols-count(idx),count(idx)
          end if
          if (count(idx) == 0) exit

          do i = 1, cols
             if (abs(w(i)) > maxcoef(i)) then
                ratio = min(maxcoef(i) / abs(w(i)),0.99)
                factor(i) = factor(i) * ratio**(1.1d0)
             end if
          end do
       end do
       deallocate(xtilde,factor)
    else
       w = lasso(x,y,t)
    end if

  end function lasso_for_acps

end module wrappers
