module lassofun
  implicit none

  private
  public :: lasso

contains

  function lasso(x,y,t,w0,maxiter,threshold,verbose) result(w)
    ! This is a Fortran re-implementation of the active set method
    ! coded by Mark Schmidt (UBC) in octave. Original source and notes:
    !
    ! https://www.cs.ubc.ca/~schmidtm/Software/lasso.html
    !
    ! Mark Schmidt. Graphical Model Structure Learning with L1-Regularization. Ph.D. Thesis, University of British Columbia, 2010.
    !
    real*8, intent(in) :: x(:,:)
    real*8, intent(in) :: y(:)
    real*8, intent(in) :: t
    real*8, intent(in), optional :: w0(:)
    integer, intent(in), optional :: maxiter
    real*8, intent(in), optional :: threshold
    logical, intent(in), optional :: verbose
    real*8 :: w(size(x,2))

    integer :: maxiter_
    real*8 :: threshold_
    integer :: i, j, k, n, p, s(1)
    integer :: iteration, qrpos(1), nsig
    integer, allocatable :: isigma(:)
    logical, allocatable :: sigma(:), sigma_old(:)
    real*8, allocatable :: xx(:,:), xy(:), xy_sigma(:), x_sigma(:,:), g(:)
    real*8, allocatable :: w_sigma(:), q(:,:), r(:,:), w_t(:), h(:), v_t(:)
    real*8, allocatable :: slope(:)
    integer*1, allocatable :: theta(:), theta_sigma(:)
    real*8 :: sum, v_denom, gamma ,min_gamma
    integer :: min_gamma_i
    logical :: verbose_

    ! process the optional parameters
    maxiter_ = 10000
    if (present(maxiter)) maxiter_ = maxiter
    threshold_ = 1d-9
    if (present(threshold)) threshold_ = threshold
    verbose_ = .false.
    if (present(verbose)) verbose_ = verbose
    if (present(w0)) then
       w = w0
    else
       w = 0d0
    end if

    ! initialize
    n = size(x,1)
    p = size(x,2)
    xx = matmul(transpose(x),x)
    xy = matmul(transpose(x),y)
    allocate(sigma(p),sigma_old(p))
    sigma_old = .true.
    sigma = (abs(w) > threshold_)
    allocate(theta(p))
    theta = xsign(w)
    allocate(q(n,n),r(n,p))
    iteration = 1

    if (verbose_) &
       write (*,'("iter   n(w)   n(step)   f(w)   opt(wi)    free")')
    ! main loop
    do while (iteration < maxiter_)
       ! get the values associated with the active set
       isigma = pack((/(i,i=1,size(sigma,1))/),sigma)
       w_sigma = w(isigma)
       theta_sigma = theta(isigma)
       ! if ((iteration > 1 .or. warm) .and. any(sigma .neqv. sigma_old)) then
       if (any(sigma .neqv. sigma_old)) then
          xy_sigma = xy(isigma)

          ! QR insert
          if (count(sigma) <= 1 .or. iteration == 1) then
             x_sigma = x(:,isigma)
             call qr(x_sigma,q,r,n,size(isigma,1))
          else
             qrpos =  findloc(sigma(isigma) .neqv. sigma_old(isigma),.true.)
             nsig = size(isigma,1) - 1
             call qrinsert(q,r,qrpos(1),x(:,s(1)),'c')
          end if
       end if

       ! Solve the local linerization
       call solveKKT2()

       ! sign infeasible case
       do while (count(xsign(w_t(isigma)) == theta_sigma) /= count(sigma))
          if (verbose_) &
             write (*,'("-> not sign feasible")')
          ! A1: Find first zero-crossing
          min_gamma = 1d0
          min_gamma_i = -1
          do k = 1, size(w,1)
             if (abs(w(k)) > threshold_) then
                gamma = -w(k)/h(k)
                if (gamma > 0 .and. gamma < min_gamma) then
                   min_gamma = gamma
                   min_gamma_i = k
                end if
             end if
          end do
          if (min_gamma_i == -1) then
             if (verbose_) &
                write (*,'("-> numerical breakdown, check for dependent columns")')
             exit ! Numerical breakdown, check for dependent columns
          end if

          ! A1: set beta to h truncated at first zero-crossing
          w = w + min_gamma * h
          ! A2: reverse sign of first zero-crossing
          theta(min_gamma_i) = -theta(min_gamma_i)
          ! A2: recompute h
          theta_sigma = theta(isigma)
          w_sigma = w(isigma)
          call solveKKT2()

          if (count(xsign(w_t(isigma)) == theta_sigma) /= count(sigma)) then
             if (verbose_) &
                write (*,'("-> still sign infeasible")')
             ! still sign infeasible
             ! A3: Remove the first zero-crossing from the active set
             sigma_old = sigma
             sigma(min_gamma_i) = .false.
             isigma = pack((/(i,i=1,size(sigma,1))/),sigma)
             w(min_gamma_i) = 0d0

             ! A3: Reset beta(min_gamma) and theta(min_gamma)
             w(min_gamma_i) = 0d0
             theta(min_gamma_i) = 0
             w_sigma = w(isigma)
             theta_sigma = theta(isigma)

             ! A3: Recompute h
             xy_sigma = xy(isigma)

             ! QR Update
             x_sigma = x(:,isigma)
             call qr(x_sigma,q,r,n,size(isigma,1)) ! could do a qr delete here instead
             call solveKKT2()

             ! final print
             if (verbose_) then
                if (count(xsign(w_t(isigma)) == theta_sigma) /= count(sigma)) then
                   write (*,'("-> it is still not sign feasible: going another round")')
                else
                   write (*,'("-> finally, it is sign feasible")')
                end if
             end if
          else
             if (verbose_) &
                write (*,'("-> now it is sign feasible")')
          end if
       end do

       ! compute violation
       v_denom = maxval(abs(xy_sigma - matmul(matmul(transpose(x(:,isigma)),x),w_t)))
       if (v_denom > 0) then
          v_t = (xy - matmul(xx,w_t)) / v_denom
       else
          if (.not.allocated(v_t)) allocate(v_t(p))
          v_t = sign(huge(1d0),xy - matmul(xx,w_t))
       end if
       j = p - count((abs(w_t) < threshold_) .and. (abs(v_t) > 1+threshold_))

       if (verbose_) then
          write (*,'(I6,2X,3(E12.4,1X),1X,2(I5,X))') iteration, sum(abs(w_t)), sum(abs(w_t-w)), &
             sum((matmul(x,w) - y)**2), j, count(sigma)
       end if

       ! check for optimality
       if (j == p .or. t == 0) exit

       ! find and add the most violating variable
       ! On the first iteration, all variables are equally violating
       ! so we're going to use the Shevade/Perkins trick to introduce
       ! a good first variable, this often means we may not have to deal with
       ! sign feasibility later issues later
       if (iteration == 1) then
          g = matmul(xx,w) - xy
          if (.not.allocated(slope)) allocate(slope(p))
          slope = 0d0
          do i = 1, p
             if (g(i) > t) then
                slope(i) = g(i) + t
             elseif (g(i) < -t) then
                slope(i) = g(i) - t
             elseif (abs(w(i)) > threshold_) then
                slope(i) = g(i) + t * ssign(w(i))
             end if
          end do
          s = maxloc(abs(slope))
       else
          s = maxloc(abs(v_t),.not.sigma)
       end if

       sigma_old = sigma
       sigma(s(1)) = .true.
       theta(s(1)) = ssign(v_t(s(1)))
       w = w_t
       iteration = iteration + 1
    end do

  contains
    subroutine solveKKT2()
      ! X = Q*R
      ! mu = max(0, theta'(X'X)^-1X'y - t) / theta'(X'X)^-1 theta)
      ! h = (X'X)^-1 * (X'(Y-Xw) - mu*theta)
      real*8 :: mu, mu_denom
      real*8, allocatable :: b(:), rred(:,:), qred(:,:)
      ! integer :: nsig

      if (.not.allocated(w_t)) then
         allocate(w_t(p))
         w_t = 0d0
      end if
      if (.not.allocated(h)) then
         allocate(h(p))
         h = 0d0
      end if
      ! if (iteration == 1) return

      ! reduce r
      nsig = size(isigma,1)
      rred = r(1:nsig,1:nsig)
      qred = q(:,1:nsig)

      ! calculate mu_denom
      b = real(theta_sigma,8)
      call solve_triangular(rred,b,'u','t')
      call solve_triangular(rred,b,'u','n')
      mu_denom = dot_product(real(theta_sigma,8),b)

      ! calcualte mu
      if (mu_denom > 0d0) then
         b = matmul(transpose(qred),y)
         call solve_triangular(rred,b,'u','n')
         mu = max((dot_product(theta_sigma,b) - t) / mu_denom,0d0)
      else
         mu = 0d0
      end if

      ! calculate h
      h = 0d0
      b = matmul(transpose(rred),matmul(transpose(qred),y) - matmul(rred,w_sigma)) - mu * theta_sigma
      call solve_triangular(rred,b,'u','t')
      call solve_triangular(rred,b,'u','n')
      h(isigma) = b
      w_t = w + h

    end subroutine solveKKT2

    function xsign(w) result(sw)
      real*8, intent(in) :: w(:)
      integer :: sw(size(w,1))

      where(w == 0)
         sw = 0
      elsewhere
         sw = int(sign(1d0,w))
      end where

    end function xsign

    function ssign(w) result(sw)
      real*8, intent(in) :: w
      integer :: sw

      if (w == 0) then
         sw = 0
      elseif (w > 0) then
         sw = 1
      else
         sw = -1
      end if

    end function ssign

  end function lasso

  subroutine qr(x,q,r,rows,cols)
    real*8, intent(inout) :: x(:,:)
    real*8, intent(inout), allocatable :: q(:,:), r(:,:)
    integer, intent(in) :: rows, cols

    integer :: lwork, info, i
    real*8, allocatable :: tau(:), work(:)
    real*8 :: xwork(1)
    real*8, allocatable :: rred(:,:), qred(:,:)

    ! initialize
    allocate(tau(min(rows,cols)))

    ! run the QR factorization
    call dgeqrf(rows,cols,x,size(x,1),tau,xwork,-1,info)
    if (info /= 0) call error()
    lwork = int(xwork(1))
    if (allocated(work)) deallocate(work)
    allocate(work(lwork))
    call dgeqrf(rows,cols,x,size(x,1),tau,work,lwork,info)
    if (info /= 0) call error()

    ! assign r
    if (allocated(r)) deallocate(r)
    allocate(r(cols,cols))
    r = 0d0
    do i = 1, cols
       r(1:min(rows,i),i) = x(1:min(rows,i),i)
    end do

    ! unpack q
    if (allocated(q)) deallocate(q)
    allocate(q(rows,cols))
    q = 0d0
    do i = 1, min(rows,cols)
       q(i,i) = 1d0
    end do
    call dormqr('L','N',rows,cols,min(rows,cols),x,size(x,1),tau,q,size(q,1),xwork,-1,info)
    if (info /= 0) call error()
    lwork = int(xwork(1))+1
    if (allocated(work)) deallocate(work)
    allocate(work(lwork))
    call dormqr('L','N',rows,cols,min(rows,cols),x,size(x,1),tau,q,size(q,1),work,lwork,info)
    if (info /= 0) call error()

  contains
    subroutine error()
      write (*,'("error in qr: ",I4)') info
      stop 1
    end subroutine error
  end subroutine qr

  subroutine qrinsert(q,r,qrpos,x,rc)
    real*8, intent(inout), allocatable :: q(:,:), r(:,:)
    integer, intent(in) :: qrpos
    real*8, intent(in) :: x(:)
    character*1, intent(in) :: rc

    real*8, allocatable :: work(:), temp(:,:)
    integer :: m, n, k

    if (rc == 'c') then
       ! update column
       ! A is mxn, Q is a square matrix mxm, R is mxn

       ! create a new row and column in r
       allocate(temp(size(r,1)+1,size(r,2)+1))
       temp = 0d0
       temp(1:size(r,1),1:size(r,2)) = r
       call move_alloc(temp,r)

       ! create a new column in q
       allocate(temp(size(q,1),size(q,2)+1))
       temp(:,1:size(q,2)) = q
       temp(:,size(q,2)+1) = 0d0
       call move_alloc(temp,q)

       ! insert the column
       m = size(q,1)
       n = size(r,2)-1
       k = size(q,2)-1
       allocate(work(k))
       call dqrinc(m,n,k,q,size(q,1),r,size(r,1),qrpos,x,work)
    else
       write (*,'("not implemented!")')
       stop 1
    end if

  end subroutine qrinsert

  subroutine solve_triangular(a,b,uplo,trans)
    ! solves a * x = b, where a is a triangular matrix
    ! uplo -> 'u' is upper triangular, 'l' is lower triangular
    ! trans -> 'n' is Ax = b, 't' or 'c' is A^Tx=b
    real*8, intent(in) :: a(:,:)
    real*8, intent(inout) :: b(:)
    character*1, intent(in) :: uplo, trans

    integer :: n

    n = size(a,1)
    if (n == 0) return
    if (size(a,2) /= n .or. size(b,1) /= n) then
       write (*,*) size(a,1), size(a,2), size(b,1)
       write (*,'("incorrect sizes in solve_triangular")')
       stop 1
    end if

    call dtrsv(uplo,trans,'n',n,a,n,b,1)

  end subroutine solve_triangular

end module lassofun
